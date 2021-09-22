/************************************************************************/
/**

   Program:    absplit
   \file       absplit.c
   
   \version    V0.1
   \date       17.09.21   
   \brief      Split an antibody PDB file into Fvs with antigens
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2021
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is distributed under the conditions of the 
   GNU GENERAL PUBLIC LICENSE (GPL) Version 3

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#ifdef __linux__
#  include <linux/limits.h>
#else
#  define PATH_MAX 4096
#endif

#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/sequtil.h"

#include "absplit.h"

/************************************************************************/
/* Defines and macros
*/
#define PROGNAME        "absplit"
#define VERSION         0.1
#define MAXBUFF         240
#define HUGEBUFF        10000
#define MAXSEQ          10000
#define ABTHRESHOLD     0.5
#define GAPOPENPENALTY  5
#define GAPEXTPENALTY   2
#define SCOREMATRIX     "BLOSUM62"
#define MAXINTERFACE    20
#define MAXCDRRES       160
#define MAXRESID        16
#define COFGDISTCUTSQ   1225.0 /* 35^2 - used to find possible VH/VL pairs */
#define INTDISTCUTSQ    900.0  /* 30^2 - used to find VH/VL interface contact */
#define CONTACTDISTSQ   36.0   /* 6^2  - used to find antigen contacts */
#define MINAGCONTACTS   14     /* Tweaked with CONTACTDISTSQ to get Ag for
                                  6o8d but not too much for 1dee */
#define MAXANTIGEN      16
#define MAXCHAINLABEL   8
#define CHAINTYPE_ATOM  (APTR)1
#define CHAINTYPE_HET   (APTR)2

typedef struct _domain
{
   int   domainNumber,
         startSeqRes,
         lastSeqRes,
         nInterface,
         nAntigenChains,
         interface[MAXINTERFACE],
         nCDRRes,
         CDRRes[MAXCDRRES];
   char  domSeq[MAXSEQ],
         chainType;
   PDB   *startRes,
         *lastRes,
         *stopRes;
   VEC3F CofG,
         IntCofG;
   REAL  pairIntDistSq,
         pairCofGDistSq;
   BOOL  used;
   
   PDBCHAIN *chain,
            *antigenChains[MAXANTIGEN];
   struct _domain *pairedDomain;
   struct _domain *next;
}
DOMAIN;

/************************************************************************/
/* Globals
*/
BOOL gVerbose = TRUE;
BOOL gNoAntigen = FALSE;

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile);
void UsageDie(void);
BOOL ProcessFile(WHOLEPDB *wpdb, char *infile, FILE *dataFp);
DOMAIN *FindVHVLDomains(PDBCHAIN *chain, FILE *dataFp, DOMAIN *domains);
void GetSequenceForChain(PDBCHAIN *chain, char *sequence);
void ExePathName(char *str, BOOL pathonly);
BOOL CheckAndMask(char *sequence, FILE *dataFp, PDBCHAIN *chain, DOMAIN **pDomains);
FILE *OpenSequenceDataFile(void);
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2);
void MaskAndAssignDomain(char *seq, PDBCHAIN *chain, char *bestMatch, char *aln1, char *aln2, DOMAIN **pDomains);
void SetChainType(DOMAIN *domain, char *header);
void SetIFResidues(DOMAIN *domain, char *header);
void SetCDRResidues(DOMAIN *domain, char *header);
void PrintDomains(DOMAIN *domains);
void SetDomainBoundaries(DOMAIN *domain);
void PairDomains(DOMAIN *domains);
void WriteDomains(DOMAIN *domains, char *filestem);
void FlagAntigens(DOMAIN *domains, PDBSTRUCT *pdbs);
void CheckAntigenContacts(DOMAIN *domain, PDBSTRUCT *pdbs);
void SetChainTypes(PDBCHAIN *chains);
BOOL inIntArray(int value, int *array, int arrayLen);
BOOL regionsMakeContact(PDB *start1, PDB *stop1, PDB *start2, PDB *stop2);





/************************************************************************/
int main(int argc, char **argv)
{
   char infile[MAXBUFF];
   
   if(ParseCmdLine(argc, argv, infile))
   {
      FILE *fp = NULL;
      
      if((fp = fopen(infile, "r"))!=NULL)
      {
         WHOLEPDB *wpdb = NULL;
         if((wpdb = blReadWholePDB(fp))!=NULL)
         {
            FILE *dataFp;
            if((dataFp=OpenSequenceDataFile())==NULL)
            {
               fprintf(stderr,"Error (%s): The antibody sequence \
datafile was not installed\n", PROGNAME);
               exit(1);
            }
            

            blReadMDM(SCOREMATRIX);
            
            if(!ProcessFile(wpdb, infile, dataFp))
            {
               fprintf(stderr,"Error (%s): Unable to split PDB into chains\n",
                       PROGNAME);
               exit(1);
            }
            
            blFreeWholePDB(wpdb);
         }
         else
         {
            fprintf(stderr,"Error (%s): Nothing read from input file (%s)\n",
                    PROGNAME, infile);
            exit(1);
         }

         fclose(fp);
      }
      else
      {
         fprintf(stderr,"Error (%s): Can't read input file (%s)\n",
                 PROGNAME, infile);
         exit(1);
      }
   }
   else
   {
      UsageDie();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \return                     Success

   Parse the command line

-  17.09.21 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile)
{
   argc--;
   argv++;
   
   infile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            gVerbose = TRUE;
            break;
         case 'q':
            gVerbose = FALSE;
            break;
         case 'n':
            gNoAntigen = TRUE;
            break;
         case 'h':
            return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there is only 1 argument left                    */
         if(argc != 1)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(FALSE);
}

/************************************************************************/
void UsageDie(void)
{
   exit(1);
}

/************************************************************************/
void GetFilestem(char *infile, char *filestem)
{
   char *start, *stop;
   
   start  = strrchr(infile, '/');
   start  = (start==NULL)?infile:start+1;
   stop   = strrchr(infile, '.');
   if(stop < start)
   {
      strcpy(filestem, start);
   }
   else
   {
      int  length;
      length = stop-start;
      strncpy(filestem, start, (int)length);
      filestem[(int)length] = '\0';
   }
}

/************************************************************************/
BOOL ProcessFile(WHOLEPDB *wpdb, char *infile, FILE *dataFp)
{
   char      filestem[MAXBUFF];
   PDBSTRUCT *pdbs;
   DOMAIN    *domains = NULL;

   GetFilestem(infile, filestem);
   if((pdbs = blAllocPDBStructure(wpdb->pdb))!=NULL)
   {
      PDBCHAIN *chain;

      SetChainTypes(pdbs->chains);
      
      for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
      {
#ifdef DEBUG
         printf("Ptr: %ld Type: %s\n", (ULONG)chain->start,
                chain->start->record_type);
#endif
         
/*         if(!strncmp(chain->start->record_type, "ATOM  ", 6)) */
         if(chain->extras == CHAINTYPE_ATOM)
         {
            printf("\n***Handling chain: %s\n", chain->chain);
            domains = FindVHVLDomains(chain, dataFp, domains);
         }
      }

      PairDomains(domains);
      FlagAntigens(domains, pdbs);
      PrintDomains(domains);
      WriteDomains(domains, filestem);
      
      FREELIST(domains, DOMAIN);
      blFreePDBStructure(pdbs);
   }
   else
   {
      return(FALSE);
   }

   return(TRUE);
}


void GetSequenceForChain(PDBCHAIN *chain, char *sequence)
{
   int i=0;
   PDBRESIDUE *r;
   if(chain->extras == CHAINTYPE_ATOM)
   {
      for(r=chain->residues; r!=NULL; NEXT(r))
      {
/*      if(!strncmp(r->start->record_type, "ATOM  ", 6)) */
         {
            sequence[i++] = blThrone(r->resnam);
         }
      }
   }
   
   sequence[i] = '\0';
}

/* From https://www.linuxquestions.org/questions/showthread.php?threadid=273614
 */
void ExePathName(char *str, BOOL pathonly)
{
  FILE *fp;
  char buf[PATH_MAX+100], *p;

  *str = '\0';
  if(!(fp = fopen("/proc/self/maps", "r")))
    return;

  fgets(buf, sizeof(buf), fp);
  fclose(fp);

  *(p = strchr(buf, '\n')) = '\0';
  while(*p != ' ')
    p--;

  strncpy(str, p+1, PATH_MAX);
  str[PATH_MAX-1] = '\0';
  if(pathonly)
  {
     p=strrchr(str, '/');
     if(p!=NULL) *p='\0';
  }
}


FILE *OpenSequenceDataFile(void)
{
   char pathname[PATH_MAX];
   
   ExePathName(pathname, TRUE);
   strcat(pathname,ABSEQFILE);
   
   return(fopen(pathname, "r"));
}


/************************************************************************/
/*>REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2)
   ---------------------------------------------------------------------
*//**
   \param[in]   theSeq   the sequence of interest
   \param[in]   seq      the database sequence
   \param[out]  align1   Alignment of our sequence
   \param[out]  align2   Alignment of database sequence
   \return               Score for alignment

   - 31.03.20 Original   By: ACRM
   - 20.09.21 Modified to use affine alignment and mutation matrix
*/
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2)
{
   int  score;
   int  alignLen;
   int  bestPossibleScore;

   bestPossibleScore = blAffinealign(seq, strlen(seq),
                                     seq, strlen(seq),
                                     FALSE,          /* verbose  */
                                     FALSE,          /* identity */
                                     GAPOPENPENALTY, /* penalty  */
                                     GAPEXTPENALTY,  /* extension  */
                                     align1,
                                     align2,
                                     &alignLen);

   score = blAffinealign(theSeq, strlen(theSeq),
                         seq, strlen(seq),
                         FALSE,          /* verbose  */
                         FALSE,          /* identity */
                         GAPOPENPENALTY, /* penalty  */
                         GAPEXTPENALTY,  /* extension  */
                         align1,
                         align2,
                         &alignLen);
   align1[alignLen] = align2[alignLen] = '\0';

#ifdef DEBUG   
   fprintf(stderr, "\n>>>%s\n", align1);
   fprintf(stderr, ">>>%s\n",   align2);
#endif
   
   return((REAL)score / (REAL)bestPossibleScore);
}


DOMAIN *FindVHVLDomains(PDBCHAIN *chain, FILE *dataFp, DOMAIN *domains)
{
   char sequence[MAXSEQ];
   
   GetSequenceForChain(chain, sequence);
#ifdef DEBUG
   printf("Chain: %s Sequence: %s\n", chain->chain, sequence);
#endif
   while(TRUE)
   {
      if(!CheckAndMask(sequence, dataFp, chain, &domains)) break;
   }
   return(domains);
   
}

BOOL CheckAndMask(char *sequence, FILE *dbFp, PDBCHAIN *chain, DOMAIN **pDomains)
{
   char        bestMatch[MAXBUFF+1];
   REAL        maxScore = 0.0;
   static char align1[HUGEBUFF+1],
               align2[HUGEBUFF+1];
   static char bestAlign1[HUGEBUFF+1],
               bestAlign2[HUGEBUFF+1];
   char header[MAXBUFF+1];
   char *seq = NULL;
   BOOL found = FALSE;
   rewind(dbFp);

   /* Find the best match in the reference sequences */
   while((seq = blReadFASTA(dbFp, header, MAXBUFF))!=NULL)
   {
      REAL score;
         
      score = CompareSeqs(sequence, seq, align1, align2);
            
      if(score > maxScore)
      {
         maxScore = score;
         strncpy(bestMatch,  header, MAXBUFF);
         strncpy(bestAlign1, align1, HUGEBUFF);
         strncpy(bestAlign2, align2, HUGEBUFF);
      }
   }

   /* If we found an antibody sequence               */
   if(maxScore > ABTHRESHOLD)
   {
      if(gVerbose)
      {
         fprintf(stderr, "Best match: %s Score: %.4f\n",
                 bestMatch, maxScore);
         fprintf(stderr, "SEQ: %s\n",   bestAlign1);
         fprintf(stderr, "REF: %s\n\n", bestAlign2);
      }
      
      MaskAndAssignDomain(sequence, chain, bestMatch, bestAlign1, bestAlign2, pDomains);
      found = TRUE;
   }
   
   return(found);
}

void SetChainType(DOMAIN *domain, char *header)
{
   char *ptr;
   domain->chainType = '?';
   if((ptr = strchr(header, '|'))!=NULL)
   {
      ptr--;
      domain->chainType = *ptr;
   }
}

void SetIFResidues(DOMAIN *domain, char *header)
{
   char *ptr1, *ptr2, *bar;
   char headerCopy[MAXBUFF+1];
   
   domain->nInterface = 0;
   strncpy(headerCopy, header, MAXBUFF);
   headerCopy[MAXBUFF] = '\0';
   if((ptr1 = strchr(headerCopy, '['))!=NULL)
   {
      ptr1++;
      bar = strchr(ptr1, '|');
      
      while((ptr2 = strchr(ptr1, ','))!=NULL)
      {
         if(ptr2 > bar)
            break;

         *ptr2 = '\0';
         sscanf(ptr1, "%d", &(domain->interface[domain->nInterface++]));
         ptr1=ptr2+1;
      }
      
      if((ptr2 = strchr(ptr1, ']'))!=NULL)
      {
         *ptr2 = '\0';
         sscanf(ptr1, "%d", &(domain->interface[domain->nInterface++]));
      }
   }
}

void SetCDRResidues(DOMAIN *domain, char *header)
{
   char *ptr1, *ptr2;
   char headerCopy[MAXBUFF+1];
   
   domain->nCDRRes = 0;

   strncpy(headerCopy, header, MAXBUFF);
   headerCopy[MAXBUFF] = '\0';
   if((ptr1 = strchr(headerCopy, '['))!=NULL)
   {
      ptr1++;
      if((ptr1 = strchr(ptr1, '['))!=NULL)
      {
         ptr1++;
      
         while((ptr2 = strchr(ptr1, ','))!=NULL)
         {
            *ptr2 = '\0';
            sscanf(ptr1, "%d", &(domain->CDRRes[domain->nCDRRes++]));
            ptr1=ptr2+1;
         }
      
         if((ptr2 = strchr(ptr1, ']'))!=NULL)
         {
            *ptr2 = '\0';
            sscanf(ptr1, "%d", &(domain->CDRRes[domain->nCDRRes++]));
         }
      }
      
   }
}

void PrintDomains(DOMAIN *domains)
{
   DOMAIN *d;

   printf("\n***Results\n");
   
   for(d=domains; d!=NULL; NEXT(d))
   {
      printf("DomNum: %d Chain: %s Start: %d Stop: %d Type: %c PairsWith: %d\n",
             d->domainNumber,
             d->chain->chain,
             d->startSeqRes,
             d->lastSeqRes,
             d->chainType,
             (d->pairedDomain==NULL)?0:d->pairedDomain->domainNumber);
      printf("%s\n", d->domSeq);
      if(d->nAntigenChains)
      {
         int i;
         printf("   AntigenChains: ");
         for(i=0; i<d->nAntigenChains; i++)
         {
            printf("%s ",d->antigenChains[i]->chain);
         }
         printf("\n");
      }
      printf("\n");
#ifdef DEBUG
      {
         int i;
         printf("IF: ");
         for(i=0; i<d->nInterface; i++)
         {
            printf("%d ", d->interface[i]);
         }
         printf("\n\n");
         printf("CDR: ");
         for(i=0; i<d->nCDRRes; i++)
         {
            printf("%d ", d->CDRRes[i]);
         }
         printf("\n\n");
      }
#endif
/*
      blWritePDBRecord(stdout,d->startRes);
      blWritePDBRecord(stdout,d->lastRes);
*/
      
   }
}

void MaskAndAssignDomain(char *seq, PDBCHAIN *chain, char *header, char *seqAln,
                         char *domAln, DOMAIN **pDomains)
{
   int    seqPos    = 0,
          alnPos    = 0,
          domSeqPos = 0;
   DOMAIN *d, *prevD;

   if(*pDomains == NULL)
   {
      INIT((*pDomains), DOMAIN);
      d = *pDomains;
      prevD = NULL;
   }
   else
   {
      d = *pDomains;
      LAST(d);
      prevD = d;
      ALLOCNEXT(d, DOMAIN);
   }
   if(d==NULL)
   {
      fprintf(stderr,"Error (%s): No memory for list of domains\n", PROGNAME);
      exit(1);
   }
   d->startSeqRes = -1;
   d->lastSeqRes  = -1;
   d->nInterface  = 0;
   d->chain       = chain;
   d->pairIntDistSq = 100000000.0;
   d->pairCofGDistSq = 100000000.0;
   d->domainNumber = (prevD==NULL)?1:prevD->domainNumber+1;
   d->pairedDomain = NULL;
   
   SetChainType(d, header);
   SetIFResidues(d, header);
   SetCDRResidues(d, header);
   
   for(seqPos=0, alnPos=0; seqPos<strlen(seqAln); seqPos++)
   {
      while(seqAln[alnPos] == '-')  /* Skip insertions */
      {
         alnPos++;
      }
      if((domAln[alnPos] != '-') &&
         (domAln[alnPos] != '\0'))
      {
         if(d->startSeqRes == (-1))
            d->startSeqRes = seqPos;
#ifdef DEBUG
         printf("Seqpos %d SeqRes %c DomRes %c\n", seqPos, seqAln[alnPos], domAln[alnPos]);
#endif
         d->lastSeqRes = seqPos;
         d->domSeq[domSeqPos++] = seq[seqPos];
         seq[seqPos] = 'X';
      }
      alnPos++;
   }
   d->domSeq[domSeqPos] = '\0';

   SetDomainBoundaries(d);
}

void SetDomainBoundaries(DOMAIN *domain)
{
   int resnum = 0,
       nCoor  = 0,
       i;
   PDB *p,
       *nextRes;

   /* Given the integer sequence positions, find the pointers to
      the PDB residues for these positions
   */
   for(p=domain->chain->start; p!=domain->chain->stop; p=nextRes)
   {
      nextRes = blFindNextResidue(p);
      if(resnum == domain->startSeqRes)
      {
         domain->startRes = p;
      }
      else if(resnum == domain->lastSeqRes)
      {
         domain->lastRes = p;
         domain->stopRes = nextRes;
         break;
      }
      
      resnum++;
   }

   /* Find the CofG of the domain                                       */
   domain->CofG.x = domain->CofG.y = domain->CofG.z = 0.0;
   nCoor  = 0;
   
   for(p=domain->startRes; p!=domain->stopRes; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ", 4))
      {
         domain->CofG.x += p->x;
         domain->CofG.y += p->y;
         domain->CofG.z += p->z;
         nCoor++;
      }
   }
   domain->CofG.x /= nCoor;
   domain->CofG.y /= nCoor;
   domain->CofG.z /= nCoor;

   /* Find the CofG of the VH/VL interface residues                     */
   domain->IntCofG.x = domain->IntCofG.y = domain->IntCofG.z = 0.0;
   nCoor  = 0;
   resnum = 0;

   for(p=domain->startRes; p!=domain->stopRes; p=nextRes)
   {
      PDB *q;
      nextRes = blFindNextResidue(p);
      for(q=p; q!=nextRes; NEXT(q))
      {
         for(i=0; i<domain->nInterface; i++)
         {
            if(resnum == domain->interface[i])
            {
               PDB *r;
               
               for(r=p; r!=nextRes; NEXT(r))
               {
                  if(!strncmp(r->atnam, "CA  ", 4))
                  {
                     domain->IntCofG.x += r->x;
                     domain->IntCofG.y += r->y;
                     domain->IntCofG.z += r->z;
                     nCoor++;
                  }
               }
               break;
            }

         }
      }
      resnum++;
   }
   domain->IntCofG.x /= nCoor;
   domain->IntCofG.y /= nCoor;
   domain->IntCofG.z /= nCoor;
}

void PairDomains(DOMAIN *domains)
{
   DOMAIN *d1, *d2;
   
   for(d1=domains; d1!=NULL; NEXT(d1))
   {
      for(d2=domains; d2!=NULL; NEXT(d2))
      {
         REAL distCofG;
         VEC3F *c1, *c2;
         c1 = &(d1->CofG);
         c2 = &(d2->CofG);

         distCofG = DISTSQ(c1, c2);
#ifdef DEBUG
         printf("CofG Distance: %.3f\n", sqrt(distCofG));
#endif
         
         if((distCofG < COFGDISTCUTSQ) && (distCofG > 1.0))
         {
            REAL distInt;
            
            c1 = &(d1->IntCofG);
            c2 = &(d2->IntCofG);

            distInt = DISTSQ(c1, c2);
            if(distInt < INTDISTCUTSQ)
            {
               if(distInt < distCofG)
               {
                  if((distCofG < d1->pairCofGDistSq) ||
                     (distInt  < d1->pairIntDistSq)  ||
                     (distCofG < d2->pairCofGDistSq) ||
                     (distInt  < d2->pairIntDistSq))
                     
                  {
                     d1->pairCofGDistSq = distCofG;
                     d1->pairIntDistSq  = distInt;
                     d2->pairCofGDistSq = distCofG;
                     d2->pairIntDistSq  = distInt;
                     d1->pairedDomain = d2;
                     d2->pairedDomain = d1;
                  }
               }
            }
            
#ifdef DEBUG
            printf("Interface Distance: %.3f\n", sqrt(distInt));
#endif
            
         }
         
         
      }
#ifdef DEBUG
      printf("\n");
#endif      
   }
   
}



void WriteDomains(DOMAIN *domains, char *filestem)
{
   DOMAIN *d, *pd;
   static int domCount = 0;
   int    i;

   
   for(d=domains; d!=NULL; NEXT(d))
      d->used = FALSE;
   
   for(d=domains; d!=NULL; NEXT(d))
   {
      if(!d->used)
      {
         FILE *fp;
         char outFile[MAXBUFF+1];
         
         sprintf(outFile, "%s_%d.pdb", filestem, domCount++);
         outFile[MAXBUFF] = '\0';

         if((fp = fopen(outFile, "w"))!=NULL)
         {
            PDB *p;
            d->used = TRUE;
            /* Write this domain */
            for(p=d->startRes; p!=d->stopRes; NEXT(p))
            {
               blWritePDBRecord(fp, p);
            }
            fprintf(fp,"TER   \n");
            /* Write partner domain */
            if((pd = d->pairedDomain) != NULL)
            {
               pd->used = TRUE;
               for(p=pd->startRes; p!=pd->stopRes; NEXT(p))
               {
                  blWritePDBRecord(fp, p);
               }
               fprintf(fp,"TER   \n");
            }
            if(!gNoAntigen)
            {
               /* Write antigen chains */
               for(i=0; i<d->nAntigenChains; i++)
               {
                  PDBCHAIN *chain = d->antigenChains[i];
                  for(p=chain->start; p!=chain->stop; NEXT(p))
                  {
                     blWritePDBRecord(fp, p);
                  }
                  fprintf(fp,"TER   \n");
               }
            }
            
            fclose(fp);
         }
      }
   }
}

void FlagAntigens(DOMAIN *domains, PDBSTRUCT *pdbs)
{
   DOMAIN *d;

   printf("\n***Looking for non-het antigens\n");

   for(d=domains; d!=NULL; NEXT(d))
   {
      d->used = FALSE;
   }
   
   for(d=domains; d!=NULL; NEXT(d))
   {
      if(!d->used)
      {
         d->used = TRUE;
         if(d->pairedDomain != NULL)
            d->pairedDomain->used = TRUE;
         
         CheckAntigenContacts(d, pdbs);
      }
   }
}


void SetChainTypes(PDBCHAIN *chains)
{
   PDBCHAIN *c;
   PDB      *p;

   for(c=chains; c!=NULL; NEXT(c))
   {
      c->extras = CHAINTYPE_HET;
      for(p=c->start; p!=c->stop; NEXT(p))
      {
         if(!strncmp(p->record_type, "ATOM  ", 6))
         {
            c->extras = CHAINTYPE_ATOM;
            break;
         }
      }
   }
}

BOOL inIntArray(int value, int *array, int arrayLen)
{
   int i;

   for(i=0; i<arrayLen; i++)
   {
      if(array[i] == value)
         return(TRUE);
   }

   return(FALSE);
}

BOOL regionsMakeContact(PDB *start1, PDB *stop1, PDB *start2, PDB *stop2)
{
   PDB *p, *q;
   
   for(p=start1; p!=stop1; NEXT(p))
   {
      for(q=start2; q!=stop2; NEXT(q))
      {
         if(DISTSQ(p,q) < CONTACTDISTSQ)
         {
#ifdef DEBUG
            printf("Contact %s%d%-s.%-s with %s%d%-s.%-s (%.3f)\n",
                   p->chain, p->resnum, p->insert, p->atnam,
                   q->chain, q->resnum, q->insert, q->atnam,
                   sqrt(DISTSQ(p,q)));
#endif
            return(TRUE);

         }
      }
   }

   return(FALSE);
}

void CheckAntigenContacts(DOMAIN *domain, PDBSTRUCT *pdbs)
{
   PDB      *p, *q, *nextResP, *nextResQ;
   PDBCHAIN *chain;
   int      nContacts = 0;
   DOMAIN   *pairedDomain = domain->pairedDomain;

   
   domain->nAntigenChains       = 0;
   if(pairedDomain != NULL)
      pairedDomain->nAntigenChains = 0;

   /* Go through each of the ATOM chains */
   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      if(chain->extras == CHAINTYPE_ATOM)
      {
         nContacts = 0;
         /* If this is not the chain for this domain or the paired domain,
            then it's a potential antigen
         */
         if((chain != domain->chain) &&
            ((pairedDomain == NULL) ||
             (chain != pairedDomain->chain)))
         {
            int resnum = 0;

#ifdef DEBUG
            printf("Checking domain %d (chain %s) against chain %s\n",
                   domain->domainNumber, domain->chain->chain, chain->chain);
#endif
            
            /* Check this domain for contacts a residue at a time */
            for(p=domain->startRes; p!=domain->stopRes; p=nextResP)
            {
               nextResP = blFindNextResidue(p);

               /* If this is a CDR residue */
               if(inIntArray(resnum, domain->CDRRes, domain->nCDRRes))
               {
                  for(q=chain->start; q!=chain->stop; q=nextResQ)
                  {
                     nextResQ = blFindNextResidue(q);
                     
                     if(regionsMakeContact(p, nextResP, q, nextResQ))
                     {
#ifdef DEBUG
                        printf("Domain %d (chain %s) makes %d contacts with chain %s\n",
                               domain->domainNumber, domain->chain->chain, nContacts, chain->chain);
#endif                           
                        if(++nContacts >= MINAGCONTACTS)
                        {
                           if(domain->nAntigenChains < MAXANTIGEN)
                              domain->antigenChains[domain->nAntigenChains++] = chain;
                           if((pairedDomain != NULL) && (pairedDomain->nAntigenChains < MAXANTIGEN))
                              pairedDomain->antigenChains[pairedDomain->nAntigenChains++] = chain;
#ifndef DEBUG
                           goto break1;
#endif
                        }
                     }
                  }
               }

               resnum++;
            }
            
            /* Check partner domain for contacts */
            if(pairedDomain != NULL)
            {
#ifdef DEBUG
               printf("Checking partner domain %d (chain %s) against chain %s\n",
                      pairedDomain->domainNumber, pairedDomain->chain->chain, chain->chain);
#endif
               resnum    = 0;
               for(p=pairedDomain->startRes; p!=pairedDomain->stopRes; p=nextResP)
               {
                  nextResP = blFindNextResidue(p);
                     
                  if(inIntArray(resnum, pairedDomain->CDRRes, pairedDomain->nCDRRes))
                  {
                     for(q=chain->start; q!=chain->stop; q=nextResQ)
                     {
                        nextResQ = blFindNextResidue(q);
                        
                        if(regionsMakeContact(p, nextResP, q, nextResQ))
                        {
#ifdef DEBUG
                           printf("Domain %d (chain %s) makes %d contacts with chain %s\n",
                                  pairedDomain->domainNumber, pairedDomain->chain->chain,
                                  nContacts, chain->chain);
#endif

                           if(++nContacts >= MINAGCONTACTS)
                           {
                              if(domain->nAntigenChains < MAXANTIGEN)
                                 domain->antigenChains[domain->nAntigenChains++] = chain;
                              if((pairedDomain != NULL) && (pairedDomain->nAntigenChains < MAXANTIGEN))
                                 pairedDomain->antigenChains[pairedDomain->nAntigenChains++] = chain;
#ifndef DEBUG
                              goto break1;
#endif
                           }
                        }
                     }
                  }
                  resnum++;
               }
            }
         break1:
            continue;
         }
      }
   }
}


