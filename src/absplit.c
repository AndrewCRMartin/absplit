/* TODO:
   WriteSEQRES needs to deal with chain labels for MODRES
*/
/************************************************************************/
/**

   Program:    absplit
   \file       absplit.c
   
   \version    V1.0
   \date       17.09.21   
   \brief      Split an antibody PDB file into Fvs with antigens
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2021-22
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
#include "bioplib/array.h"

#include "absplit.h"

/************************************************************************/
/* Defines and macros
*/
#define PROGNAME        "abYsplit"
#define VERSION         "V1.0"
#define MAXBUFF         240
#define HUGEBUFF        10000
#define MAXSEQ          10000
#define ABTHRESHOLD     0.45
#define GAPOPENPENALTY  5
#define GAPEXTPENALTY   2
#define SCOREMATRIX     "BLOSUM62"
#define MAXINTERFACE    20
#define MAXCDRRES       160
#define MAXRESID        16
#define MINSEQLEN       50
#define COFGDISTCUTSQ   1225.0 /* 35^2 - for possible VH/VL pairs       */
#define INTDISTCUTSQ    400.0  /* 20^2 - for VH/VL interface contact    */
#define CONTACTDISTSQ   36.0   /* 6^2  - for antigen contacts           */
#define MINAGCONTACTS   14     /* Tweaked with CONTACTDISTSQ to get Ag for
                                  6o8d but not too much for 1dee
                                  Currently includes xtal packing for 
                                  1a6v                                  */
#define MINHETATOMS     8
#define MAXANTIGEN      16
#define MAXCHAINS       80
#define MAXCHAINLABEL   blMAXCHAINLABEL
#define MAXHETANTIGEN   160
#define CHAINTYPE_PROT  (APTR)1
#define CHAINTYPE_NUCL  (APTR)2
#define CHAINTYPE_HET   (APTR)3
#define RES_WRITTEN_NO  (APTR)0
#define RES_WRITTEN_YES (APTR)1

typedef struct _domain
{
   int   domainNumber,
         startSeqRes,
         lastSeqRes,
         interface[MAXINTERFACE],
         nInterface,
         nAntigenChains,
         nHetAntigen,
         nCDRRes,
         CDRRes[MAXCDRRES];
   char  domSeq[MAXSEQ],
         newAbChainLabel[8],
         newAgChainLabels[MAXANTIGEN][8],
         chainType;
   PDB   *startRes,
         *lastRes,
         *stopRes;
   VEC3F CofG,
         IntCofG;
   REAL  pairIntDistSq,
         pairCofGDistSq;
   BOOL  used;
   PDBRESIDUE *hetAntigen[MAXHETANTIGEN];
   PDBCHAIN   *chain,
              *antigenChains[MAXANTIGEN];
   struct _domain *pairedDomain;
   struct _domain *next;
}  DOMAIN;

/************************************************************************/
/* Globals
*/
BOOL gVerbose   = FALSE;
BOOL gQuiet     = FALSE;
BOOL gNoAntigen = FALSE;


/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile);
void UsageDie(void);
BOOL ProcessFile(WHOLEPDB *wpdb, char *infile, FILE *dataFp);
DOMAIN *FindVHVLDomains(WHOLEPDB *wpdb, PDBCHAIN *chain, FILE *dataFp,
                        DOMAIN *domains);
void GetSequenceForChain(WHOLEPDB *wpdb, PDBCHAIN *chain, char *sequence);
void ExePathName(char *str, BOOL pathonly);
BOOL CheckAndMask(char *sequence, FILE *dataFp, PDBCHAIN *chain,
                  DOMAIN **pDomains);
FILE *OpenSequenceDataFile(void);
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2);
void MaskAndAssignDomain(char *seq, PDBCHAIN *chain, char *bestMatch,
                         char *aln1, char *aln2, DOMAIN **pDomains);
void SetChainAsLightOrHeavy(DOMAIN *domain, char *header);
void SetIFResidues(DOMAIN *domain, char *header, char *seqAln, char *refAln);
void SetCDRResidues(DOMAIN *domain, char *header, char *seqAln, char *refAln);
void PrintDomains(DOMAIN *domains);
void SetDomainBoundaries(DOMAIN *domain);
void PairDomains(DOMAIN *domains);
void WriteDomains(WHOLEPDB *wpdb, DOMAIN *domains, char *filestem);
BOOL FlagProteinAntigens(DOMAIN *domains, PDBSTRUCT *pdbs);
BOOL IsNonPeptideHet(WHOLEPDB *wpdb, PDBRESIDUE *res);
BOOL CheckAntigenContacts(DOMAIN *domain, PDBSTRUCT *pdbs);
void SetChainAsAtomOrHetatm(PDBCHAIN *chains);
BOOL inIntArray(int value, int *array, int arrayLen);
void GetSequenceForChainSeqres(WHOLEPDB *wpdb, PDBCHAIN *chain,
                               char *sequence);
BOOL RegionsMakeContact(PDB *start1, PDB *stop1, PDB *start2, PDB *stop2);
void FlagHetAntigenChains(DOMAIN *domains, PDBSTRUCT *pdbs);
void FlagHetAntigenResidues(WHOLEPDB *wpd, DOMAIN *domains,
                            PDBSTRUCT *pdbs);
PDB *RelabelAntibodyChain(DOMAIN *domain,
                          BOOL *lowerCaseLight, BOOL *lowerCaseHeavy,
                          char *remark950);
PDB *RelabelAntigenChains(DOMAIN *domain, char *remark950);
void WriteSeqres(FILE *fp, WHOLEPDB *wpdb, DOMAIN *d);
int CountResidueAtoms(PDBRESIDUE *res);
char *blFixSequenceWholePDB(WHOLEPDB *wpdb, char **outChains,
                            BOOL ignoreSeqresForMissingChains,
                            BOOL upper, BOOL quiet, char *label);
int TransferResnum(int refResnum, char *seqAln, char *refAln);
int RealSeqLen(char *seq);
BOOL IsStandardResidue(PDBRESIDUE *res);
int FindLastAlignmentPosition(char *refAln);
int IsKeyResidue(int seqPos, int *refKeys,
                 char *seqAln, char *refAln);

BOOL DomainSequenceMatchesChainSequence(DOMAIN *domain, PDBCHAIN *chain);





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
            
            /* Read the mutation matrix                                 */
            blReadMDM(SCOREMATRIX);
            
            /* Do the real work of processing this file                 */
            if(!ProcessFile(wpdb, infile, dataFp))
            {
               fprintf(stderr,"Error (%s): Unable to split PDB into \
chains\n", PROGNAME);
               exit(1);
            }
            
            blFreeWholePDB(wpdb);
         }
         else
         {
            fprintf(stderr,"Error (%s): Nothing read from input file \
(%s)\n", PROGNAME, infile);
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
            gQuiet   = FALSE;
            break;
         case 'q':
            gVerbose = FALSE;
            gQuiet   = TRUE;
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
   printf("%s %s (c) UCL, Prof. Andrew C.R. Martin\n", PROGNAME, VERSION);

   printf("\nUsage: abysplit [-v][-q][-n] file.pdb       \n");
   printf("           -v Verbose\n");
   printf("           -q Quiet\n");
   printf("           -n Do not include the antigen in the output\n");
   printf("\nTakes a PDB file containing one or more antibodies and \
splits it into\n");
   printf("separate antibody files, retaining the antigen in each. \
Handles both\n");
   printf("ATOM (protein) antigens and HETATM (hapten) antigens. Deals \
with scFvs\n");
   printf("and other complex non-standard structures.\n\n");

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
      int length = stop-start;
      strncpy(filestem, start, (int)length);
      filestem[(int)length] = '\0';
   }
}


/************************************************************************/
BOOL ProcessFile(WHOLEPDB *wpdb, char *infile, FILE *dataFp)
{
   char      filestem[MAXBUFF],
             *sequence;
   PDBSTRUCT *pdbs;
   PDB       *pdb;
   int       nAtoms;

   GetFilestem(infile, filestem);

   /* Remove any waters from the PDB linked list                        */
   pdb = blStripWatersPDBAsCopy(wpdb->pdb, &nAtoms);
   FREELIST(wpdb->pdb, PDB);
   wpdb->pdb = pdb; wpdb->natoms = nAtoms;
   
   if((pdbs = blAllocPDBStructure(wpdb->pdb))!=NULL)
   {
      char **outChains = NULL;
      
      if((outChains = (char **)blArray2D(sizeof(char),
                                         MAXCHAINS,
                                         blMAXCHAINLABEL))==NULL)
      {  /* TODO This could be the actual number of chains instead of 
            MAXCHAINS 
            Also needs to be freed at the end.
         */
         fprintf(stderr,"Error: No memory for outChains array\n");
         return(FALSE);
      }
      
      if((sequence = blFixSequenceWholePDB(wpdb, outChains,
                                           TRUE, FALSE,
                                           TRUE, NULL))!=NULL)
      {
         /* NOTE! The output of blFixSequenceWholePDB isn't used        */
         PDBCHAIN *chain;
         DOMAIN   *domains = NULL;

#ifdef DEBUG
         fprintf(stderr, "Sequence:\n%s\n", sequence);
#endif
         
         SetChainAsAtomOrHetatm(pdbs->chains);
         
         for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
         {
#ifdef DEBUG
            printf("Ptr: %ld Type: %s\n", (ULONG)chain->start,
                   chain->start->record_type);
#endif
            if(chain->extras == CHAINTYPE_PROT)
            {
               printf("\n***Handling chain: %s\n", chain->chain);
               domains = FindVHVLDomains(wpdb, chain, dataFp, domains);
            }
         }
         
         if(domains != NULL)
         {
            PairDomains(domains);
            
            FlagProteinAntigens(domains, pdbs);
            FlagHetAntigenChains(domains, pdbs);
            FlagHetAntigenResidues(wpdb, domains, pdbs);
            
            PrintDomains(domains);
            WriteDomains(wpdb, domains, filestem);
            
            FREELIST(domains, DOMAIN);
            blFreePDBStructure(pdbs);
         }
         else
         {
            fprintf(stderr,"Error (abYsplit): no antibody domains \
found\n");
            return(FALSE);
         }
      }
      else
      {
         return(FALSE);
      }
   }
   else
   {
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
void GetSequenceForChain(WHOLEPDB *wpdb, PDBCHAIN *chain, char *sequence)
{
   int        i=0;
   PDBRESIDUE *r;

   if((chain->extras == CHAINTYPE_PROT) ||
      (chain->extras == CHAINTYPE_NUCL))
   {
      for(r=chain->residues; r!=NULL; NEXT(r))
      {
         if(IsStandardResidue(r))
         {
#ifdef DEBUG
            fprintf(stderr, "%s\n", r->resnam);
#endif
            sequence[i++] = blThrone(r->resnam);
         }
      }
   }
   
   sequence[i] = '\0';
}

/************************************************************************/
/* 
   From 
   https://www.linuxquestions.org/questions/showthread.php?threadid=273614
*/
void ExePathName(char *str, BOOL pathonly)
{
  FILE *fp;
  char buf[PATH_MAX+100],
       *p;

  *str    = '\0';
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
     if(p!=NULL)
        *p='\0';
  }
}


/************************************************************************/
FILE *OpenSequenceDataFile(void)
{
   char pathname[PATH_MAX];
   
   ExePathName(pathname, TRUE);
   strcat(pathname,ABSEQFILE);
   
   return(fopen(pathname, "r"));
}


/************************************************************************/
/*>REAL CompareSeqs(char *seqresSeq, char *refSeq, 
                    char *alignSeqres, char *alignRef)
   -----------------------------------------------
*//**
   \param[in]   seqresSeq     the sequence of interest
   \param[in]   refSeq        the database sequence
   \param[out]  alignSeqres   Alignment of our sequence
   \param[out]  alignRef      Alignment of database sequence
   \return                    Score for alignment

   - 31.03.20 Original   By: ACRM
   - 20.09.21 Modified to use affine alignment and mutation matrix
*/
REAL CompareSeqs(char *seqresSeq, char *refSeq,
                 char *alignSeqres, char *alignRef)
{
   int  score;
   int  alignLen;
   int  bestPossibleScore,
        bestPossibleScore1 = 100000,
        bestPossibleScore2 = 100000;

   bestPossibleScore1 = blAffinealign(refSeq, strlen(refSeq),
                                      refSeq, strlen(refSeq),
                                      FALSE,          /* verbose        */
                                      FALSE,          /* identity       */
                                      GAPOPENPENALTY, /* penalty        */
                                      GAPEXTPENALTY,  /* extension      */
                                      alignSeqres,
                                      alignRef,
                                      &alignLen);

   if(strlen(seqresSeq) >= 50)
   {
      bestPossibleScore2 = blAffinealign(seqresSeq, strlen(seqresSeq),
                                         seqresSeq, strlen(seqresSeq),
                                         FALSE,          /* verbose     */
                                         FALSE,          /* identity    */
                                         GAPOPENPENALTY, /* penalty     */
                                         GAPEXTPENALTY,  /* extension   */
                                         alignSeqres,
                                         alignRef,
                                         &alignLen);
   }
   
   bestPossibleScore = MIN(bestPossibleScore1, bestPossibleScore2);

   if(bestPossibleScore <= 0)
      return(0);
   
   score = blAffinealign(seqresSeq, strlen(seqresSeq),
                         refSeq, strlen(refSeq),
                         FALSE,          /* verbose                     */
                         FALSE,          /* identity                    */
                         GAPOPENPENALTY, /* penalty                     */
                         GAPEXTPENALTY,  /* extension                   */
                         alignSeqres,
                         alignRef,
                         &alignLen);

   alignSeqres[alignLen]  = alignRef[alignLen] = '\0';

#ifdef DEBUG   
   fprintf(stderr, "\n>>>%s\n", alignSeqres);
   fprintf(stderr, ">>>%s\n",   alignRef);
#endif

   return((REAL)score / (REAL)bestPossibleScore);
}


/************************************************************************/
DOMAIN *FindVHVLDomains(WHOLEPDB *wpdb, PDBCHAIN *chain, FILE *dataFp,
                        DOMAIN *domains)
{
   char sequence[MAXSEQ];
   
/*   GetSequenceForChainSeqres(wpdb, chain, sequence); */
   
   GetSequenceForChain(wpdb, chain, sequence);
#ifndef DEBUG
   printf("Chain: %s Sequence: %s\n", chain->chain, sequence);
#endif
   while(TRUE)
   {
      if(!CheckAndMask(sequence, dataFp, chain, &domains)) break;
   }

   return(domains);
}

/************************************************************************/
int RealSeqLen(char *seq)
{
   char *chp;
   int seqlen = 0;

   for(chp=seq; *chp != '\0'; chp++)
   {
      if((*chp != 'X') && (*chp != '-'))
         seqlen++;
   }

   return(seqlen);
}

/************************************************************************/
BOOL CheckAndMask(char *seqresSeq, FILE *dbFp, PDBCHAIN *chain,
                  DOMAIN **pDomains)
{
   char        bestMatchFastaHeader[MAXBUFF+1];
   REAL        maxScore = 0.0;
   static char alignSeqres[HUGEBUFF+1],
               alignRef[HUGEBUFF+1];
   static char bestAlignSeqres[HUGEBUFF+1],
               bestAlignRef[HUGEBUFF+1];
   char        header[MAXBUFF+1];
   char        *refSeq = NULL;
   BOOL        found = FALSE;

   if(RealSeqLen(seqresSeq) < MINSEQLEN)
      return(FALSE);

   /* Find the best match in the reference sequences                    */
   rewind(dbFp);
   while((refSeq = blReadFASTA(dbFp, header, MAXBUFF))!=NULL)
   {
      REAL score;
         
      score = CompareSeqs(seqresSeq, refSeq, alignSeqres, alignRef);
      if(score > maxScore)
      {
         maxScore = score;
         strncpy(bestMatchFastaHeader,  header, MAXBUFF);
         strncpy(bestAlignSeqres, alignSeqres, HUGEBUFF);
         strncpy(bestAlignRef, alignRef, HUGEBUFF);
      }
      free(refSeq);
   }

   /* If we found an antibody sequence                                  */
#ifndef DEBUG
   printf("MaxScore : %f\n", maxScore);
   printf("Sequence : %s\n", bestAlignSeqres);
   printf("Reference: %s\n", bestAlignRef);
#endif
   
   if(maxScore > ABTHRESHOLD)
   {
      if(gVerbose)
      {
         fprintf(stderr, "Best match: %s Score: %.4f\n",
                 bestMatchFastaHeader, maxScore);
         fprintf(stderr, "SEQ: %s\n",   bestAlignSeqres);
         fprintf(stderr, "REF: %s\n\n", bestAlignRef);
      }
      
      MaskAndAssignDomain(seqresSeq, chain, bestMatchFastaHeader,
                          bestAlignSeqres, bestAlignRef, pDomains);
#ifndef DEBUG
      printf("Masked   : %s\n", seqresSeq);
#endif
      found = TRUE;
   }
   
   return(found);
}


/************************************************************************/
/* Uses the label in the FASTA file to set the chain type as heavy
   or light
*/
void SetChainAsLightOrHeavy(DOMAIN *domain, char *fastaHeader)
{
   char *ptr;
   domain->chainType = '?';
   if((ptr = strchr(fastaHeader, '|'))!=NULL)
   {
      ptr--;
      domain->chainType = *ptr;
   }
}


/************************************************************************/
/* TODO!
   This needs to take the alignment and translate the residues numbers
   to the PDB sequential number instead of what is in the FASTA file
   header
*/
void SetIFResidues(DOMAIN *domain, char *fastaHeader, char *seqAln,
                   char *refAln)
{
   char *ptr1, *ptr2, *bar;
   char headerCopy[MAXBUFF+1];
   int  refResnum,
        seqResnum,
        nIFRes = 0,
        IFResidues[MAXINTERFACE];
   int  seqLen = RealSeqLen(seqAln);
      
   domain->nInterface = 0;
   strncpy(headerCopy, fastaHeader, MAXBUFF);
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
         sscanf(ptr1, "%d", &refResnum);
         IFResidues[nIFRes++] = refResnum;
         ptr1  = ptr2+1;
      }
      
      if((ptr2 = strchr(ptr1, ']'))!=NULL)
      {
         *ptr2 = '\0';
         sscanf(ptr1, "%d", &refResnum);
         IFResidues[nIFRes++] = refResnum;
      }
      IFResidues[nIFRes++] = -1;
      for(seqResnum=1; seqResnum<=seqLen; seqResnum++)
      {
         if(IsKeyResidue(seqResnum, IFResidues, seqAln, refAln))
            domain->interface[domain->nInterface++] = seqResnum;
      }
   }
}

/************************************************************************/
void SetCDRResidues(DOMAIN *domain, char *fastaHeader, char *seqAln,
                    char *refAln)
{
   char *ptr1, *ptr2,
        headerCopy[MAXBUFF+1];
   int  refResnum,
        seqResnum,
        nCDRRes = 0,
        CDRResidues[MAXCDRRES];
   int  seqLen = RealSeqLen(seqAln);

   domain->nCDRRes = 0;

   strncpy(headerCopy, fastaHeader, MAXBUFF);
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
            sscanf(ptr1, "%d", &refResnum);
            CDRResidues[nCDRRes++] = refResnum;
            ptr1=ptr2+1;
         }
      
         if((ptr2 = strchr(ptr1, ']'))!=NULL)
         {
            *ptr2 = '\0';
            sscanf(ptr1, "%d", &refResnum);
            CDRResidues[nCDRRes++] = refResnum;
         }
         CDRResidues[nCDRRes] = -1;

         for(seqResnum=1; seqResnum<=seqLen; seqResnum++)
         {
            if(IsKeyResidue(seqResnum, CDRResidues, seqAln, refAln))
               domain->CDRRes[domain->nCDRRes++] = seqResnum;
         }
      }
   }
}


/* seqPos numbered from 1 */
int IsKeyResidue(int seqPos, int *refKeys,
                 char *seqAln, char *refAln)
{
   int i, count, alnPos=0, refPos;
   
   /* Correct seqPos to the position in the alignment (numbered from 1) */
   for(i=0, count=0; i<strlen(seqAln); i++)
   {
      if(seqAln[i] != '-')
         count++;
      if(count == seqPos)
      {
         alnPos = i+1;
         break;
      }
   }

   /* Now count to this position in the reference sequence              */
   for(i=0, refPos=0; (i<alnPos && i<strlen(seqAln)); i++)
   {
      if(refAln[i] != '-')
         refPos++;
   }

   /* Now see if this position is in the list of key reference positions*/
   for(i=0; refKeys[i] >=0; i++)
   {
      if(refPos == refKeys[i])
         return(refPos);
   }
   
   return(0);
}

/************************************************************************/
int TransferResnum(int refResnum, char *seqAln, char *refAln)
{
   int atomResnum = 0,
       i;

   /* The position of interest is missing in the atom sequence          */
   if(seqAln[refResnum] == '-')
      return(-1);
   
   for(i=0; i<refResnum; i++)
   {
      if(seqAln[i] != '-')
         atomResnum++;
   }

   return(atomResnum);
}



/************************************************************************/
void PrintDomains(DOMAIN *domains)
{
   DOMAIN *d;

   printf("\n***Results\n");
   
   for(d=domains; d!=NULL; NEXT(d))
   {
      printf("Domain: %d Chain: %s Start: %d Stop: %d Type: %c \
PairsWithDomain: %d (Chain: %s)\n",
             d->domainNumber,
             d->chain->chain,
             d->startSeqRes,
             d->lastSeqRes,
             d->chainType,
             (d->pairedDomain==NULL)?0:d->pairedDomain->domainNumber,
             (d->pairedDomain==NULL)?"none":d->pairedDomain->chain->chain
         );
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
      if(d->startRes) blWritePDBRecord(stdout,d->startRes);
      if(d->lastRes)  blWritePDBRecord(stdout,d->lastRes);
#endif
   }
}


/************************************************************************/
void MaskAndAssignDomain(char *seq, PDBCHAIN *chain, char *fastaHeader,
                         char *seqAln, char *refAln, DOMAIN **pDomains)
{
   int    seqPos      = 0,
          alnPos      = 0,
          domSeqPos   = 0,
          firstAlnPos = 0,
          lastAlnPos  = 0;
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
      fprintf(stderr,"Error (%s): No memory for list of domains\n",
              PROGNAME);
      exit(1);
   }
   d->startSeqRes    = -1;
   d->lastSeqRes     = -1;
   d->chain          = chain;
   d->pairIntDistSq  = 100000000.0;
   d->pairCofGDistSq = 100000000.0;
   d->domainNumber   = (prevD==NULL)?1:prevD->domainNumber+1;
   d->pairedDomain   = NULL;
   d->nCDRRes        = 0;
   d->nInterface     = 0;
   d->nHetAntigen    = 0;
   d->nAntigenChains = 0;


#ifndef DEBUG_SET_CDR
   printf("SEQ      : %s\n", seqAln);
   printf("REF      : %s\n", refAln);
#endif

   lastAlnPos = FindLastAlignmentPosition(refAln);

   SetChainAsLightOrHeavy(d, fastaHeader);
   SetIFResidues(d,          fastaHeader, seqAln, refAln);
   SetCDRResidues(d,         fastaHeader, seqAln, refAln);
   
   for(seqPos=0, alnPos=0;
       ((seqPos<strlen(seqAln)) && (alnPos < lastAlnPos));
       seqPos++)
   {
      while(seqAln[alnPos] == '-')  /* Skip insertions                  */
      {
         alnPos++;
      }
      
      if(refAln[alnPos] != '\0')
      {
         if(d->startSeqRes == (-1))
            d->startSeqRes = seqPos;
#ifdef DEBUG
         printf("Seqpos %d SeqRes %c DomRes %c\n",
                seqPos, seqAln[alnPos], refAln[alnPos]);
#endif
         d->lastSeqRes = seqPos;
         d->domSeq[domSeqPos++] = seq[seqPos];
/*         seq[seqPos] = 'X'; */
      }
      alnPos++;
   }
   d->domSeq[domSeqPos] = '\0';

   /* Mask the sequence */
   firstAlnPos = -1;
   lastAlnPos  = -1;
   for(seqPos=0, alnPos=0;
       seqPos<strlen(seq) && alnPos<strlen(seqAln);
       alnPos++)
   {
      /* If this is an aligned position */
      if((seqAln[alnPos] != '-') &&
         (refAln[alnPos] != '-'))
      {
         seq[seqPos] = 'X';
      }
      if(seqAln[alnPos] != '-')
         seqPos++;
   }

   SetDomainBoundaries(d);
}

/************************************************************************/
void SetDomainBoundaries(DOMAIN *domain)
{
   int resnum = 0,
       nCoor  = 0,
       i;
   PDB *p,
       *nextRes;

   domain->startRes = NULL;
   domain->stopRes  = NULL;
   domain->lastRes  = NULL;
   
   /* Given the integer sequence positions, find the pointers to
      the PDB residues for these positions
   */
   for(p=domain->chain->start;
       p!=domain->chain->stop && p!=NULL;
       p=nextRes)
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
   nCoor          = 0;
   
   for(p=domain->startRes; p!=NULL && p!=domain->stopRes; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ", 4))
      {
#ifdef FUBAR
         fprintf(stderr, "Calc CofG %d: %s%d%s\n",
                 domain->domainNumber, p->chain, p->resnum, p->insert);
#endif
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


/************************************************************************/
void PairDomains(DOMAIN *domains)
{
   DOMAIN *d1, *d2;
   
   for(d1=domains; d1!=NULL; NEXT(d1))
   {
      for(d2=domains; d2!=NULL; NEXT(d2))
      {
         if(d1 != d2)
         {
            REAL distCofGSq;
            VEC3F *c1, *c2;
            c1 = &(d1->CofG);
            c2 = &(d2->CofG);
            
            distCofGSq = DISTSQ(c1, c2);

#ifdef FUBAR
            fprintf(stderr, "Chain1 %s; CofG1 %.3f %.3f %.3f; \
Chain2 %s; CofG2 %.3f %.3f %.3f; \
Dist %.3f\n",
                    d1->lastRes->chain,
                    c1->x, c1->y, c1->z, 
                    d2->startRes->chain,
                    c2->x, c2->y, c2->z,
                    sqrt(distCofGSq));
#endif
            if((distCofGSq < COFGDISTCUTSQ) && (distCofGSq > 1.0))
            {
               REAL distIntSq;
               
               c1 = &(d1->IntCofG);
               c2 = &(d2->IntCofG);
               
               distIntSq = DISTSQ(c1, c2);
               if(distIntSq < INTDISTCUTSQ)
               {
                  if(distIntSq < distCofGSq)
                  {
                     if(((distCofGSq < d1->pairCofGDistSq) ||
                         (distIntSq  < d1->pairIntDistSq)) &&
                        ((distCofGSq < d2->pairCofGDistSq) ||
                         (distIntSq  < d2->pairIntDistSq)))
                     {
                        d1->pairCofGDistSq = distCofGSq;
                        d1->pairIntDistSq  = distIntSq;
                        d2->pairCofGDistSq = distCofGSq;
                        d2->pairIntDistSq  = distIntSq;
                        d1->pairedDomain = d2;
                        d2->pairedDomain = d1;
#ifdef DEBUG
                        printf("*Paired domain %d with %d\n",
                               d1->domainNumber, d2->domainNumber);
#endif
                     }
                  }
               }
               
#ifdef DEBUG
               printf("CofG Distance (domain %d to %d): %.3f\n",
                      d1->domainNumber, d2->domainNumber,
                      sqrt(distCofGSq));
               printf("Interface Distance (domain %d to %d): %.3f\n\n",
                      d1->domainNumber, d2->domainNumber,
                      sqrt(distIntSq));
#endif
            }
         }
      }
#ifdef DEBUG
      printf("\n");
#endif      
   }
}


/************************************************************************/
void WriteDomains(WHOLEPDB *wpdb, DOMAIN *domains, char *filestem)
{
   DOMAIN     *d, *pd;
   static int domCount = 0;
   int        i;
   
   for(d=domains; d!=NULL; NEXT(d))
      d->used = FALSE;
   
   for(d=domains; d!=NULL; NEXT(d))
   {
      if(!d->used)
      {
         FILE *fp;
         char outFile[MAXBUFF+1],
              complex[8];
         BOOL written = FALSE;

         /* Assume not a complex                                        */
         complex[0] = '\0';

         for(i=0;i<d->nAntigenChains;i++)
         {
            if(d->antigenChains[i]->extras == CHAINTYPE_PROT)
            {
               strncat(complex, "P", 7-strlen(complex));
               break;
            }
         }
         
         for(i=0;i<d->nAntigenChains;i++)
         {
            if(d->antigenChains[i]->extras == CHAINTYPE_NUCL)
            {
               strncat(complex, "N", 7-strlen(complex));
               break;
            }
         }
         
         if(d->nHetAntigen)
            strncat(complex, "H", 7-strlen(complex));
         
         sprintf(outFile, "%s_%d%s.pdb", filestem, domCount++, complex);
         outFile[MAXBUFF] = '\0';

         if((fp = fopen(outFile, "w"))!=NULL)
         {
            PDB *pdb1, *pdb2, *pdb3, *p;
            char remark950Domain[100],
                 remark950Partner[100],
                 remark950Antigen[MAXANTIGEN * 100];
            BOOL lowerCaseLight = FALSE,
                 lowerCaseHeavy = FALSE;

            
            d->used = TRUE;

            pdb1 = RelabelAntibodyChain(d,
                                        &lowerCaseLight, &lowerCaseHeavy,
                                        remark950Domain);
            pdb2 = RelabelAntibodyChain(d->pairedDomain,
                                        &lowerCaseLight, &lowerCaseHeavy,
                                        remark950Partner);
            pdb3 = RelabelAntigenChains(d,
                                        remark950Antigen);

            fprintf(fp, "REMARK 950 CHAIN-TYPE  LABEL ORIGINAL\n");
            fprintf(fp, remark950Domain);
            fprintf(fp, remark950Partner);
            fprintf(fp, remark950Antigen);

            WriteSeqres(fp, wpdb, d);

            /* Write this domain                                        */
            for(p=pdb1; p!=NULL; NEXT(p))
            {
               blWritePDBRecord(fp, p);
            }
            fprintf(fp,"TER   \n");
            FREELIST(pdb1, PDB);
            
            /* Write partner domain                                     */
            if((pd = d->pairedDomain) != NULL)
            {
               pd->used = TRUE;
               for(p=pdb2; p!=NULL; NEXT(p))
               {
                  blWritePDBRecord(fp, p);
               }
               fprintf(fp,"TER   \n");
               FREELIST(pdb2, PDB);
            }
            
            if(!gNoAntigen)
            {
               PDB *prev;
               
               /* Write antigen chains                                  */
#ifdef OLD
               for(i=0; i<d->nAntigenChains; i++)
               {
                  PDBCHAIN *chain = d->antigenChains[i];
                  for(p=chain->start; p!=chain->stop; NEXT(p))
                  {
                     blWritePDBRecord(fp, p);
                  }
                  fprintf(fp,"TER   \n");
               }
#endif
               if(pdb3!=NULL)
               {
                  prev = pdb3;
                  for(p=pdb3; p!=NULL; NEXT(p))
                  {
                     if(!CHAINMATCH(p->chain, prev->chain))
                     {
                        fprintf(fp,"TER   \n");
                     }
                     blWritePDBRecord(fp, p);
                     prev=p;
                  }
                  fprintf(fp,"TER   \n");
               }
               
               
               /* Write any HET chains                                  */

               /* Clear flags to say a residue has been written         */
               for(i=0; i<d->nHetAntigen; i++)
               {
                  PDBRESIDUE *res = d->hetAntigen[i];
                  res->extras = RES_WRITTEN_NO;
               }
               if(pd!=NULL)
               {
                  for(i=0; i<pd->nHetAntigen; i++)
                  {
                     PDBRESIDUE *res = pd->hetAntigen[i];
                     res->extras = RES_WRITTEN_NO;
                  }
               }
               
               /* First domain                                          */
               for(i=0; i<d->nHetAntigen; i++)
               {
                  PDBRESIDUE *res = d->hetAntigen[i];
                  res->extras = RES_WRITTEN_YES;
                  fprintf(stderr,"Writing domain %d HET residue %s\n",
                          d->domainNumber, res->resid);
                  for(p=res->start; p!=res->stop; NEXT(p))
                  {
                     blWritePDBRecord(fp, p);
                  }
               }
               if(d->nHetAntigen)
                  fprintf(fp,"TER   \n");

               /* Partner domain                                        */
               if(pd!=NULL)
               {
                  for(i=0; i<pd->nHetAntigen; i++)
                  {
                     PDBRESIDUE *res = pd->hetAntigen[i];
                     if(res->extras == RES_WRITTEN_NO)
                     {
                        fprintf(stderr,
                                "Writing domain %d HET residue %s\n",
                                pd->domainNumber, res->resid);
                        written = TRUE;
                        for(p=res->start; p!=res->stop; NEXT(p))
                        {
                           blWritePDBRecord(fp, p);
                        }
                     }
                  }
                  if(written)
                     fprintf(fp,"TER   \n");
               }
               
            }
            
            fclose(fp);
         }
      }
   }
}

/************************************************************************/
BOOL FlagProteinAntigens(DOMAIN *domains, PDBSTRUCT *pdbs)
{
   DOMAIN *d;
   BOOL   foundAntigen = FALSE;

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
         
         foundAntigen |= CheckAntigenContacts(d, pdbs);
      }
   }
   return(foundAntigen);
}


/************************************************************************/
void SetChainAsAtomOrHetatm(PDBCHAIN *chains)
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
            if(!strncmp(p->resnam, "  U", 3) ||
               !strncmp(p->resnam, "  A", 3) ||
               !strncmp(p->resnam, "  C", 3) ||
               !strncmp(p->resnam, "  G", 3) ||
               !strncmp(p->resnam, " DT", 3) ||
               !strncmp(p->resnam, " DA", 3) ||
               !strncmp(p->resnam, " DC", 3) ||
               !strncmp(p->resnam, " DG", 3))
            {
               c->extras = CHAINTYPE_NUCL;
            }
            else
            {
               c->extras = CHAINTYPE_PROT;
            }
            
            break;
         }
      }
   }
}


/************************************************************************/
BOOL InIntArray(int value, int *array, int arrayLen)
{
   int i;

   for(i=0; i<arrayLen; i++)
   {
      if(array[i] == value)
         return(TRUE);
   }

   return(FALSE);
}


/************************************************************************/
BOOL RegionsMakeContact(PDB *start1, PDB *stop1, PDB *start2, PDB *stop2)
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


/************************************************************************/
BOOL CheckAntigenContacts(DOMAIN *domain, PDBSTRUCT *pdbs)
{
   PDB      *p, *q,
            *nextResP, *nextResQ;
   PDBCHAIN *chain;
   int      nContacts     = 0;
   DOMAIN   *pairedDomain = domain->pairedDomain;
   BOOL     foundAntigen  = FALSE;
   
   domain->nAntigenChains = 0;
   if(pairedDomain != NULL)
      pairedDomain->nAntigenChains = 0;

   /* Go through each of the ATOM chains                                */
   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      if((chain->extras == CHAINTYPE_PROT) ||
         (chain->extras == CHAINTYPE_NUCL))
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

#ifdef DEBUG_AG_CONTACTS
            printf("Checking domain %d (chain %s) against chain %s\n",
                   domain->domainNumber, domain->chain->chain,
                   chain->chain);
#endif

            /* Now check whether the sequence of this chain matches
               that of the domain or of the paired domain. If so this 
               must be crystal packing rather than an antigen interaction.
            */
            if(DomainSequenceMatchesChainSequence(domain, chain))
            {
               printf("Crystal packing: Chain %s = Domain %d (chain %s)\n",
                      chain->chain, domain->domainNumber,
                      domain->chain->chain);
               continue;
            }
            
            if(DomainSequenceMatchesChainSequence(pairedDomain, chain))
            {
               printf("Crystal packing: Chain %s = Domain %d (chain %s)\n",
                      chain->chain, pairedDomain->domainNumber,
                      pairedDomain->chain->chain);
               continue;
            }
            
            /* Check this domain for contacts a residue at a time       */
            for(p=domain->startRes; p!=domain->stopRes; p=nextResP)
            {
               nextResP = blFindNextResidue(p);

               /* If this is a CDR residue                              */
               if(InIntArray(resnum, domain->CDRRes, domain->nCDRRes))
               {
                  for(q=chain->start; q!=chain->stop; q=nextResQ)
                  {
                     nextResQ = blFindNextResidue(q);
                     
                     if(RegionsMakeContact(p, nextResP, q, nextResQ))
                     {
#ifdef DEBUG_AG_CONTACTS
                        printf("Domain %d (chain %s) makes %d contacts \
with chain %s\n",
                               domain->domainNumber,
                               domain->chain->chain,
                               nContacts, chain->chain);
#endif                           
                        if(++nContacts >= MINAGCONTACTS)
                        {
                           foundAntigen = TRUE;
                           if(domain->nAntigenChains < MAXANTIGEN)
                              domain->antigenChains[domain->nAntigenChains++] = chain;
                           if((pairedDomain != NULL) &&
                              (pairedDomain->nAntigenChains < MAXANTIGEN))
                              pairedDomain->antigenChains[pairedDomain->nAntigenChains++] = chain;
                           goto break1;
                        }
                     }
                  }
               }

               resnum++;
            }
            
            /* Check partner domain for contacts                        */
            if(pairedDomain != NULL)
            {
#ifdef DEBUG_AG_CONTACTS
               printf("Checking partner domain %d (chain %s) against \
chain %s\n",
                      pairedDomain->domainNumber,
                      pairedDomain->chain->chain, chain->chain);
#endif
               resnum = 0;
               for(p=pairedDomain->startRes;
                   p!=pairedDomain->stopRes;
                   p=nextResP)
               {
                  nextResP = blFindNextResidue(p);
                     
                  if(InIntArray(resnum, pairedDomain->CDRRes,
                                pairedDomain->nCDRRes))
                  {
                     for(q=chain->start; q!=chain->stop; q=nextResQ)
                     {
                        nextResQ = blFindNextResidue(q);
                        
                        if(RegionsMakeContact(p, nextResP, q, nextResQ))
                        {
#ifdef DEBUG_AG_CONTACTS
                           printf("Domain %d (chain %s) makes %d \
contacts with chain %s\n",
                                  pairedDomain->domainNumber,
                                  pairedDomain->chain->chain,
                                  nContacts, chain->chain);
#endif

                           if(++nContacts >= MINAGCONTACTS)
                           {
                              foundAntigen = TRUE;
                              if(domain->nAntigenChains < MAXANTIGEN)
                                 domain->antigenChains[domain->nAntigenChains++] = chain;
                              if((pairedDomain != NULL) &&
                                 (pairedDomain->nAntigenChains < MAXANTIGEN))
                                 pairedDomain->antigenChains[pairedDomain->nAntigenChains++] = chain;
                              goto break1;
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

   return(foundAntigen);
}

void GetSequenceForChainSeqres(WHOLEPDB *wpdb, PDBCHAIN *chain,
                               char *sequence)
{
   static MODRES *modres = NULL;
   static char   *seqres = NULL,
                 *seqresCopy = NULL,
                 **seqchains = NULL;
   int chainNum = 0;
   BOOL foundChain = FALSE;
   

   if(seqres == NULL)
   {
      if((seqchains =
          (char **)blArray2D(sizeof(char), MAXCHAINS,
                             MAXCHAINLABEL))!=NULL)
      {
         modres = blGetModresWholePDB(wpdb);
         seqres = blGetSeqresAsStringWholePDB(wpdb, seqchains,
                                              modres, FALSE);
         if(seqres)
         {
            if((seqresCopy = (char *)malloc((strlen(seqres)+1) *
                                            sizeof(char)))==NULL)
            {
               fprintf(stderr, "No memory for copy the sequence\n");
               exit(1);
            }
         }
      }
   }

   if(seqres)
   {
      strcpy(seqresCopy, seqres);
      for(chainNum=0; seqchains[chainNum][0] != '\0'; chainNum++)
      {
         if(CHAINMATCH(chain->chain, seqchains[chainNum]))
         {
            foundChain = TRUE;
            break;
         }
      }
   }

   if(foundChain)
   {
      int i;
      char *chp = seqresCopy;
      for(i=0; i<chainNum; i++)
      {
         if(chp)
         {
            if((chp = strchr(chp, '*'))!=NULL)
            {
               *chp = '\0';
               chp++;
            }
         }
      }
      if(chp)
      {
         char *star;
         if((star = strchr(chp, '*'))!=NULL)
         {
            *star = '\0';
         }
      }

      for(i=0; chp[i]; i++)
      {
         sequence[i] = chp[i];
      }
      sequence[i] = '\0';
   }
   else
   {
      GetSequenceForChain(wpdb, chain, sequence);
   }
}
   
   
/************************************************************************/
void FlagHetAntigenChains(DOMAIN *domains, PDBSTRUCT *pdbs)
{
   DOMAIN   *d;
   PDBCHAIN *c;
   
   printf("\n***Looking for HET antigen chains\n");

   /* Step through the chains                                           */
   for(c=pdbs->chains; c!=NULL; NEXT(c))
   {
      /* If it's a HET chain                                            */
      if(c->extras == CHAINTYPE_HET)
      {
         PDBRESIDUE *r;

         /* Go through the residues in this HET chain                   */
         for(r=c->residues; r!=NULL; NEXT(r))
         {
            /* If it isn't a water and it has more than 3 atoms         */
            if(!ISWATER(r) && (CountResidueAtoms(r) >= MINHETATOMS))
            {
               /* Go through the antibody domains                       */
               for(d=domains; d!=NULL; NEXT(d))
               {
                  PDB *pc, *pd, *pdRes, *nextResP;
                  int resnum = 0;

                  /* Go through the atoms in this non-water HET residue */
                  for(pc=r->start; pc!=r->stop; NEXT(pc))
                  {
                     /* Go through the residues in this antibody domain */
                     for(pdRes=d->startRes;
                         pdRes!=d->stopRes;
                         pdRes=nextResP)
                     {
                        nextResP = blFindNextResidue(pdRes);
                        resnum++;
                        
                        /* If this is a CDR residue                     */
                        if(InIntArray(resnum, d->CDRRes, d->nCDRRes))
                        {
                           for(pd=pdRes; pd!=nextResP; NEXT(pd))
                           {
                              /* If this atom is close enough to the HET 
                                 residue 
                              */
                              if(DISTSQ(pc, pd) < CONTACTDISTSQ)
                              {
                                 printf("HET group %s%d%s contacts \
Domain %d\n",
                                        pc->chain, pc->resnum, pc->insert,
                                        d->domainNumber);
#ifndef DEBUG
                                 printf("%s%d%s.%s contacts %s%d%s.%s\n",
                                        pc->chain, pc->resnum,
                                        pc->insert, pc->atnam,
                                        pd->chain, pd->resnum,
                                        pd->insert, pd->atnam);
#endif
                                 /* Flag this HET residue as making 
                                    contact   
                                 */
                                 if(d->nHetAntigen < MAXHETANTIGEN)
                                 {
                                    d->hetAntigen[d->nHetAntigen++] = r;
#ifdef DEBUG
                                    printf("Stored domain %d residue \
%s\n",
                                           d->domainNumber, r->resid);
#endif
                                 }
                                 goto lastatom;
                              }  /* In distance range                   */
                           }  /* Step through atoms in this residue     */
                        } /* Is a CDR residue                           */
                     } /* Step through residues in this domain          */
                  }  /* step through atoms in the non-water HET residue */
               lastatom:
                  continue;
               }  /* step through domains                               */
            }  /* not water                                             */
         }  /* step through residues in this HET chain                  */
      }  /* Is a HET chain                                              */
   }  /* Step through chains                                            */
}


/************************************************************************/
void FlagHetAntigenResidues(WHOLEPDB *wpdb, DOMAIN *domains,
                            PDBSTRUCT *pdbs)
{
   DOMAIN   *d;
   PDBCHAIN *c;
   
   printf("\n***Looking for HET antigen residues\n");

   /* Step through the chains                                           */
   for(c=pdbs->chains; c!=NULL; NEXT(c))
   {
      /* If it's not a HET chain                                        */
      if(c->extras != CHAINTYPE_HET)
      {
         PDBRESIDUE *r;
         /* Go through the residues in this non-HET chain               */
         for(r=c->residues; r!=NULL; NEXT(r))
         {
            /* If it's a HETATM group and has no peptide backbone and 
               isn't just an ion
            */
            if(IsNonPeptideHet(wpdb,r))
            {
               if(!ISWATER(r)) /* If it isn't a water                   */
               {
                  /* Go through the antibody domains                    */
                  for(d=domains; d!=NULL; NEXT(d))
                  {
                     PDB *pc, *pd, *pdRes, *nextResP;
                     int resnum = 0;
                     /* Go through the atoms in this non-water, 
                        non-peptide HET residue 
                     */
                     for(pc=r->start; pc!=r->stop; NEXT(pc))
                     {
                        /* Go through the residues in this antibody domain */
                        for(pdRes=d->startRes;
                            pdRes!=d->stopRes;
                            pdRes=nextResP)
                        {
                           nextResP = blFindNextResidue(pdRes);
                           resnum++;
                           
                           /* If this is a CDR residue                     */
                           if(InIntArray(resnum, d->CDRRes, d->nCDRRes))
                           {
                              
                              for(pd=pdRes; pd!=nextResP; NEXT(pd))
                              {
                                 /* If this atom is close enough to the HET 
                                    residue 
                                 */
                                 if(DISTSQ(pc, pd) < CONTACTDISTSQ)
                                 {
#ifdef DEBUG
                                    printf("HET group %s%d%s contacts Domain \
%d\n",
                                           pc->chain, pc->resnum, pc->insert,
                                           d->domainNumber);
                                    printf("%s%d%s.%s contacts %s%d%s.%s\n",
                                           pc->chain, pc->resnum,
                                           pc->insert, pc->atnam,
                                           pd->chain, pd->resnum,
                                           pd->insert, pd->atnam);
#endif
                                    /* Flag this HET residue as making contact*/
                                    if(d->nHetAntigen < MAXHETANTIGEN)
                                    {
                                       d->hetAntigen[d->nHetAntigen++] = r;
#ifndef DEBUG
                                       printf("Stored domain %d residue %s\n",
                                              d->domainNumber, r->resid);
#endif
                                    }
                                    goto lastatom;
                                 }  /* In range                               */
                              }  /* Step through atoms in this antibody residue*/
                           }  /* If it's a CDR residue */
                        }  /* Step through residues in this antibody domain*/
                     }  /* step through atoms in the non-water HET res  */
                  lastatom:
                     continue;
                  }  /* step through domains                            */
               }  /* not water                                          */
            }  /* Is a non-peptide HET                                  */
         }  /* step through residues in this non-HET chain              */
      }  /* Is a non-HET chain                                          */
   }  /* Step through chains                                            */
}

int CountResidueAtoms(PDBRESIDUE *res)
{
   PDB *p;
   int nAtoms = 0;
   
   
   /* Step through atoms in this residue                                */
   for(p=res->start; p!=res->stop; NEXT(p))
   {
      nAtoms++;
   }
   return(nAtoms);
}


BOOL IsNonPeptideHet(WHOLEPDB *wpdb, PDBRESIDUE *res)
{
   PDB *p;
   int hasBackbone = 0,
       nAtoms      = 0;
   BOOL isAllHet   = TRUE;
   STRINGLIST *s;
   
   
   /* Step through atoms in this residue                                */
   for(p=res->start; p!=res->stop; NEXT(p))
   {
      if(!strncmp(p->record_type, "HETATM", 6))
      {
         if(!strncmp(p->atnam, "N   ", 4) ||
            !strncmp(p->atnam, "CA  ", 4) ||
            !strncmp(p->atnam, "C   ", 4) ||
            !strncmp(p->atnam, "O   ", 4) ||
            !strncmp(p->atnam, "P   ", 4) ||
            !strncmp(p->atnam, "OP1 ", 4) ||
            !strncmp(p->atnam, "OP2 ", 4))
            hasBackbone++;
      }
      else
      {
         isAllHet = FALSE;
      }
      nAtoms++;
   }
   if((hasBackbone >= 3) || (!isAllHet) || (nAtoms < MINHETATOMS))
      return(FALSE);
   
   /* Check it isn't just an ion                                        */
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      /* Find HETNAM records in the header                              */
      if(!strncmp(s->string, "HETNAM", 6))
      {
         char hetnam[8], *chp;

         /* Grab this HETNAM record residue name                        */
         strncpy(hetnam, s->string+11, 3);
         hetnam[4] = '\0';
         KILLLEADSPACES(chp, hetnam);

         /* If it matches the residue in question                       */
         if(!strncmp(chp, res->start->resnam, 3))
         {
            /* If the record contains the word ' ION' it's an ion       */
            if(strstr(s->string, " ION"))
               return(FALSE);
         }
      }
   }

   return(TRUE);
   
}


BOOL IsStandardResidue(PDBRESIDUE *res)
{
   PDB        *p;
   int        hasBackbone = 0;
   BOOL       isAllHet    = TRUE;
   
   /* Step through atoms in this residue                                */
   for(p=res->start; p!=res->stop; NEXT(p))
   {
      /* As soon as we find an ATOM record it's standard                */
      if(!strncmp(p->record_type, "ATOM  ", 6))
         return(TRUE);
      
      if(!strncmp(p->record_type, "HETATM", 6))
      {
         if(!strncmp(p->atnam, "N   ", 4) ||
            !strncmp(p->atnam, "CA  ", 4) ||
            !strncmp(p->atnam, "C   ", 4) ||
            !strncmp(p->atnam, "O   ", 4) ||
            !strncmp(p->atnam, "P   ", 4) ||
            !strncmp(p->atnam, "OP1 ", 4) ||
            !strncmp(p->atnam, "OP2 ", 4))
            hasBackbone++;
      }
      else
      {
         isAllHet = FALSE;
      }
   }

   if((hasBackbone <= 2) || (isAllHet))
      return(FALSE);
   
   return(TRUE);
}



PDB *RelabelAntibodyChain(DOMAIN *domain, BOOL *lowerCaseLight,
                          BOOL *lowerCaseHeavy, char *remark950)
{
   PDB *p, *q;
   PDB *pdb = NULL;

   remark950[0] = '\0';
   
   if(domain)
   {
      /* 12.07.22 Support light chain and heavy chain dimers with upper 
         and lower case labels.
      */
      char chainLabel = domain->chainType;
      if(chainLabel == 'L')
      {
         if(*lowerCaseLight)
            chainLabel = 'l';
         *lowerCaseLight = TRUE;
      }
      else if(chainLabel == 'H')
      {
         if(*lowerCaseHeavy)
            chainLabel = 'h';
         *lowerCaseHeavy = TRUE;
      }
      
      domain->newAbChainLabel[0] = chainLabel;
      domain->newAbChainLabel[1] = '\0';

      sprintf(remark950, "REMARK 950 CHAIN %c     %c%6s\n",
              domain->chainType, chainLabel,
              domain->startRes->chain);

      for(p=domain->startRes; p!=domain->stopRes; NEXT(p))
      {
         if(pdb==NULL)
         {
            INIT(pdb, PDB);
            q=pdb;
         }
         else
         {
            ALLOCNEXT(q, PDB);
         }
         if(q==NULL)
         {
            FREELIST(pdb, PDB);
            return(NULL);
         }
         
         blCopyPDB(q, p);
         q->chain[0] = chainLabel;
         q->chain[1] = '\0';
      }
   }
   return(pdb);
}

   
PDB *RelabelAntigenChains(DOMAIN *domain, char *remark950)
{
   PDB *p, *q;
   PDB *pdb = NULL;
   int i;
   
   remark950[0] = '\0';

   /* Write antigen chains                                              */
   for(i=0; i<domain->nAntigenChains; i++)
   {
      PDBCHAIN *chain = domain->antigenChains[i];
      char     chainLabel[8],
               record[MAXBUFF];

      strcpy(chainLabel, chain->start->chain);
      if(CHAINMATCH(chainLabel, "L") || CHAINMATCH(chainLabel, "H"))
      {
         LOWER(chainLabel);
      }

      sprintf(record, "REMARK 950 CHAIN A%6s%6s\n",
              chainLabel, chain->start->chain);
      strcat(remark950, record);
      domain->newAgChainLabels[i][0] = chainLabel[0];
      domain->newAgChainLabels[i][1] = '\0';

      for(p=chain->start; p!=chain->stop; NEXT(p))
      {
         if(pdb==NULL)
         {
            INIT(pdb, PDB);
            q=pdb;
         }
         else
         {
            ALLOCNEXT(q, PDB);
         }
         if(q==NULL)
         {
            FREELIST(pdb, PDB);
            return(NULL);
         }
         
         blCopyPDB(q, p);
         strcpy(q->chain, chainLabel);
      }
   }
   return(pdb);
}

void WriteSeqres(FILE *fp, WHOLEPDB *wpdb, DOMAIN *domain)
{
   STRINGLIST *s;
   int        i;
   char       buffer[MAXBUFF];

   /* Print any MODRES records                                          */
   /************** TODO: need to fix the chain name          ************/
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "MODRES", 6))
      {
         fprintf(fp, "%s", s->string);
      }
   }
   
   /* Print SEQRES for this domain's chain                              */
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "SEQRES", 6))
      {
         if(s->string[11] == domain->startRes->chain[0])
         {
            strcpy(buffer, s->string);
            buffer[11] = domain->newAbChainLabel[0];
            fprintf(fp, "%s", buffer);
         }
      }
   }

   /* Print SEQRES for the partner domain's chain                       */
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "SEQRES", 6))
      {
         if(domain->pairedDomain != NULL)
         {
            if(s->string[11] == domain->pairedDomain->startRes->chain[0])
            {
               strcpy(buffer, s->string);
               buffer[11] = domain->pairedDomain->newAbChainLabel[0];
               fprintf(fp, "%s", buffer);
            }
         }
      }
   }
   
   /* Print SEQRES for the antigen chains                               */
   for(i=0; i<domain->nAntigenChains; i++)
   {
      PDBCHAIN *chain = domain->antigenChains[i];

      for(s=wpdb->header; s!=NULL; NEXT(s))
      {
         if(!strncmp(s->string, "SEQRES", 6))
         {
            if(s->string[11] == chain->chain[0])
            {
               strcpy(buffer, s->string);
               buffer[11] = domain->newAgChainLabels[i][0];
               fprintf(fp, "%s", buffer);
            }
         }
      }
   }
}


/* Note that we don't actually use the output from this!                */
char *blFixSequenceWholePDB(WHOLEPDB *wpdb, char **outChains,
                            BOOL ignoreSeqresForMissingChains,
                            BOOL upper, BOOL quiet, char *label)
{
   PDB    *pdb;
   MODRES *modres         = NULL;
   char   *seqresSequence = NULL,
          *atomSequence   = NULL,
          *fixedSequence  = NULL,
          **seqresChains  = NULL,
          **atomChains    = NULL;
   int    nAtomChains,
          len1;
   
   pdb = wpdb->pdb;
   
   if((seqresChains = (char **)blArray2D(sizeof(char),
                                         MAXCHAINS,
                                         blMAXCHAINLABEL))==NULL)
   {
      fprintf(stderr,"Error: No memory for seqresChains array\n");
      return(NULL);
   }
   
   /* Read MODRES and SEQRES records                                    */
   modres         = blGetModresWholePDB(wpdb);
   seqresSequence = blGetSeqresAsStringWholePDB(wpdb,
                                                seqresChains,
                                                modres, TRUE);
   
   /* Get list of chains from the PDB linked list                       */
   if((atomChains = blGetPDBChainLabels(pdb, &nAtomChains))
      == NULL)
   {
      fprintf(stderr,"Error: No memory for atom chain labels\n");
      return(NULL);
   }
   
   /* Convert PDB linked list to a sequence                             */
   if((atomSequence = blPDB2SeqX(pdb))==NULL)
   {
      fprintf(stderr,"Error: No memory for sequence data\n");
      return(NULL);
   }
   /* Append a * since blPDB2Seq() doesn't do this; note that
      this will have to change if we fix blPDB2Seq() in future
   */
   len1 = strlen(atomSequence);
   if((atomSequence=(char *)realloc((void *)atomSequence,
                                    (len1+2)*sizeof(char)))==NULL)
   {
      fprintf(stderr,"Error: No memory to expand sequence data\n");
      return(NULL);
   }
   strcat(atomSequence,"*");
   
#ifdef DEBUG
   fprintf(stderr,"\nSEQRES sequence:\n");
   fprintf(stderr,"%s\n", seqresSequence);
   fprintf(stderr,"\nATOM sequence:\n");
   fprintf(stderr,"%s\n", atomSequence);
#endif
   
   /* Fiddle with sequences to combine information from SEQRES
      and ATOM records
   */
   if((fixedSequence = blFixSequence(seqresSequence,atomSequence,
                                     seqresChains,atomChains,outChains,
                                     ignoreSeqresForMissingChains,
                                     nAtomChains,upper,quiet,label))
      ==NULL)
      return(NULL);
   
#ifdef DEBUG
   fprintf(stderr, "Fixed sequence:\n");
   fprintf(stderr, "%s\n", fixedSequence);
#endif

   /* Free allocated memory                                             */
   if(seqresChains!=NULL)
      blFreeArray2D(seqresChains, MAXCHAINS, blMAXCHAINLABEL);
   if(atomChains!=NULL)
      blFreeArray2D(atomChains, nAtomChains, blMAXCHAINLABEL);

   if(modres!=NULL)
      FREELIST(modres, MODRES);

   FREE(seqresSequence);
   FREE(atomSequence);
   
   return(fixedSequence);
}


int FindLastAlignmentPosition(char *refAln)
{
   int pos;

   for(pos = strlen(refAln) - 1;
       ((pos >= 0) && (refAln[pos] == '-'));
       pos--);
   return(pos);
}

BOOL DomainSequenceMatchesChainSequence(DOMAIN *domain, PDBCHAIN *chain)
{
   if(domain != NULL)
   {
      PDBRESIDUE *r;
      char       chainSeq[MAXSEQ+1],
                 alignChainSeq[MAXSEQ+1],
                 alignDomSeq[MAXSEQ+1];
      int        nRes     = 0,
                 nAligned = 0,
                 nMatched = 0,
                 alignLen,
                 i;
      

      /* Assemble the chain sequence                                    */
      for(r=chain->residues; r!=NULL; NEXT(r))
      {
         chainSeq[nRes++] = blThrone(r->resnam);
      }
      chainSeq[nRes] = '\0';

      /* Align the sequences                                            */
      blAffinealign(chainSeq, strlen(chainSeq),
                    domain->domSeq, strlen(domain->domSeq),
                    FALSE,          /* verbose                  */
                    TRUE,           /* identity                 */
                    2,              /* penalty                  */
                    1,              /* extension                */
                    alignChainSeq,
                    alignDomSeq,
                    &alignLen);
      /* Score the matched residues                                     */
      for(i=0; i<alignLen; i++)
      {
         if((alignChainSeq[i] != '-') &&
            (alignDomSeq[i]   != '-'))
         {
            nAligned++;
            if(alignChainSeq[i] == alignDomSeq[i])
            {
               nMatched++;
            }
         }
      }

      if(((REAL)nMatched/(REAL)nAligned) >= 0.98)
      {
         return(TRUE);
      }
   }
   
   return(FALSE);
}
