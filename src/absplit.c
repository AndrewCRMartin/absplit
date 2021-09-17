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
#define PROGNAME "absplit"
#define VERSION  0.1
#define MAXBUFF  240
#define HUGEBUFF 10000
#define MAXSEQ   10000

/************************************************************************/
/* Globals
*/
BOOL gVerbose = TRUE;

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile);
void UsageDie(void);
BOOL ProcessFile(WHOLEPDB *wpdb, char *infile, FILE *dataFp);
void FindVHVLDomains(PDBCHAIN *chain, FILE *dataFp);
void GetSequenceForChain(PDBCHAIN *chain, char *sequence);
void ExePathName(char *str, BOOL pathonly);
BOOL CheckAndMask(char *sequence, FILE *dataFp);
FILE *OpenSequenceDataFile(void);
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2);
int CalcShortSeqLen(char *align1, char *align2);
void Mask(char *seq, char *aln1, char *aln2);







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
   char filestem[MAXBUFF];
   PDBSTRUCT *pdbs;
   
   GetFilestem(infile, filestem);
   if((pdbs = blAllocPDBStructure(wpdb->pdb))!=NULL)
   {
      PDBCHAIN *chain;
      
      for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
      {
#ifdef DEBUG
         printf("Ptr: %ld Type: %s\n", (ULONG)chain->start,
                chain->start->record_type);
#endif
         
         if(!strncmp(chain->start->record_type, "ATOM", 4))
         {
            printf("***Handling chain: %s\n", chain->chain);
            FindVHVLDomains(chain, dataFp);
         }
         
      }
      
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
   for(r=chain->residues; r!=NULL; NEXT(r))
   {
      sequence[i++] = blThrone(r->resnam);
   }
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
*/
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2)
{
   int  score;
   int  alignLen;
   int  shortSeqLen = MIN(strlen(theSeq), strlen(seq));
   
   score = blAlign(theSeq, strlen(theSeq),
                   seq, strlen(seq),
                   FALSE, /* verbose  */
                   TRUE,  /* identity */
                   5,     /* penalty  */
                   align1,
                   align2,
                   &alignLen);
   align1[alignLen] = align2[alignLen] = '\0';

   shortSeqLen = CalcShortSeqLen(align1, align2);

#ifdef DEBUG   
   fprintf(stderr, "\n>>>%s\n", align1);
   fprintf(stderr, ">>>%s %d\n", align2, shortSeqLen);
#endif
   
   return((REAL)score / (REAL)shortSeqLen);
}


/************************************************************************/
int CalcShortSeqLen(char *align1, char *align2)
{
   int start, stop, i, len1=0, len2=0,
      alnLen = strlen(align1);
   

   /* Find where the alignment starts                                   */
   for(start=0; start<alnLen; start++)
   {
      if((align1[start] != '-') && (align2[start] != '-'))
         break;
   }
   /* Find where the alignment stops                                    */
   for(stop=alnLen-1; stop>start; stop--)
   {
      if((align1[stop] != '-') && (align2[stop] != '-'))
         break;
   }
   /* Step between start and stop and calculate the number of residues 
      in each sequence
   */
   for(i=start; i<=stop; i++)
   {
      if(align1[i] != '-')
         len1++;
      if(align2[i] != '-')
         len2++;
   }

   return(MIN(len1, len2));
}
            
            
            
void FindVHVLDomains(PDBCHAIN *chain, FILE *dataFp)
{
   char sequence[MAXSEQ];
   GetSequenceForChain(chain, sequence);
#ifdef DEBUG
   printf("Chain: %s Sequence: %s\n", chain->chain, sequence);
#endif
   while(TRUE)
   {
      if(!CheckAndMask(sequence, dataFp)) break;
   }
}

BOOL CheckAndMask(char *sequence, FILE *dbFp)
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
   if(maxScore > 0.5)
   {
      fprintf(stderr, "Best match: %s *** %.4f\n",
              bestMatch, maxScore);
      fprintf(stderr, "%s\n", bestAlign1);
      fprintf(stderr, "%s\n", bestAlign2);
      Mask(sequence, bestAlign1, bestAlign2);
      found = TRUE;
   }
   
   return(found);
}

void Mask(char *seq, char *aln1, char *aln2)
{
   int i, j;
   for(i=0, j=0; i<strlen(seq); i++)
   {
      while(aln1[j] == '-')  /* Skip insertions */
      {
         j++;
      }
      if(aln2[j] != '-')
      {
         seq[i] = 'X';
      }
      j++;
   }
}


