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

#include "absplit.h"

/************************************************************************/
/* Defines and macros
*/
#define PROGNAME "absplit"
#define VERSION  0.1
#define MAXBUFF  240
#define MAXSEQ   10000

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile);
void UsageDie(void);
BOOL ProcessFile(WHOLEPDB *wpdb, char *infile, FILE *dataFp);
void FindVHVLDomains(PDBCHAIN *chain, FILE *dataFp);
void GetSequenceForChain(PDBCHAIN *chain, char *sequence);
void ExePathName(char *str, BOOL pathonly);
void CheckAndMask(char *sequence, FILE *dataFp);
FILE *OpenSequenceDataFile(void);





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
#ifdef DEBUG
      int count = 1;
#endif
      
      for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
      {
#ifdef DEBUG
         printf("Ptr: %ld Type: %s\n", (ULONG)chain->start,
                chain->start->record_type);
#endif
         
         if(!strncmp(chain->start->record_type, "ATOM", 4))
         {
#ifdef DEBUG
            printf("Chain: %s %d\n", chain->chain, count++);
#endif
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

            
            
            
void FindVHVLDomains(PDBCHAIN *chain, FILE *dataFp)
{
   char sequence[MAXSEQ];
   GetSequenceForChain(chain, sequence);
#ifdef DEBUG
   printf("Chain: %s Sequence: %s\n", chain->chain, sequence);
#endif
   CheckAndMask(sequence, dataFp);
   
   
}

void CheckAndMask(char *sequence, FILE *dataFp)
{
}



