#include <string.h>
#include <stdio.h>

int main(int argc, char **argv);
int IsKeyResidue(int seqPos, int *refKeys, char *seqAln, char *refAln);
int RealSeqLen(char *seq);


int main(int argc, char **argv)
{
   char *seqAln = "-VQLQTIGRNLVGPSD";
   char *refAln = "AVQLQRI-RTK--PNN";
   int  refKeys[] = {3,4,5,8,9,10,11,-1};
   int  seqPos, refPos;
   int seqLen = RealSeqLen(seqAln);

   for(seqPos=1; seqPos<=seqLen; seqPos++)
   {
      refPos = IsKeyResidue(seqPos, refKeys, seqAln, refAln);
      if(refPos)
         printf("Seq %d Ref %d\n", seqPos, refPos);
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

