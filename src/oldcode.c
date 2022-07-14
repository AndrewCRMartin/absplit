#ifdef OLD
/* seqPos numbered from 1 */
int IsKeyResidue(int seqPos, int *refKeys,
                 char *seqAln, char *refAln)
{
   int i, count, alnPos=0, refPos;
   
   /* Correct seqPos to the position in the alignment (numbered from 1) */
   for(i=0, count=0; i<strlen(seqAln); i++)
   {
      if((seqAln[i] != '-') || (seqAln[i] != 'X'))
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
      if((refAln[i] != '-') && (refAln[i] != 'X'))
         refPos++;
   }

   /* Now see if this position is in the list of key reference positions*/
   for(i=0; refKeys[i] >=0; i++)
   {
#ifndef DEBUG
      printf("seqPos %d refPos %d refKeys[%d] %d\n", seqPos, refPos, i, refKeys[i]);
#endif
      if(refPos == refKeys[i])
         return(1);
   }
   
   return(0);
}
#endif

/************************************************************************/
void MaskAndAssignDomain(char *seq, PDBCHAIN *chain, char *fastaHeader,
                         char *seqAln, char *refAln, DOMAIN **pDomains)
{
   int    seqPos      = 0,
          alnPos      = 0,
#ifdef OLD
          firstAlnPos = 0,
          lastAlnPos  = 0;
#endif
          domSeqPos   = 0;
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


#ifdef DEBUG_SET_CDR
   printf("SEQ      : %s\n", seqAln);
   printf("REF      : %s\n", refAln);
#endif

   SetChainAsLightOrHeavy(d, fastaHeader);

#ifdef OLD
   lastAlnPos = FindLastAlignmentPosition(refAln);

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
#endif

   /* Mask the sequence */
   for(seqPos=0, alnPos=0;
       seqPos<strlen(seq) && alnPos<strlen(seqAln);
       alnPos++)
   {
      /* If this is an aligned position */
      if((seqAln[alnPos] != '-') &&
         (seqAln[alnPos] != 'X') &&
         (refAln[alnPos] != '-'))
      {
         d->domSeq[domSeqPos++] = seq[seqPos];
         seq[seqPos] = 'X';
         if(d->startSeqRes < 0)
            d->startSeqRes = seqPos;
         if(seqPos > d->lastSeqRes)
            d->lastSeqRes = seqPos;
      }
      if(seqAln[alnPos] != '-')
         seqPos++;
   }

   d->domSeq[domSeqPos] = '\0';

   SetIFResidues(d,          fastaHeader, seqAln, refAln);
   SetCDRResidues(d,         fastaHeader, seqAln, refAln);
#ifndef DEBUG
   {
      int i;
      
      printf("IF Residue numbers: \n");
      for(i=0;i<d->nInterface; i++)
         printf("%d ", d->interface[i]);
      printf("\n");
      printf("CDR Residue numbers: \n");
      for(i=0;i<d->nCDRRes; i++)
         printf("%d ", d->CDRRes[i]);
      printf("\n");
   }
#endif
   
   SetDomainBoundaries(d);
}

