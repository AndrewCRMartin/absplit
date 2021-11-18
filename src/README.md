6bpe
----

- doesn't recognize chain I as an antibody (it has an NTer truncation)

The code in the UseSeqres branch is a start to fix this and reads the
sequence from the SEQRES records. However it needs to sort out the
correct offsets into the PDB linked list (assumed to be the same as
the sequence data)

Reverted the change to use SEQRES data as too much else breaks. The
answer may be to use SEQRES to find the Ab region and then redo the
alignment with the ATOM sequence.

Scoring is now based on a percentage of the shorter sequence which
largely fixes this.

- Note that chains H and I are still not paired presumably because of
  all the missing residues and the interface residues not being correctly
  identified.

- Note that chain H finds chain I as an 'antigen' with 14 contacts and
  finds G with just 9 contacts. This is probably because of the
  problem with not correcting the CDR residue numbers based on the
  alignment


1a6v
----

- 3 VH/VL pairs
- a HET chain residue for each
- first pair contacts a chain from the second pair
- second pair contacts a chain from the third pair

5jor
----

A good example where chains B/C are indicated as interacting with
chain L as an antigen - it really is just touching the edge and
shouldn't be listed as an antigen

This is also an example where the CDR residue numbers need to be
corrected by the alignment



pdbabnum needs to take a flag to stop it changing the chain label

