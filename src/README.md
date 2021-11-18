Currently fails on
../../AbDb2/samples/pdb6bpe.ent
as it doesn't recognize chain I as an antibody (it has an NTer truncation)

The code in the UseSeqres branch is a start to fix this and reads the
sequence from the SEQRES records. However it needs to sort out the
correct offsets into the PDB linked list (assumed to be the same as
the sequence data)
Reverted the change to use SEQRES data

pdbabnum needs to take a flag to stop it changing the chain label


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



