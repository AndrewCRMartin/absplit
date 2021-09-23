Currently fails on

../../AbDb2/samples/pdb6bpe.ent
as it doesn't recognize chain I as an antibody (it has an NTer truncation)

The code in the UseSeqres branch is a start to fix this and reads the
sequence from the SEQRES records. However it needs to sort out the
correct offsets into the PDB linked list (assumed to be the same as
the sequence data)


