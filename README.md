# absplit
## (c) 2021 UCL, Andrew C.R. Martin

Code to take an antibody PDB file and split it into separate Fvs

## Known issues

1. The code needs to retain the PDB header for downstream steps
2. Code needs to apply symmetry operators (e.g. 1b0w, 2rhe)
   - decided to do this first with pdbsymm
3. Some truncated sequences won't score highly enough to be flagged
   as antibodies - just score over the aligned region.
4. Heavy chain only not being done properly (at least in header)
   e.g. 1t2j_0, 1shm_1
5. Bence-Jones dimers are labelling both chains as L - should be L
   and l e.g. 1rei
   

## Algorithm

1. Read the PDB file and split into chains
2. For each chain, scan with a library of VH and VL domains to locate
these. This is done by alignment and the standard interface positions
are recorded.
3. For each VH/VL domain, find the CoG
4. For each VH/VL domain measure the CofG distance to all others. If
within a cutoff, then this is a potential VH/VL pair.
5. For each potential VH/VL pair (note that these can be VL/VL or
VH/VL and may be within the same chain) check the standard interface
residues (assigned from the alignment) and ensure that they are close
together
6. Add back HETATMs (other than water)
7. We now have assigned VH/VL pairs, so check those against *other*
chains to look for antibody/antigen interactions. Note that an
'antigen' may be a different antibody.
8. Now check the HETATMs for HET antigens
9. Output the VH/VL domains together with any contacted antigen

