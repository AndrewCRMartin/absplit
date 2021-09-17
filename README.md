# absplit
## (c) 2021 UCL, Andrew C.R. Martin

Code to take an antibody PDB file and split it into separate Fvs

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

