input=$1
absplit=${HOME}/git/absplit/src/absplit
numberabpdb=${HOME}/git/absplit/src/numberabpdb.pl
combinefaa=${HOME}/git/absplit/src/combinefaa.pl
getfooter=${HOME}/git/absplit/src/getHETAndFooterRecords.pl

# Split the file into component Fvs and antigens
$absplit $input

stem=`basename $input .ent`
stem=`basename $stem  .pdb`

# Repair the PDB files adding missing residues and trimming the SEQRES
for file in ${stem}_*.pdb
do
    echo $file
    pdbrepair -t $file > `basename $file .pdb`.fix
done

# Extract the sequence for L and H and combine into one FASTA sequence
for file in ${stem}_*.fix
do
    tmpfaa=`basename $file .fix`.tmp
    faa=`basename $file .fix`.faa
    id=`basename $file .fix`
    id=`echo $id | sed 's/^pdb//'`
    pdbgetchain L,H $file | pdb2pir -s -i -c -f $file > $tmpfaa
    $combinefaa -l=$id $tmpfaa > $faa
done

# Extract any HETATM footer, number the antibodies and add back the footer.
# Finally renumber the atoms to reconstruct the MASTER and CONECT records
for file in ${stem}_*.fix
do
    footer=`basename $file .fix`.foot
    $getfooter $file > $footer

    $numberabpdb -k $file `basename $file .fix`.tmp
    cat $footer >> `basename $file .fix`.tmp
    pdbrenum -d `basename $file .fix`.tmp > `basename $file .fix`.kab
    
    $numberabpdb -c $file `basename $file .fix`.tmp
    cat $footer >> `basename $file .fix`.tmp
    pdbrenum -d `basename $file .fix`.tmp > `basename $file .fix`.cho
    
    $numberabpdb -m $file `basename $file .fix`.tmp
    cat $footer >> `basename $file .fix`.tmp
    pdbrenum -d `basename $file .fix`.tmp > `basename $file .fix`.mar
    rm $footer
done

# Cleanup
rm ${stem}_*.fix
rm ${stem}_*.pdb
rm ${stem}_*.tmp

# pdbrepair -t pdb7mfa_0.num | pdbdummystrip | pdb2pir -s 
