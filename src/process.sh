input=$1
absplit=${HOME}/git/absplit/bin/absplit
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
#    pdbgetchain L,l,H,h $file | pdb2pir -s -i -c -f $file > $tmpfaa
    pdbgetchain L,l,H,h $file | pdb2pir -s -i -c -f > $tmpfaa
    $combinefaa -l=$id $tmpfaa > $faa
done

# Extract any HETATM footer, number the antibodies and add back the footer.
# Finally renumber the atoms to reconstruct the MASTER and CONECT records
badFile=0
for file in ${stem}_*.fix
do
#    footer=`basename $file .fix`.foot
#    $getfooter $file > $footer

    tmpfile=`basename $file .fix`.tmp

    $numberabpdb -k $file $tmpfile
    line1=`head -1 $tmpfile | awk '{print $1}'`
    if [ "X$line1" == "XMASTER" ]; then
        badFile=1
    else
#       cat $footer >> `basename $file .fix`.tmp
       pdbrenum -d `basename $file .fix`.tmp > `basename $file .fix`.kab
    fi
    
    $numberabpdb -c $file $tmpfile
    line1=`head -1 $tmpfile | awk '{print $1}'`
    if [ "X$line1" == "XMASTER" ]; then
        badFile=1
    else
#        cat $footer >> `basename $file .fix`.tmp
        pdbrenum -d `basename $file .fix`.tmp > `basename $file .fix`.cho
    fi
    
    $numberabpdb -m $file $tmpfile
    line1=`head -1 $tmpfile | awk '{print $1}'`
    if [ "X$line1" == "XMASTER" ]; then
        badFile=1
    else
#        cat $footer >> `basename $file .fix`.tmp
        pdbrenum -d `basename $file .fix`.tmp > `basename $file .fix`.mar
    fi
#    rm $footer

    if [ $badFile == 1 ]; then
        mv $file `basename $file .fix`.bad
    else
        rm $file
    fi
done

# Cleanup
#rm ${stem}_*.fix
rm ${stem}_*.pdb
rm ${stem}_*.tmp

# pdbrepair -t pdb7mfa_0.num | pdbdummystrip | pdb2pir -s 
