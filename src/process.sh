input=$1
absplit=${HOME}/git/absplit/src/absplit
numberabpdb=${HOME}/git/absplit/src/numberabpdb.pl
combinefaa=${HOME}/git/absplit/src/combinefaa.pl
getfooter=${HOME}/git/absplit/src/getHETAndFooterRecords.pl

$absplit $input

for file in *.pdb
do
    echo $file
    pdbrepair -t $file > `basename $file .pdb`.fix
done

for file in *.fix
do
    tmpfaa=`basename $file .fix`.tmp
    faa=`basename $file .fix`.faa
    id=`basename $file .fix`
    id=`echo $id | sed 's/^pdb//'`
    pdbgetchain L,H $file | pdb2pir -s -i -c -f $file > $tmpfaa
    $combinefaa -l=$id $tmpfaa > $faa
done

for file in *.fix
do
    footer=`basename $file .fix`.foot
    $getfooter $file > $footer
    $numberabpdb -k $file `basename $file .fix`.kab
    cat $footer >> `basename $file .fix`.kab
    $numberabpdb -c $file `basename $file .fix`.cho
    cat $footer >> `basename $file .fix`.cho
    $numberabpdb -m $file `basename $file .fix`.mar
    cat $footer >> `basename $file .fix`.mar
    rm $footer
done

rm *.fix
rm *.pdb
rm *.tmp

# pdbrepair -t pdb7mfa_0.num | pdbdummystrip | pdb2pir -s 
