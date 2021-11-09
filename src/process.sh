input=$1
absplit=${HOME}/git/absplit/src/absplit
numberabpdb=${HOME}/git/absplit/src/numberabpdb.pl

$absplit /serv/data/pdb/pdb7mfa.ent

for file in *.pdb
do
    echo $file
    pdbrepair -t $file > `basename $file .pdb`.fix
done

for file in *.fix
do
    pdb2pir -s -c $file > `basename $file .fix`.pir
done

for file in *.fix
do
    $numberabpdb -k $file `basename $file .fix`.kab
    $numberabpdb -k $file `basename $file .fix`.cho
    $numberabpdb -k $file `basename $file .fix`.mar
done

rm *.fix
rm *.pdb

# pdbrepair -t pdb7mfa_0.num | pdbdummystrip | pdb2pir -s 
