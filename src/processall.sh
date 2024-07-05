ablist=$1
process=${HOME}/git/absplit/src/process.sh
nr=${HOME}/git/absplit/src/nr.pl
pdbdir=/serv/data/pdb

for ab in `awk -F_ '{print $1}' $ablist | sort -u`
do
    pdbfile="$pdbdir/pdb${ab}.ent"
    echo $pdbfile
    $process $pdbfile
done

$nr `pwd` > AbClusters.txt

