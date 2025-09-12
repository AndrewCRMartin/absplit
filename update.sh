sacslist=/serv/www/html_abybank/sacs/list/antibodies.txt
absplitdir=${HOME}/git/absplit/
builddirtop=/data/abdbbuild
webdir=/serv/www/html_abybank/abdb/snapshots/

datestamp=`date +%4Y%m%d`
abdbdir=abdb_${datestamp}
builddir="${builddirtop}/${abdbdir}"
mkdir -p $builddir
if [ ! -d $builddir ]; then
   echo "Couldn't create build directory: $builddir"
   exit
fi

cd $builddir
pwd=`pwd`
if [ "X$pwd" != "X$builddir" ]; then
    echo "Unable to change to build directory"
else
    nice -10 $absplitdir/src/processall.sh $sacslist &> update-${datestamp}.log
    cd $builddirtop
    zip -r ${abdbdir}.zip $abdbdir &>/dev/null
    mv ${abdbdir}.zip $webdir
    rm -rf $abdbdir
fi



