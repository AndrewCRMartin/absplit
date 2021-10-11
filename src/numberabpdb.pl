#!/usr/bin/perl

use strict;

my $inFile  = shift(@ARGV);
my $outFile = shift(@ARGV);

my $header = `grep REMARK $inFile`;
my $lightChain = `pdbgetchain L $inFile | egrep '^(ATOM|HETATM)'`;
my $heavyChain = `pdbgetchain H $inFile | egrep '^(ATOM|HETATM)'`;
my $antigen    = `egrep '^(ATOM|HETATM)' $inFile | grep -v ' L ' | grep -v ' H '`;

my $fileLH  = "/var/tmp/numberpdb_LH_$$"  . '_' . time();
my $fileNum = "/var/tmp/numberpdb_Num_$$" . '_' . time();
my $outTemp = "/var/tmp/numberpdb_Out_$$" . '_' . time();
if(open(my $fp, '>', $fileLH))
{
    print($fp $lightChain);
    print($fp $heavyChain);
    close($fp);
}
else
{
    die "Can't write $fileLH";
}

`pdbabnum $fileLH | egrep -v '^(MASTER|END)' > $fileNum`;
WriteToFile($outTemp, $header, 0);
`cat $fileNum >> $outTemp`;
WriteToFile($outTemp, $antigen, 1);
`pdbrenum -d $outTemp $outFile`;

unlink $fileLH;
unlink $fileNum;
unlink $outTemp;

sub WriteToFile
{
    my($fileName, $content, $append) = @_;
    my $type = ($append)?'>>':'>';
    if(open(my $fp, $type, $fileName))
    {
        print($fp $content);
        close($fp);
    }
    else
    {
        die "Can't write to $fileName";
    }
}


