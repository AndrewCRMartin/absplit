#!/usr/bin/perl -s

use strict;

my $inFile  = shift(@ARGV);
my $outFile = shift(@ARGV);

my $scheme = '-c';
$scheme = '-k' if(defined($::k));
$scheme = '-m' if(defined($::m));

my $header = `egrep '(REMARK|SEQRES)' $inFile`;
my $lightChain = `pdbgetchain L $inFile | egrep '^(ATOM|HETATM)'`;
my $heavyChain = `pdbgetchain H $inFile | egrep '^(ATOM|HETATM)'`;
my $antigen    = `egrep '^(ATOM|HETATM)' $inFile | grep -v ' L ' | grep -v ' H '`;

my $fileLH  = "/var/tmp/numberpdb_LH_$$"  . '_' . time();
my $fileNum = "/var/tmp/numberpdb_Num_$$" . '_' . time();
my $outTemp = "/var/tmp/numberpdb_Out_$$" . '_' . time();

# Grab just light and heavy chains into a temporary file
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

# Apply numbering to the temp file and save in another temp file
`pdbabnum $scheme $fileLH | egrep -v '^(MASTER|END)' > $fileNum`;

# Write the header to the final output tempfile
WriteToFile($outTemp, $header, 0);
# Add the numbered coordinates
`cat $fileNum >> $outTemp`;
# Add the antigen
WriteToFile($outTemp, $antigen, 1);
# Renumber the atoms
`pdbdummystrip $outTemp | pdbrenum -d > $outFile`;

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


