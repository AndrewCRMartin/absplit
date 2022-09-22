#!/usr/bin/perl -s

use strict;

my $inFile  = shift(@ARGV);
my $outFile = shift(@ARGV);

my $scheme = '-c';
$scheme = '-k' if(defined($::k));
$scheme = '-m' if(defined($::m));

my $header = `egrep '(REMARK|MODRES)' $inFile`;
my $lightSeqres   = `fgrep SEQRES $inFile | egrep '( L | l )'`;
my $heavySeqres   = `fgrep SEQRES $inFile | egrep '( H | h )'`;
my $antigenSeqres = `fgrep SEQRES $inFile | egrep -v -i '( H | L )'`;
my $lightChain   = `pdbgetchain L,l $inFile | egrep '^(ATOM|HETATM)'`;
my $heavyChain   = `pdbgetchain H,h $inFile | egrep '^(ATOM|HETATM)'`;
my $antigen      = `egrep '^(ATOM|HETATM)' $inFile | grep -v -i ' L ' | grep -v -i ' H '`;

my $fileLH  = "/var/tmp/numberpdb_LH_$$"  . '_' . time();
my $fileNum = "/var/tmp/numberpdb_Num_$$" . '_' . time();
my $outTemp = "/var/tmp/numberpdb_Out_$$" . '_' . time();

# Grab just light and heavy chains into a temporary file
if(open(my $fp, '>', $fileLH))
{
    print($fp $lightSeqres);
    print($fp $heavySeqres);
    print($fp $lightChain);
    print($fp $heavyChain);
    close($fp);
}
else
{
    die "Can't write $fileLH";
}

my $errFile = $inFile;
$errFile =~ s/^.*\///;
$errFile .= $scheme . ".err";

# Apply numbering to the temp file and save in another temp file
`pdbabnum $scheme $fileLH 2> $errFile | egrep -v '^(MASTER|END|SEQRES)' > $fileNum`;

# Write the header to the final output tempfile
WriteToFile($outTemp, $header, 0);
# Add the SEQRES data
WriteToFile($outTemp, $lightSeqres, 1);
WriteToFile($outTemp, $heavySeqres, 1);
WriteToFile($outTemp, $antigenSeqres, 1);

# Add the numbered coordinates
`cat $fileNum >> $outTemp`;
# Add the antigen
WriteToFile($outTemp, $antigen, 1);
# Renumber the atoms
#`pdbdummystrip $outTemp | pdbrenum -d > $outFile`;
`pdbdummystrip $outTemp > $outFile`;
FixChainLabels($outFile, $outTemp);
`pdbrenum -d $outFile > $outTemp`;
`pdbconect $outTemp $outFile`;

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

sub FixChainLabels
{
    my($pdbFile, $tmpFile) = @_;

    my $text = `grep 'REMARK 950 CHAIN ' $pdbFile | awk '{print \$5}'`;
    my @chains = split(/\s+/, $text);
    if(($chains[0] eq 'L') && ($chains[1] eq 'l'))
    {
        FixSecondChain($pdbFile, $tmpFile, 'L','l');
    }
    elsif(($chains[0] eq 'H') && ($chains[1] eq 'h'))
    {
        FixSecondChain($pdbFile, $tmpFile, 'H','h');
    }
}


# Relabel second L or H chain as l or h
sub FixSecondChain
{
    my($pdbFile, $tmpFile, $oldChain, $newChain) = @_;
    my $maxResnum = -100;
    my $relabel   = 0;
    
    if(open(my $in, '<', $pdbFile))
    {
        if(open(my $out, '>', $tmpFile))
        {
            while(<$in>)
            {
                if(/^ATOM/)
                {
                    my $atomLine = $_;
                    my $chain    = substr($_, 21, 1);
                    my $resnum   = substr($_, 22, 4);
                    if($chain eq $oldChain)
                    {
                        if(!$relabel)
                        {
                            if($resnum > $maxResnum)
                            {
                                $maxResnum = $resnum;
                            }
                            elsif($resnum < $maxResnum)
                            {
                                $relabel = 1;
                                print $out "TER   \n";
                            }
                        }

                        if($relabel)
                        {
                            $atomLine = substr($_, 0, 21) . $newChain . substr($_, 22);
                        }
                        print $out $atomLine;
                    }
                    else
                    {
                        print $out $_;
                    }
                }
                else
                {
                    print $out $_;
                }
            }
            close $out;
        }
        else
        {
            die "Can't write file: $tmpFile\n";
        }
        close $in;
        `cp $tmpFile $pdbFile`;    
    }
    else
    {
        die "Can't read file: $pdbFile\n";
    }   
}
