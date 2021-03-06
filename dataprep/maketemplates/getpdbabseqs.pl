#!/usr/bin/perl -s
use strict;
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../..");
use config;


my $parentDir  = abs_path("$FindBin::Bin/../..");
my $configFile = "$parentDir/absplit.conf";
my $dataDir    = "$parentDir/data";

my %config  = config::ReadConfig($configFile);
my $abdir   = $config{'abpdbdir'};
my $bindir  = $config{'bindir'};
my $datasub = $config{'datasub'};
my $allseq  = "$dataDir/allabs.faa";
my $repseq  = "$dataDir/cdhit.faa";

if((! -e "$repseq") || defined($::f) || defined($::force))
{
    `mkdir -p $dataDir`;
    unlink $allseq;
    print STDERR "Getting list of antibody PDB files...";
    my @pdbfiles = GetFileList($abdir, '.pdb');
    print STDERR "done\n";

    print STDERR "Extracting sequence data and concatenating files...";
    foreach my $file (@pdbfiles)
    {
        my $pdbcode = $file;
        $pdbcode =~ s/.*\///;
        $pdbcode =~ s/\..*?$//;
        `pdbgetchain L,H ${abdir}/$file | pdb2pir -c -f -l ${pdbcode}_ >>$allseq`;
        print STDERR '.';
    }
    print STDERR "done\n";

    # Grab and compile CD-HIT if not there
    if( ! -d "cdhit" )
    {
        print STDERR "Building cd-hit\n";
        `git clone git\@github.com:weizhongli/cdhit.git`;
        `(cd cdhit; make)`;
        print STDERR "done\n";
    }
    
    print STDERR "Running CD-HIT...";
    `cdhit/cd-hit -c 0.6 -n 3 -i $allseq -o $repseq`;
    print STDERR "done\n";
    
    #SplitFastaFiles($repseq, $tplDir);
    #print STDERR "Antibody template sequences are in $tplDir\n";
    unlink $allseq;
    unlink "${repseq}.clstr";
}
else
{
    print STDERR "Non-redundant sequences already installed\n";
}


#*************************************************************************
sub SplitFastaFiles
{
    my($repseq, $tplDir) = @_;

    if(open(my $fpin, '<', $repseq))
    {
        my $seq = '';
        while(<$fpin>)
        {
            if(/\>/)
            {
                if($seq ne '')
                {
                    PrintSeq($seq, $tplDir);
                    $seq = '';
                }
            }
            $seq .= $_;
        }
        
        if($seq ne '')
        {
            PrintSeq($seq, $tplDir);
            $seq = '';
        }
    }
}

sub PrintSeq
{
    my($seq, $tplDir) = @_;
    my @lines = split(/\n/, $seq);
    my $filename = $lines[0];
    $filename =~ s/\>//;
    $filename = "${tplDir}/${filename}.faa";
    if(open(my $fpOut, '>', $filename))
    {
        print $fpOut $seq;
        close $fpOut;
    }
    else
    {
        print STDERR "Error: Unable to write $filename\n";
    }
}



#*************************************************************************
sub GetFileList
{
    my($dir, $ext) = @_;
    my @files = ();
    if(opendir(my $dirFp, $dir))
    {
        @files = grep(/$ext/, readdir($dirFp));
        closedir($dirFp);
    }
    return(@files);
}
