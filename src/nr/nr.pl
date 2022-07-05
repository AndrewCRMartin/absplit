#!/usr/bin/perl

use strict;

$::gCDHit = '../../dataprep/cdhit/cd-hit';

my $fastaDir = shift @ARGV;
my $tmpDir   = MakeTempDir();
Die("Cannot create temporary directory $tmpDir") if($tmpDir eq '');
my $faaFile = "$tmpDir/all.faa";

CombineFAA($fastaDir, $faaFile);
print "Temp dir: $tmpDir\n";
Cluster($faaFile, $tmpDir);


sub Cluster
{
    my($faaFile, $tmpDir) = @_;
    `$::gCDHit -c 1.0 -i $faaFile -o $tmpDir/clusters.faa`;
    my $clsFile = "$tmpDir/clusters.faa.clstr";
    my %clusters = ReadClusters($clsFile);
}




sub CombineFAA
{
    my($faaDir, $outFile) = @_;
    
    if(opendir(my $fpDir, $faaDir))
    {
        if(open(my $fpOut, '>', $outFile))
        {
            my @files = grep /\.faa/, readdir($fpDir);
            foreach my $file (@files)
            {
                ProcessFASTAFile($fpOut, "$faaDir/$file");
            }
            close($fpOut);
        }
        else
        {
            printf STDERR "Error: Cannot write combined .faa file ($outFile)\n";
            exit 1;
        }
        closedir($fpDir);
    }
    else
    {
        printf STDERR "Error: Cannot open directory of .faa files ($faaDir)\n";
        exit 1;
    }
}

sub ProcessFASTAFile
{
    my($fpOut, $file) = @_;

    if(open(my $fp, '<', $file))
    {
        my $header = '';
        my $seq1   = '';
        my $seq2   = '';
        my $line   = 0;
        while(<$fp>)
        {
            chomp;
            if($line == 0)
            {
                $header = $_;
            }
            elsif($line == 1)
            {
                $seq1 = $_;
            }
            elsif($line == 2)
            {
                $seq2 = $_;
            }
            $line++;
        }
        close $fp;
        my $id = $header;
        $id =~ s/\|.*$//;
        if($seq1 ne '')
        {
            print $fpOut "${id}_1\n";
            print $fpOut "$seq1\n";
        }
        if($seq2 ne '')
        {
            print $fpOut "${id}_2\n";
            print $fpOut "$seq2\n";
        }
    }
}














































sub ReadClusters
{
    my ($clsFile) = @_;
    my %clusters = ();
    my @members  = ();
    my $index    = 0;

    if(open(my $fp, '<', $clsFile))
    {
        while(<$fp>)
        {
            chomp;
            if(/^>/)
            {
                if(scalar(@members))
                {
                    foreach my $key (@members)
                    {
                        $index = 0;
                        foreach my $member (@members)
                        {
                            if($member ne $key)
                            {
                                $clusters{$key}[$index++] = $member;
                            }
                        }
                    }
                }
                @members = ();
            }
            else
            {
                my @fields = split;
                my $member = $fields[2];
                $member =~ s/^\>//;
                $member =~ s/\|.*//;
#                print "Storing $member\n";
                push @members, $member;
            }
        }
        close $fp;

        # and the last one
        if(scalar(@members))
        {
            foreach my $key (@members)
            {
                $index = 0;
                foreach my $member (@members)
                {
                    if($member ne $key)
                    {
                        $clusters{$key}[$index++] = $member;
                    }
                }
            }
        }
    }
    return(%clusters);
}

sub MakeFAA
{
    my($fastaDir, $tmpDir, $chain) = @_;

    my $outFile = "$tmpDir/all_${chain}.faa";

    if(open(my $fh, '>', $outFile))
    {
        if(opendir(my $fd, $fastaDir))
        {
            my @files = grep /\.faa/, readdir($fd);
            closedir($fd);
            foreach my $inFile (@files)
            {
                my $content = GetLorH("$fastaDir/$inFile", $chain);
                print $fh $content;
            }
        }
        else
        {
            return('');
        }
        close($fh);
    }
    else
    {
        return('');
    }
    return($outFile);
}

sub GetLorH
{
    my($inFile, $chain) = @_;
    my $content = '';
    if(open(my $fp, '<', $inFile))
    {
        my $header = <$fp>;
        my $light  = <$fp>;
        my $heavy  = <$fp>;
        close $fp;

        $header =~ s/L_H\s*$/$chain/;
        $content = "$header\n";
        my $sequence = ($chain eq 'L')?$light:$heavy;
        $sequence =~ s/\s//g;
        return('') if($sequence eq '');
        $content .= "$sequence\n"; 
    }
    return($content);
}

sub Die
{
    my($msg) = @_;
    print STDERR "$msg\n";
    exit 1;
}

sub MakeTempDir
{
    my $tmpDir = "/var/tmp/abdb" . $$ . time();
    `mkdir -p $tmpDir`;
    return '' if(! -d $tmpDir);
    return ($tmpDir);
}


