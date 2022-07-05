#!/usr/bin/perl

use strict;

my $dir = shift(@ARGV);

if(opendir(my $fpDir, $dir))
{
    my @files = grep /\.faa/, readdir($fpDir);
    foreach my $file (@files)
    {
        ProcessFile($file);
    }
    closedir($fpDir);
}
else
{
    printf STDERR "Error: Cannot open directory of .faa files ($dir)\n";
    exit 1;
}

sub ProcessFile
{
    my($file) = @_;
    #    print "$file\n";
    if(open(my $fp, '<', $file))
    {
        my $header = '';
        my $seq1   = '';
        my $seq2   = '';
        my $line   = 0;
        while(<>)
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
            print ">${id}_1\n";
            print "$seq1\n";
        }
        if($seq2 ne '')
        {
            print ">${id}_2\n";
            print "$seq2\n";
        }
    }
}
