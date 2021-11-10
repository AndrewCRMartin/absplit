#!/usr/bin/perl -s

use strict;

$::l="Undefined" if(!defined($::l));

# Read the sequence data
my %seqs  = ReadFAA();
# And write as a single sequence
WriteFAA($::l, %seqs);

sub WriteFAA
{
    my ($label, %seqs) = @_;
    print ">$label|L_H\n";
    print "\U$seqs{'ChainL'}\n";
    print "\U$seqs{'ChainH'}\n";
}

sub ReadFAA
{
    my $label = '';
    my %seqs  = ();
    while(<>)
    {
        chomp;
        
        if(/^\>(.*)/)
        {
            $label = $1;
        }
        else
        {
            if(defined($seqs{$label}))
            {
                $seqs{$label} .= $_;
            }
            else
            {
                $seqs{$label} = $_;
            }
        }
    }
    return(%seqs);
}
