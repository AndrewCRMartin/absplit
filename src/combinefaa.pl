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
    my $parts = '';
    $parts = AddPart($parts, 'L', 'ChainL', %seqs);
    $parts = AddPart($parts, 'l', 'Chainl', %seqs);
    $parts = AddPart($parts, 'H', 'ChainH', %seqs);
    $parts = AddPart($parts, 'h', 'Chainh', %seqs);
    
    print ">$label|$parts\n";
    print "\U$seqs{'ChainL'}\n" if(defined($seqs{'ChainL'}));
    print "\U$seqs{'Chainl'}\n" if(defined($seqs{'Chainl'}));
    print "\U$seqs{'ChainH'}\n" if(defined($seqs{'ChainH'}));
    print "\U$seqs{'Chainh'}\n" if(defined($seqs{'Chainh'}));
}

sub AddPart
{
    my($parts, $label, $key, %seqs) = @_;

    if(defined($seqs{$key}))
    {
        if($parts eq '')
        {
            $parts = $label;
        }
        else
        {
            $parts .= "_$label";
        }
    }
    return($parts);
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
