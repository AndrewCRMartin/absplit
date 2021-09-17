#!/usr/bin/perl

use lib '.';
use strict;
use config;

my %config = config::ReadConfig('absplit.conf');

my $bindir  = $config{'bindir'};
my $datasub = $config{'datasub'};
my $datadir = "$bindir/$datasub";              # Data for absplit
$config{'datadir'} = $datadir;                 # Put it in the config

if(! -e "$datadir/cdhit.faa")
{
    `(cd dataprep; ./maketemplates/getpdbabseqs.pl)`;
}

MakeMake(%config);
`(cd src; make; make install)`;

sub MakeMake
{
    my(%config) = @_;

    my $infile  = "./src/Makefile.src";
    my $outfile = "./src/Makefile";

    if(open(my $in, '<', $infile))
    {
        if(open(my $out, '>', $outfile))
        {
            while(<$in>)
            {
                Process(\$_, %config);
                print $out $_;
            }
            close $out;
        }
        else
        {
            die "Can't write $outfile";
        }
        close $in;
    }
    else
    {
        die "Can't write $infile";
    }
}

sub Process
{
    my ($pLine, %config) = @_;

    foreach my $key (keys %config)
    {
        my $value = $config{$key};
        $$pLine =~ s/\{$key\}/$value/g;
    }
}
