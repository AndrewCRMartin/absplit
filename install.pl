#!/usr/bin/perl -s

use lib '.';
use strict;
use config;

my $force = (defined($::f)||defined($::force))?'-force':'';

my %config = config::ReadConfig('absplit.conf');

my $bindir  = $config{'bindir'};
my $datasub = $config{'datasub'};
my $datadir = "$bindir/$datasub";              # Data for absplit
$config{'datadir'} = $datadir;                 # Put it in the config

# Build the template data if not present
if((! -e "./data/templates.faa") || ($force ne ''))
{
    `(cd dataprep; ./maketemplates/getpdbabseqs.pl $force)`;
    `(cd dataprep; ./maketemplates/findinteractingresidues.pl $force)`;
}
# Install the template data
`mkdir -p $datadir`;
`cp -R data/* $datadir`;

# Build and install the executable
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
