#!/usr/bin/perl

use lib '.';
use strict;
use config;

my %config = config::ReadConfig('absplit.conf');

if($ARGV[0] ne '')
{
    $config{'bindir'} = $ARGV[0];
}
my $bindir  = $config{'bindir'};
my $datasub = $config{'datasub'};
my $datadir = "$bindir/$datasub";              # Data for absplit


if(! -e "$datadir/cdhit.faa")
{
    `(cd dataprep; ./maketemplates/getpdbabseqs.pl)`;
}
