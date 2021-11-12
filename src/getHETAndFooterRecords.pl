#!/usr/bin/perl

use strict;
my $output = '';

while(<>)
{
    if(/^TER/)
    {
        $output = '';
    }
    else
    {
        $output .= $_;
    }
}

print $output;
