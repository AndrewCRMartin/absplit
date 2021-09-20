#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) Prof. Andrew C. R. Martin, UCL, 2020
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
# Add the path of the executable to the library path
use FindBin;
use lib $FindBin::Bin;
use strict;
use fasta;

@::keyResidues = qw/L38 L40 L41 L44 L46 L87 H33 H42 H45 H60 H62 H91 H105/;

use Cwd qw(abs_path);
my $dataDir      = abs_path("$FindBin::Bin/../../data");
my $cdhitFile    = "$dataDir/cdhit.faa";
my $templateFile = "$dataDir/templates.faa";

if((! -e $templateFile) || defined($::f) || defined($::force))
{
    if(open(my $in, '<', $cdhitFile))
    {
        if(open(my $out, '>', $templateFile))
        {
            my($id, $info, $sequence);
            print STDERR "Finding interface residues";
            while((($id, $info, $sequence) = fasta::ReadFasta($in)) && ($id ne ""))
            {
                print STDERR '.';
                my $header = FindIFResidues($info, $sequence);
                print $out "$header\n";
                print $out "$sequence\n";
            }
            print STDERR "done\n";
            close($out);
        }
        else
        {
            print STDERR "Error: unable to write template file: $templateFile\n";
            exit 1;
        }
        close($in);
    }
    else
    {
        print STDERR "Error: Unable to read CD-HIT file: $cdhitFile\n";
        exit 1;
    }
}

sub FindIFResidues
{
    my($header, $sequence) = @_;

    my $ifResidues   = '';
    my $tmpFastaFile = "/tmp/fifr_faa_$$" . time();
    my $tmpSeqFile   = "/tmp/fifr_seq_$$" . time();
    if(open(my $fp, '>', $tmpFastaFile))
    {
        print $fp "$header\n";
        print $fp "$sequence\n";
        close($fp);

        my $exe = "abnum -f -c $tmpFastaFile > $tmpSeqFile";
        `$exe`;

        $ifResidues = doFindIFResidues($tmpFastaFile, $tmpSeqFile);
        
        unlink $tmpFastaFile;
        unlink $tmpSeqFile;
    }
    return("$header|[$ifResidues]");
}

sub doFindIFResidues
{
    my($faaFile, $seqFile) = @_;

    my $pos        = 0;
    my $ifResidues = '';
    if(open(my $fp, '<', $seqFile))
    {
        while(<$fp>)
        {
            chomp;
            my($resID, $aa) = split;
            my @matches = grep(/^$resID$/, @::keyResidues);
            if($matches[0] ne '')
            {
#                print "> $resID\n";
                $ifResidues .= "$pos,";
            }
            $pos++;
        }
        chop $ifResidues;
        close $fp;
    }
    return($ifResidues);
}
