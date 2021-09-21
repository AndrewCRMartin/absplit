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
%::cdrDefs = ('L1' => ['L24', 'L34'],
              'L2' => ['L50', 'L56'],
              'L3' => ['L89', 'L97'],
              'H1' => ['H26', 'H35'],
              'H2' => ['H50', 'H65'],
              'H3' => ['H95', 'H102']);

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
            print STDERR "Finding interface and CDR residues";
            while((($id, $info, $sequence) = fasta::ReadFasta($in)) && ($id ne ""))
            {
                print STDERR '.';
                my $header = FindInteractingResidues($info, $sequence);
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
else
{
    print STDERR "Annotated sequence data already exists. Use -force to rewrite\n";
}

sub FindInteractingResidues
{
    my($header, $sequence) = @_;

    my $ifResidues   = '';
    my $cdrResidues  = '';
    my $tmpFastaFile = "/tmp/fifr_faa_$$" . time();
    my $tmpSeqFile   = "/tmp/fifr_seq_$$" . time();
    if(open(my $fp, '>', $tmpFastaFile))
    {
        print $fp "$header\n";
        print $fp "$sequence\n";
        close($fp);

        my $exe = "abnum -f -c $tmpFastaFile > $tmpSeqFile";
        `$exe`;

        $ifResidues  = FindIFResidues($tmpFastaFile, $tmpSeqFile);
        $cdrResidues = FindCDRResidues($tmpFastaFile, $tmpSeqFile);
        
        unlink $tmpFastaFile;
        unlink $tmpSeqFile;
    }
    return("$header|[$ifResidues]|[$cdrResidues]");
}

sub FindIFResidues
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

sub FindCDRResidues
{
    my($faaFile, $seqFile) = @_;

    my $pos         = 0;
    my $cdrResidues = '';
    if(open(my $fp, '<', $seqFile))
    {
        while(<$fp>)
        {
            chomp;
            my($resID, $aa) = split;
            foreach my $cdr (keys %::cdrDefs)
            {
                if(resGE($resID, $::cdrDefs{$cdr}[0]) &&
                   resLE($resID, $::cdrDefs{$cdr}[1]))
                {
#                    print "> $resID\n";
                    $cdrResidues .= "$pos,";
                }
            }
            $pos++;
        }
        chop $cdrResidues;
        close $fp;
    }
    return($cdrResidues);
}



#*************************************************************************
#> ($chain, $resnum, $insert) = ParseResID($resid)
#  -----------------------------------------------
#  Input:   string   $resid    Residue identifier (e.g. L27A)
#  Returns: string   $chain    The chain label (or undef)
#           int      $resnum   The residue number (or undef)
#           string   $insert   The insert code (or undef)
#
#  Parse a residue identifier to extract the chain, residue number and
#  insert code
#
#  19.09.13  Original  By: ACRM
sub ParseResID
{
    my ($resid) = @_;

    if(!($resid =~ /([a-zA-Z]?)(\d+)([a-zA-Z]?)/))
    {
        return(undef, undef, undef);
    }
    return($1, $2, $3);
}

#*************************************************************************
#> BOOL resLE($res1, $res2)
#  ------------------------
#  Input:   string   $res1    First residue ID
#           string   $res2    Second residue ID
#  Return:  BOOL              Is $res1 <= $res2
#
#  Tests whether a residue ID is <= another residue ID
#
#  19.09.13  Original  By: ACRM
sub resLE
{
    my($res1, $res2) = @_;

    my @resA = ParseResID($res1);
    my @resB = ParseResID($res2);

    if($resA[0] eq $resB[0])         # Chain matches
    {
        if($resA[1] < $resB[1])      # ResNum1 less than ResNum2
        {
            return(1);
        }
        elsif($resA[1] == $resB[1])  # ResNum1 == ResNum2
        {
            if($resA[2] le $resB[2]) # Insert1 <= Insert2
            {
                return(1);
            }
        }
    }

    return(0);
}

#*************************************************************************
#> BOOL resGE($res1, $res2)
#  ------------------------
#  Input:   string   $res1    First residue ID
#           string   $res2    Second residue ID
#  Return:  BOOL              Is $res1 >= $res2
#
#  Tests whether a residue ID is >= another residue ID
#
#  19.09.13  Original  By: ACRM
sub resGE
{
    my($res1, $res2) = @_;

    my @resA = ParseResID($res1);
    my @resB = ParseResID($res2);

    if($resA[0] eq $resB[0])
    {
        if($resA[1] > $resB[1])
        {
            return(1);
        }
        elsif($resA[1] == $resB[1])
        {
            if($resA[2] ge $resB[2])
            {
                return(1);
            }
        }
    }

    return(0);
}

