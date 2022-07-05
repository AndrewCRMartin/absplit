#!/usr/bin/perl -s

use strict;

$::gCDHit = '../dataprep/cdhit/cd-hit';

my $fastaDir = shift @ARGV;
my $tmpDir   = MakeTempDir();
Die("Cannot create temporary directory $tmpDir") if($tmpDir eq '');
my $faaFile = "$tmpDir/all.faa";

CombineAbsplitFastaFiles($fastaDir, $faaFile);
print "Temp dir: $tmpDir\n" if(defined($::d));
my $sequenceClusterFile = ClusterSequences($faaFile, $tmpDir);
my @abClusters = CreateFinalAntibodyClusters($sequenceClusterFile);
print "# Free Antibody:Complexed Antibody\n";
foreach my $abCluster (@abClusters)
{
    my $reordered = Reorder($abCluster);
    print "$reordered\n";
}

unlink($tmpDir) if(!defined($::d));

#-----------------------------------------------------------------------
sub Reorder
{
    my($cluster) = @_;
    my $debug = '';

    # If we have the debugging line numbers, split them off
    if($cluster =~ /(.*)\s*:\s*(.*)/)
    {
        $debug   = $1;
        $cluster = $2;
    }

    # Split into individual fields
    my @items = split(/\s+/, $cluster);

    # Reassemble into separate free and complex lists
    my @free = ();
    my @complexed = ();
    for my $item (@items)
    {
        if($item =~ /\d$/) # ID ends with a number, so not a complex
        {
            push @free, $item;
        }
        else
        {
            push @complexed, $item;
        }
    }

    # Assemble the final description
    my $result = '';
    if($debug ne '')
    {
        $result = "$debug:";
    }
    $result .= join(',',@free);
    $result .= ':';
    $result .= join(',',@complexed);
    
    return($result);
}


#-----------------------------------------------------------------------
sub ClusterSequences
{
    my($faaFile, $tmpDir) = @_;
    `$::gCDHit -c 1.0 -i $faaFile -o $tmpDir/clusters.faa`;
    my $CDHitClsFile = "$tmpDir/clusters.faa.clstr";
    my @clusters = ReadCDHitClusters($CDHitClsFile);

    my $clsFile = "$tmpDir/seqClusters.dat";
    if(open(my $clsFp, '>', $clsFile))
    {
        foreach my $cluster (@clusters)
        {
            print $clsFp "$cluster\n";
        }
        close $clsFp;
    }
    else
    {
        Die("Cannot write sequence cluster file $clsFile");
    }

    return($clsFile);
}


#-----------------------------------------------------------------------
sub CombineAbsplitFastaFiles
{
    my($faaDir, $outFile) = @_;
    
    if(opendir(my $fpDir, $faaDir))
    {
        if(open(my $fpOut, '>', $outFile))
        {
            my @files = grep /\.faa/, readdir($fpDir);
            foreach my $file (@files)
            {
                ProcessAbsplitFastaFile($fpOut, "$faaDir/$file");
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


#-----------------------------------------------------------------------
sub ProcessAbsplitFastaFile
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


#-----------------------------------------------------------------------
sub ReadCDHitClusters
{
    my ($clsFile) = @_;
    my @clusters = ();
    my $index    = 0;
    my @members  = ();

    if(open(my $fp, '<', $clsFile))
    {
        while(<$fp>)
        {
            chomp;
            if(/^>/)
            {
                if(scalar(@members))
                {
                    my $output = join(' ', @members);
                    push @clusters, $output;
                }
                @members = ();
            }
            else
            {
                my @fields = split;
                my $member = $fields[2];
                $member =~ s/^\>//;
                $member =~ s/\.\.\..*//;
                print STDERR "Storing $member\n" if(defined($::d));
                push @members, $member;
            }
        }
        close $fp;

        # and the last one
        if(scalar(@members))
        {
            my $output = join(' ', @members);
            push @clusters, $output;
        }
    }
    return(@clusters);
}


#-----------------------------------------------------------------------
sub Die
{
    my($msg) = @_;
    print STDERR "$msg\n";
    exit 1;
}


#-----------------------------------------------------------------------
sub MakeTempDir
{
    my $tmpDir = "/var/tmp/abdb" . $$ . time();
    `mkdir -p $tmpDir`;
    return '' if(! -d $tmpDir);
    return ($tmpDir);
}


#-----------------------------------------------------------------------
sub CreateFinalAntibodyClusters
{
    my ($seqClusFile) = @_;
    
    my @data       = ();
    my @clusters   = ();
    my @scClusters = ();
    my @results    = ();

    if(open(my $fp, '<', $seqClusFile))
    {
        while(<$fp>)
        {
            chomp;
            my(@ids) = split;
            push @data, \@ids;
        }
        close($fp);
    }
    else
    {
        Die("Cannot read sequence cluster file $seqClusFile\n");
    }

    my $row = 0;
    foreach my $aRow (@data)
    {
        while(scalar(@{$aRow}))
        {
            # Take the last item from this row
            my $item = pop(@$aRow);
            # Construct the name of its partner
            my ($partner, $stem) = GetPartnerName($item);
            # Find which row the partner is in
            my $partnerRow = FindPartnerRow($partner, @data);
            print("$item $partner $partnerRow\n") if(defined($::d));
            if($partnerRow >= 0)
            {
                push(@{$clusters[$row][$partnerRow]}, $stem);
                RemovePartner($partner, \@{$data[$partnerRow]});
            }
            else
            {
                push(@{$scClusters[$row]}, $stem);
            }
        }
        
        $row++;
    }
    
    my $nRows = $row;
    
    for(my $i=0; $i<$nRows; $i++)
    {
        for(my $j=0; $j<$nRows; $j++)
        {
            if(defined($clusters[$i][$j]))
            {
                my $data = '';
                $data .= "$i $j: " if(defined($::d));
                foreach my $item (@{$clusters[$i][$j]})
                {
                    $data .= "$item ";
                }
                push @results, $data;
            }
        }
    }
    
    my $i=0;
    foreach my $aRow (@scClusters)
    {
        if(defined($aRow))
        {
            my $data = '';
            $data .= "$i: " if(defined($::d));
            foreach my $item (@$aRow)
            {
                $data .= "$item ";
            }
            push @results, $data;
        }
        $i++;
    }

    return(@results);
}


#-----------------------------------------------------------------------
sub RemovePartner
{
    my($partner, $aRow) = @_;
    @$aRow = grep ! /$partner/, @$aRow;
}


#-----------------------------------------------------------------------
sub FindPartnerRow
{
    my($partner, @data) = @_;
    my $row = 0;
    foreach my $aRow (@data)
    {
        foreach my $item (@$aRow)
        {
            if($item eq $partner)
            {
                return($row);
            }
        }
        $row++;
    }
    return(-1);
}


#-----------------------------------------------------------------------
sub GetPartnerName
{
    my($item) = @_;
    if($item =~ '_1$')
    {
        $item =~ s/_1$/_2/;
    }
    else
    {
        $item =~ s/_2$/_1/;
    }
    my $stem = $item;
    $stem =~ s/_[12]$//;
    return($item, $stem);
}

