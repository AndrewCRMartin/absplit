#!/usr/bin/perl

use strict;

$::gCDHit = '../dataprep/cdhit/cd-hit';

my $fastaDir = shift @ARGV;
my $tmpDir   = MakeTempDir();
Die("Cannot create temporary directory $tmpDir") if($tmpDir eq '');

my $lightFAA = MakeFAA($fastaDir, $tmpDir, 'L');
my $heavyFAA = MakeFAA($fastaDir, $tmpDir, 'H');
if(($lightFAA eq '') || ($heavyFAA eq ''))
{
    Die("Unable to read FASTA files from $fastaDir");
}

print "Temp dir: $tmpDir\n";

my %lightClusters = Cluster($lightFAA, $tmpDir, 'L');
my %heavyClusters = Cluster($heavyFAA, $tmpDir, 'H');
# These are now indexed by each ID and contain an array of other members
# of the cluster

my %clusters = MergeClusters(\%lightClusters, \%heavyClusters);

PrintClusters(%clusters);
    

sub MergeClusters
{
    my ($hLightClusters, $hHeavyClusters) = @_;

    my %clusters = ();
    my @deleted  = ();

    # Create clusters where they span heavy and light
    foreach my $key (sort keys %$hLightClusters)
    {
        if(defined($$hHeavyClusters{$key}))
        {
            foreach my $member (sort @{$$hLightClusters{$key}})
            {
                if(inArray($member, @{$$hHeavyClusters{$key}}))
                {
                    push(@{$clusters{$key}}, $member);
                    delete ($$hHeavyClusters{$member});
                    delete ($$hLightClusters{$member});
                    push(@deleted, $member);
                }
            }
            delete ($$hHeavyClusters{$key});
            delete ($$hLightClusters{$key});
            push(@deleted, $key);
        }
    }

    AddSingleChainClusters($hLightClusters, \%clusters, \@deleted);
    AddSingleChainClusters($hHeavyClusters, \%clusters, \@deleted);

    AddSingletons($hLightClusters, \%clusters, \@deleted);
    AddSingletons($hHeavyClusters, \%clusters, \@deleted);

    return(%clusters);
}

sub AddSingletons
{
    my($hChainClusters, $hClusters, $aDeleted) = @_;
    foreach my $key (sort keys %$hChainClusters)
    {
        if(!defined($$hClusters{$key}))
        {
            $$hClusters{$key} = [];
        }
    }
    
}

sub AddSingleChainClusters
{
    my($hChainClusters, $hClusters, $aDeleted) = @_;

    foreach my $key (sort keys %$hChainClusters)
    {
        my $hasMembers = 0;
        foreach my $member (sort @{$$hChainClusters{$key}})
        {
            if(!inArray($member, @$aDeleted))
            {
                push(@{$$hClusters{$key}}, $member);
                delete $$hClusters{$member};
                $hasMembers = 1;
            }
            delete ($$hChainClusters{$member});
            push(@$aDeleted, $member);
        }
        if($hasMembers)
        {
            delete ($$hChainClusters{$key});
            push(@$aDeleted, $key);
        }
    }
}

sub PrintClusters
{
    my(%clusters) = @_;
    foreach my $key (keys %clusters)
    {
        my $i = 0;
        print "$key";
        foreach my $member (sort @{$clusters{$key}})
        {
            print " $member";
            $i++;
        }
        print "\n";
    }
}

sub inArray
{
    my($value, @array) = @_;
    foreach my $item (@array)
    {
        return(1) if($item eq $value);
    }
    return(0);
}

sub Cluster
{
    my($faaFile, $tmpDir, $chain) = @_;
    `$::gCDHit -c 1.0 -i $faaFile -o $tmpDir/$chain.faa`;
    my $clsFile = "$tmpDir/$chain.faa.clstr";
    my %clusters = ReadClusters($clsFile);
}


sub ReadClusters
{
    my ($clsFile) = @_;
    my %clusters = ();
    my @members  = ();
    my $index    = 0;

    if(open(my $fp, '<', $clsFile))
    {
        while(<$fp>)
        {
            chomp;
            if(/^>/)
            {
                if(scalar(@members))
                {
                    foreach my $key (@members)
                    {
                        $index = 0;
                        foreach my $member (@members)
                        {
                            if($member ne $key)
                            {
                                $clusters{$key}[$index++] = $member;
                            }
                        }
                    }
                }
                @members = ();
            }
            else
            {
                my @fields = split;
                my $member = $fields[2];
                $member =~ s/^\>//;
                $member =~ s/\|.*//;
#                print "Storing $member\n";
                push @members, $member;
            }
        }
        close $fp;

        # and the last one
        if(scalar(@members))
        {
            foreach my $key (@members)
            {
                $index = 0;
                foreach my $member (@members)
                {
                    if($member ne $key)
                    {
                        $clusters{$key}[$index++] = $member;
                    }
                }
            }
        }
    }
    return(%clusters);
}

sub MakeFAA
{
    my($fastaDir, $tmpDir, $chain) = @_;

    my $outFile = "$tmpDir/all_${chain}.faa";

    if(open(my $fh, '>', $outFile))
    {
        if(opendir(my $fd, $fastaDir))
        {
            my @files = grep /\.faa/, readdir($fd);
            closedir($fd);
            foreach my $inFile (@files)
            {
                my $content = GetLorH("$fastaDir/$inFile", $chain);
                print $fh $content;
            }
        }
        else
        {
            return('');
        }
        close($fh);
    }
    else
    {
        return('');
    }
    return($outFile);
}

sub GetLorH
{
    my($inFile, $chain) = @_;
    my $content = '';
    if(open(my $fp, '<', $inFile))
    {
        my $header = <$fp>;
        my $light  = <$fp>;
        my $heavy  = <$fp>;
        close $fp;

        $header =~ s/L_H\s*$/$chain/;
        $content = "$header\n";
        my $sequence = ($chain eq 'L')?$light:$heavy;
        $sequence =~ s/\s//g;
        return('') if($sequence eq '');
        $content .= "$sequence\n"; 
    }
    return($content);
}

sub Die
{
    my($msg) = @_;
    print STDERR "$msg\n";
    exit 1;
}

sub MakeTempDir
{
    my $tmpDir = "/var/tmp/abdb" . $$ . time();
    `mkdir -p $tmpDir`;
    return '' if(! -d $tmpDir);
    return ($tmpDir);
}


