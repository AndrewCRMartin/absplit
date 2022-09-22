#!/usr/bin/perl

my %lightClusters = ();
$lightClusters{'1r70_0'} = [
    '1qtj_2',
    '1r70_1',
    '1qtj_3',
    '2esg_0'];
$lightClusters{'1qtj_2'} = [
    '1r70_0',
    '1r70_1',
    '1qtj_3',
    '2esg_0'];
$lightClusters{'1r70_1'} = [
    '1qtj_2',
    '1r70_0',
    '1qtj_3',
    '2esg_0'];
$lightClusters{'1qtj_3'} = [
    '1r70_0',
    '1r70_1',
    '1qtj_2',
    '2esg_0'];
$lightClusters{'2esg_0'} = [
    '1qtj_2',
    '1r70_1',
    '1qtj_3',
    '1r70_0'];
$lightClusters{'1lll_0'} = [
    '1lll_2'];
$lightClusters{'1lll_2'} = [
    '1lll_0'];
$lightClusters{'0lll_0'} = [];

my %heavyClusters = ();
$heavyClusters{'1qtj_2'} = [
    '1r70_1',
    '1qtj_3',
    '2esg_0'];
$heavyClusters{'1r70_1'} = [
    '1qtj_2',
    '1qtj_3',
    '2esg_0'];
$heavyClusters{'1qtj_3'} = [
    '1r70_1',
    '1qtj_2',
    '2esg_0'];
$heavyClusters{'2esg_0'} = [
    '1qtj_2',
    '1r70_1',
    '1qtj_3'];
$heavyClusters{'1hhh_0'} = [
    '1hhh_2'];
$heavyClusters{'1hhh_2'} = [
    '1hhh_0'];
$heavyClusters{'0hhh_0'} = [];

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
