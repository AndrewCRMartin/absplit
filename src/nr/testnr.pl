#!/usr/bin/perl -s
use strict;
my @data = ();
my @clusters = ();
my @scClusters = ();
while(<>)
{
    chomp;
    my(@ids) = split;
    push @data, \@ids;
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
            print "$i $j: " if(defined($::d));
            foreach my $item (@{$clusters[$i][$j]})
            {
                print "$item ";
            }
            print "\n";
        }
    }
}

printf("Single chain clusters...\n") if(defined($::d));
my $i=0;
foreach my $aRow (@scClusters)
{
    if(defined($aRow))
    {
        print "$i: " if(defined($::d));
        foreach my $item (@$aRow)
        {
            print "$item ";
        }
        print "\n";
    }
    $i++;
}



sub RemovePartner
{
    my($partner, $aRow) = @_;
    @$aRow = grep ! /$partner/, @$aRow;
}

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

