#!/usr/bin/perl
use strict;
my @data = ();
my @clusters = ();
#my @scClusters = ();
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
#        printf("%d\n", scalar(@{$aRow}));
#        print "@{$aRow}\n";
        # Take the last item from this row
        my $item = pop(@$aRow);
        # Construct the name of its partner
        my ($partner, $stem) = GetPartnerName($item);
        # Find which row the partner is in
        my $partnerRow = FindPartnerRow($partner, @data);
        if($partnerRow >= 0)
        {
            push(@{$clusters[$row][$partnerRow]}, $stem);
            RemovePartner($partner, $data[$row]);
        }
#        else
#        {
#            push(@{$scClusters[$row]}, $stem);
#        }
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
            print "$i $j: ";
            foreach my $item (@{$clusters[$i][$j]})
            {
                print "$item ";
            }
            print "\n";
        }
    }
}

#printf("Single chain clusters...\n");
#foreach my $aRow (@scClusters)
#{
#    if(defined($aRow))
#    {
#        foreach my $item (@$aRow)
#        {
#            print "$item ";
#        }
#        print "\n";
#    }
#}



sub RemovePartner
{
    my($partner, $pRow) = @_;
    @{$pRow} = grep ! /$partner/, @{$pRow};
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
    if($item =~ '_1')
    {
        $item =~ s/_1/_2/;
    }
    else
    {
        $item =~ s/_2/_1/;
    }
    my $stem = $item;
    $stem =~ s/_[12]//;
    return($item, $stem);
}

