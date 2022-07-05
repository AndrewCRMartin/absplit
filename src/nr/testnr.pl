#!/usr/bin/perl
use strict;
my @data = ();
my $clusters = [];
while(<>)
{
    chomp;
    my(@ids) = split;
    push @data, \@ids;
}

my $row = 0;
foreach my $aRow (@data)
{
    # Take the last item from this row
    my $item = pop(@$aRow);
    # Construct the name of its partner
    my ($partner, $stem) = GetPartnerName($item);
    # Find which row the partner is in
    my $partnerRow = FindPartner($partner, @data);
    push(@($clusters[$row][$partnerRow]), $stem);
    RemovePartner($partner, $data[$row]);

    $row++;
}



foreach my $aRow (@data)
{
    foreach my $item (@$aRow)
    {
        print "$item ";
    }
    print "\n";
}


sub RemovePartner
{
    my($partner, $pRow) = @_;
    @{$pRow} = grep ! /$partner/, @{$pRow};
}

sub FindPartner
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

