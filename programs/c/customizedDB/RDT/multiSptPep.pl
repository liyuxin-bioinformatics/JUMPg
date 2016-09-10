use strict;
use warnings;

if ($#ARGV != 2)
{
	die "Usage: multiSptPep.pl file.readCounts cutoff out\n";
}

my $cutoff=$ARGV[1];

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[2]");
while(<IN>)
{
	chomp;
	my @t=split(/\t/,$_);
	if ($t[2]>=$cutoff)
	{
		print OUT "\>$t[1]\n$t[0]\n";
	}
}
close IN;
close OUT;
