#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)==0)
{
	die "perl scanCounts_BED.pl .BED\n";
}

my %lines;
open(IN,$ARGV[0]) || die "Cannot open $ARGV[0]\n";
my $line=<IN>; print $line;
while(<IN>)
{
	chomp;
	my ($chr,$start,$end,$name,$reads,$strand,$s1,$e1,$col,$blkN,$blkSize,$blkStart)=split(/\t/,$_);
	if ($blkN>1) {$col='0,128,0'; } # junc.
	else {$col='0,0,0'; } # body
	print "$chr\t$start\t$end\t$name\t$reads\t$strand\t$s1\t$e1\t$col\t$blkN\t$blkSize\t$blkStart\n"; 
}
close IN;

