#!/usr/bin/perl
use strict;
use warnings;

my $star_junc=shift @ARGV;
open(IN,$star_junc) || die "Cannot open STAR junction file $star_junc\n";
while(<IN>)
{
        chomp;
        my ($chr,$start,$end,$strand,$motifCode,$annotated,$uniReads,$multiReads,$maxOverhang)=split /\t/,$_;

	next unless ($chr =~ m/chr\d+/ || $chr =~ m/chr[XYM]/); # only consider major chrmosomes 
								# could expand to any chromosome in the future
	$chr =~ s/chr//;
	$start--;
	$end++;

	# 6:15533626:+,6:15585950:+
	print "$chr\:$start\:\+\,$chr\:$end\:\+\t$strand\n";
}

