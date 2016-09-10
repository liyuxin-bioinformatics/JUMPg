#!/usr/bin/perl
use strict;
use warnings;
my %hash;
open(IN,"$ARGV[0]") or die "cannot open $ARGV[0]\n"; # multistage_RDT_pep_db125_reads.alg
my $line=<IN>;print $line;
while(<IN>)
{
	my @t=split /\t/,$_;
	my $pep=$t[3];

	if (!defined($hash{$pep})) { $hash{$pep}{psm}=0; $hash{$pep}{score}=0; }
	if (!defined($hash{$pep}) or $hash{$pep}{score}<$t[1]) {
		$hash{$pep}{score}=$t[1];
		$hash{$pep}{line}=$_;
	}
	$hash{$pep}{psm}++;
}
close IN;

foreach my $h (values %hash) 
{
	my @t=split /\t/, $$h{line};
	$t[4]=$$h{psm};
	print join("\t",@t);
}
