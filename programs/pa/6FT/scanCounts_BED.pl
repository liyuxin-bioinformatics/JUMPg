#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)==0)
{
	die "perl scanCounts_BED.pl scan_based.BED\n";
}

my %lines;
# %lines{$line without outfile}{$outfile}
open(IN,$ARGV[0]) || die "Cannot open $ARGV[0]\n";
my $line=<IN>; print $line;
while(<IN>)
{
	chomp;
	my @t=split /\t/, $_;
	my $outfile=pop @t;
	$line=join("\t",@t);
	$lines{$line}{$outfile}='';
	#if (!defined($lines{$line})) { $lines{$_}=0; }
	#$lines{$_}++;
}
close IN;

foreach my $line (keys %lines)
{
	my @t=split /\t/, $line;
	#my $count=$lines{$line};
	my $count=scalar(keys %{$lines{$line}});
	$t[3]="$t[3]\[$count\]";
=head
	if ($count==1) { $t[4]=200; }
	elsif ($count==2) { $t[4]=300; }
	elsif (3<=$count and $count<=5) { $t[4]=400; }
	elsif (6<=$count and $count<=10) { $t[4]=550; }
	elsif (11<=$count and $count<=30) { $t[4]=650; }
	elsif (31<=$count and $count<=100) { $t[4]=750; }
	elsif (101<=$count) { $t[4]=900; }
=cut
	for (my $i=0; $i<=$#t-1; $i++) { print "$t[$i]\t"; }
	print "$t[$#t]\n";
}
