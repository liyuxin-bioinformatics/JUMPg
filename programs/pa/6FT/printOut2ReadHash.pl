#!/usr/bin/perl -w
use strict;
use Storable;

if (scalar(@ARGV)<2)
{
	die "usage: perl  printOut2ReadHash.pl out2read.hash output\n";
}


#load hash
my $out2read=retrieve($ARGV[0]);

#output
open(OUT,">$ARGV[1]");
print OUT "outfile\tXCorr\tdCn\tpeptide\tscanCount\treadCount\treadIDs\n";
foreach my $out (keys %$out2read)
{
	print OUT "$out\t$$out2read{$out}{XCorr}\t$$out2read{$out}{dCn}\t$$out2read{$out}{peptide}\t$$out2read{$out}{scanCount}\t",scalar(keys %{$$out2read{$out}{readID}}),"\t";
	foreach my $readID (keys %{$$out2read{$out}{readID}}) { print OUT "$readID,"; }
	print OUT "\n";
}
close OUT;
