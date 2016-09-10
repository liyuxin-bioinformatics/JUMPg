#!/usr/bin/perl -w -I /home/yli4/lib/perl5
use strict;
use Storable;

if (scalar(@ARGV)!=3)
{
	die "usage: perl extract_MS_reads.pl reads.fas out2read.hash output\n";
}

#load hash
my $out2read=retrieve($ARGV[1]);

#build %reads
my %reads;
foreach my $out (keys %{$out2read})
{
	foreach my $read (keys %{$$out2read{$out}}) { $reads{$read}=$out; }
}

#
open(IN,"$ARGV[0]");
open(OUT,">$ARGV[2]");

while(<IN>)
{
	chomp;
	my $id=$_;
	$id=substr($id,1,length($id)-1);

	my $seq=<IN>;
	chomp($seq);

	if (defined($reads{$id}))
	{
		print OUT ">$id\n$seq\n";
	}
}

close IN;
close OUT;
