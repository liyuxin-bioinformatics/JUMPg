#!/usr/bin/perl -w 

use strict;
use Storable;

if (scalar(@ARGV)!=2)
{
        die "Usage: perl  fragment.fas out_prefix\n";
}

my $infile=$ARGV[0];
open(FAS,">$ARGV[1].fas");
open(RDC,">$ARGV[1].readCounts");
open(ID,">$ARGV[1].seqIDs");

###############################################
my $maxsize=100000;
###############################################
system("cp $infile in.tmp");
my (%uniqf);
my $intmpSize= -s "in.tmp";

while ($intmpSize>0)
{
	undef %uniqf; my %firstID;

	open(IN,"in.tmp");
	open(OUT,">out.tmp");

	while(<IN>)
	{
		chomp;
		my $id=$_; $id =~ s/^>//;
		
		my $seq=<IN>;
		chomp($seq);

		if (exists($uniqf{$seq})) 
		{
			$uniqf{$seq}{$id}='';
		}
		else
		{
			if (scalar(keys %uniqf)<=$maxsize)
			{
				$uniqf{$seq}{$id}='';
				$firstID{$seq}=$id;
			}
			else
			{
				print OUT '>',"$id\n$seq\n";
			}
		}
		
	}

	close IN;
	close OUT;

	foreach my $seq (keys  %uniqf)
	{
		print FAS '>',$firstID{$seq},"\n$seq\n";
		print RDC "$seq\t$firstID{$seq}\t",scalar(keys %{$uniqf{$seq}}),"\n";
		print ID "$firstID{$seq}\t$seq\t";
		foreach my $seqID (keys %{$uniqf{$seq}}) { print ID "$seqID,"; }  print ID "\n";
	}


	system("mv out.tmp in.tmp");
	$intmpSize= -s "in.tmp";
	print "in.tmp file size = $intmpSize\n";
}

close FAS;
close RDC;
close ID;
