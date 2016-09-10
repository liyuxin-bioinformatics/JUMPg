#!/usr/bin/perl -w 

use strict;

if (scalar(@ARGV)!=2)
{
        die "Usage: perl pepToReads.pl file.readCounts out \n";
}

my $infile=$ARGV[0];
open(RST,">$ARGV[1]");

###############################################
my $maxsize=5000000;
###############################################
system("cp $infile in.tmp");
my (%uniqf);
my $intmpSize= -s "in.tmp";

while ($intmpSize>0)
{
	undef %uniqf; #my %firstID;

	open(IN,"in.tmp");
	open(OUT,">out.tmp");

	while(<IN>)
	{
		chomp;my $line=$_;
		my @tmp=split(/\t/,$_);
		my @aa=split(//,$tmp[0]); pop(@aa); shift(@aa);
		my $seq=join('',@aa);
		
		if (exists($uniqf{$seq})) 
		{
			$uniqf{$seq}+=$tmp[2];
		}
		else
		{
			if (scalar(keys %uniqf)<=$maxsize)
			{
				$uniqf{$seq}=$tmp[2];
				#$firstID{$seq}=$id;
			}
			else
			{
				print OUT "$line\n";
			}
		}
		
	}

	close IN;
	close OUT;

	foreach my $seq (keys  %uniqf)
	{
		print RST "$seq\t$uniqf{$seq}\n";
	}

	system("mv out.tmp in.tmp");
	$intmpSize= -s "in.tmp";
}

close RST;
