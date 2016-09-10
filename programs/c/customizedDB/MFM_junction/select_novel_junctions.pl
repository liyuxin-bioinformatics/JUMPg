#!/usr/bin/perl

unless (scalar(@ARGV)==4)
{
	die "perl select_novel_junctions.pl Peng_SJMM_RNAseq_Exon_Junction.txt oneSample_readCut multiSample_readCut out\n";
}

#----------------------------------------------------------
my $single_ReadCut=$ARGV[1];
my $multi_ReadCut=$ARGV[2];
#----------------------------------------------------------

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[3]");
my ($line);

$line=<IN>;  print OUT $line;
while($line=<IN>)
{
	chomp($line);
	my @t=split(/\t/,$line);

	next if ($t[2] eq 'core');

	# check criteria
	my ($sampleN,$maxR)=(0,0);
	for (my $i=7;$i<=$#t;$i++)
	{
		if ($t[$i]) { $sampleN++; }
		if ( $t[$i]>$maxR ) { $maxR=$t[$i]; }
	}

	#print "$t[0]\t$sampleN\t$maxR\t$single_ReadCut\n";

	if ( ($sampleN>=2 and $maxR>=$multi_ReadCut) || $maxR>=$single_ReadCut )
	{
		print OUT "$line\n";
	}
}

close IN;
close OUT;
