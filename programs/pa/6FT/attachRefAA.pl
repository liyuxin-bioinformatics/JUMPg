#!/usr/bin/perl -I /home/yli4/lib/perl5

use PrimarySeq;

$tr=PrimarySeq->new();

open(IN,"$ARGV[0]"); # MM_RDT_baseQC_db137_test1_reads.out1.alg.updated2.pep.gnSb
$line=<IN>; chomp($line);
print "$line\trefAA\tAAchangeCount\tAAchangeRepresent\tAAchanges\n";
while(<IN>)
{
        chomp;
        @t=split /\t/,$_;
	$strand=$t[$#t-6];
	$nt=$t[$#t];

	if ($strand eq '-') { $nt=$tr->revcom($nt); }
	$refAA=$tr->translate($nt);

	$altAA=$t[3];
	$altAA =~ s/[\*\@\#]//g;
	
	$k=0; $AAchanges=$AAchangeRepresent='';
	if (length($altAA) eq length($refAA))
	{
		@alt=split //,$altAA;
		@ref=split //,$refAA;

		for ($i=0;$i<=$#ref;$i++)
		{
			if ($ref[$i] ne $alt[$i])
			{
				$k++;
				$AAchangeRepresent="$ref[$i]$i$alt[$i]";
				$AAchanges.=$AAchangeRepresent;
				$AAchanges.=',';

				$AAchangeRepresent="$ref[$i]\_$alt[$i]";
			}
		}
	}

	print "$_\t$refAA\t$k\t$AAchangeRepresent\t$AAchanges\n";
}

