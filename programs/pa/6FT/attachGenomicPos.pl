open(IN,"$ARGV[0]"); # mouseBrain_RDT_peptides.sc.bed.check
$line=<IN>;
while(<IN>)
{
	chomp;
	@t=split /\t/,$_;
	$t[3] =~ s/\[\d+\]$//;;
	$hash{$t[3]}{strand}=$t[5];
	$hash{$t[3]}{chr}=$t[0];
	$hash{$t[3]}{start}=$t[1];
	$hash{$t[3]}{end}=$t[2];

	$hash{$t[3]}{exonCount}=$t[9];
	$hash{$t[3]}{exonSizes}=$t[10];
	$hash{$t[3]}{exonStarts}=$t[11];
}
close IN;

open(IN,"$ARGV[1]"); # mouseBrain_RDT_test1_reads.out3.alg.updated2.pep.check
$line=<IN>; chomp($line);
print "$line\tgenome_chr\tgenome_strand\tgenome_start\tgenome_end\texonCount\texonSizes\texonStarts\n";
while(<IN>)
{
	chomp;
	@t=split /\t/,$_;
	if (defined($hash{$t[3]}))
	{
		print "$_\t$hash{$t[3]}{chr}\t$hash{$t[3]}{strand}\t$hash{$t[3]}{start}\t$hash{$t[3]}{end}\t$hash{$t[3]}{exonCount}\t$hash{$t[3]}{exonSizes}\t$hash{$t[3]}{exonStarts}\n";
	}
	else { die "$t[3] not found in $ARGV[0]\n"; }
}
close IN;

