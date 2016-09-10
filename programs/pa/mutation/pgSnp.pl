open(IN,$ARGV[0]) || die "Cannot open $ARGV[0]\n"; #mutPep5_mut_all_dbSNP.txt
#print 'track type=pgSnp visibility=3 db=hg19 name="pgSnp" description="SJMM_MS_spt_variants"',"\n";
print 'track type=pgSnp visibility=3 db=hg19 name="pgSnp" description="',$ARGV[1],'"',"\n";
$line=<IN>;
while(<IN>)
{
	chomp;
	($pos,$other)=split /\t/, $_;
	($pos,$alleles)=split /\./, $pos;
	($chr,$start,$end)=split /[\:\-]/, $pos;
	($ref,$mut)=split /\_/,$alleles;
        #my ($line,$type, $genes,$chr,$start,$end,$ref,$mut)=split /\t/, $_;
	#next if ($type =~ m/^synonymous/);
	print "$chr\t$start\t$end\t$ref\/$mut\t2\t0,0\t0,0\n";
}
close IN;
