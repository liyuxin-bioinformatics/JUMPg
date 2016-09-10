open(IN,$ARGV[0]) || die "cannot open $ARGV[0]\n";
while(<IN>)
{
	chomp;
	#>BCCIP|uc001ljb.4|p.E182D|chr10:127520123-127520123.A_T
	s/>//; 
	my $header=$_;
	my @t=split /\|/,$_;
	my $seq=<IN>;
	chomp($seq);
	$hash{$t[$#t]}{$seq}=$header;
}
close IN;

open(OUT,">$ARGV[1]") || die "cannot open $ARGV[1]\n";
foreach my $h (values %hash)
{
	foreach my $seq (keys %{$h})
	{
		print OUT ">$$h{$seq}\n$seq\n";
	}
}
close OUT;
