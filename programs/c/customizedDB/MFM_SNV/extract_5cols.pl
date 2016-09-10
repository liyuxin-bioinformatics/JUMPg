open(IN,"$ARGV[0]");
while(<IN>)
{
	next if /^NormalSample/;
	chomp;
	@t=split(/\t/,$_);
	if ($t[5] eq 'SNP') {print "$t[3]\t$t[4]\t$t[4]\t$t[9]\t$t[10]\n";}
	elsif ($t[5] eq 'insertion') {print "$t[3]\t$t[4]\t$t[4]\t-\t$t[10]\n";}
	elsif ($t[5] eq 'deletion') {print "$t[3]\t$t[4]\t",$t[4]+length($t[9])-1,"\t$t[9]\t-\n";}

}
close IN;
