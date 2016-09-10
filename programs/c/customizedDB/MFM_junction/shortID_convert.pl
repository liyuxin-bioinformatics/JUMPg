open(IN,"$ARGV[0]");
open(IDS,">$ARGV[1].ids");
open(FAS,">$ARGV[1].fas");
$rd=int(rand(99999));
$k=0;
while (<IN>)
{
	$k++;
	chomp;
	s/>//;
	print IDS "$rd\_$k\t$_\n";
	print FAS ">$rd\_$k\n";

	$line=<IN>;
	print FAS $line;
}
close IN;
close IDS;
close FAS;
