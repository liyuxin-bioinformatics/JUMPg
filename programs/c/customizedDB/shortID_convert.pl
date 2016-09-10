open(IN,"$ARGV[0]");
open(IDS,">$ARGV[1].ids");
open(FAS,">$ARGV[1].fas");
$rd=int(rand(99999));
$k=0;
$starter='cu|';
$ender='|cu';
while (<IN>)
{
	$k++;
	chomp;
	s/>//;
	#print IDS "$rd\_$k\t$_\n";
	#print FAS ">$rd\_$k\n";
	$id="$starter$rd\_$k$ender";
	print IDS "$id\t$_\n";
	#print FAS ">$id\n";
	print FAS ">$id $_\n";

	$line=<IN>;
	print FAS $line;
}
close IN;
close IDS;
close FAS;
