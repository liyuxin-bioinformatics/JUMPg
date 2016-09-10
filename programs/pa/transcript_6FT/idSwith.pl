=head
use IDtxt_parser;
$idp=IDtxt_parser->new;
$o2n=$idp->parseTab($ARGV[0],0); # .ids
=cut
open(IN,$ARGV[0]); # .ids
while(<IN>) {
	chomp;
	($o,$n)=split /\t/,$_;
	$o2n->{$o}=$n;
}
close IN;


open(IN,$ARGV[1]); # ID.txt
#$line=<IN>;
#$line=<IN>;
while(<IN>) {
	chomp;
	next if (/^Database/ or /^Peptide/);
	@t=split /\;/,$_;
	$t[1]=$o2n->{$t[1]};
	print join(";",@t),"\n";
}
close IN;
