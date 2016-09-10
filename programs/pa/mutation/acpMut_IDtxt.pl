#!/usr/bin/perl

use strict;
use warnings;

my %mutpep;
open(IN,$ARGV[0]) || die "cannot open $ARGV[0]\n"; # mutPep4_pep.txt
while(<IN>)
{
	next if (/^peptide/);
	my @t=split /\t/,$_;
	$mutpep{$t[0]}='';
}
close IN;

open(IN,$ARGV[1]) || die "cannot open $ARGV[1]\n"; # ID_mutation_level_3_origIDs.txt
while(<IN>)
{
	if (/^Database/) {print $_; next;}
	if (/^Peptide/) {print $_; next;}

	my @t=split /\;/,$_;
	my $pep=$t[0];
	#my ($pep,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos,$precursor_peak_intensity_percentage)=split(/\;/,$_);
	chop($pep);chop($pep);$pep=reverse($pep);
	chop($pep);chop($pep);$pep=reverse($pep);

	if (defined($mutpep{$pep})) 
	{
		# peptide pos: relative pos + offset (=max(0, mutation pos - 31))
		my $pos=$t[$#t-1];
		my $Protein=$t[1];
		my $crtPos=proPos($pos,$Protein);
		$t[$#t-1]=$crtPos;
		print join(";",@t);
	}
}
close IN;

#--------------------------------------------------------------------------------------
sub proPos
{
	my ($pos,$mut)=@_;
	$pos =~ m/^AA(\d+)toAA(\d+)$/; my ($aa1,$aa2)=($1,$2);
	# ATP5B|uc001slr.3|p.A170T|chr12:57037720-57037720.C_T
	my ($gene,$id,$proChg,$gPos)=split /\|/,$mut;
	my $vPos;
	# p.A170T   p.241_242del  
	if ($proChg =~ m/p\.(\d+)\_(\d+)del/)  {$vPos=$1;}
	elsif ($proChg =~ m/p\.[A-Z\-\*]+(\d+)[A-Z\-\*]+/) {$vPos=$1;}
	else {die "not recognized patter for proChg: $proChg\n";}
	#my $vPos=$1;
	my $offset=($vPos-31>=0)?($vPos-31):0;
	$aa1=$aa1+$offset;
	$aa2=$aa2+$offset;
	return join("",'AA',$aa1,'toAA',$aa2);
}
