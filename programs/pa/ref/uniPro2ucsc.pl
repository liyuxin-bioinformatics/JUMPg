#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=2)
{
	die "Usage: perl uniPro2refSeq.pl ID.txt uniProt2015toUCSChg19_v1.0.txt\n";
}

#----------------------------------------------------------------
my $idtxt=$ARGV[0];
my $up2ucsc=$ARGV[1];
#----------------------------------------------------------------

my (%unipro2gene,%prohash);
open(IN,$up2ucsc);
my $line=<IN>;
while(<IN>)
{
	chomp;
	next if (/^#/);
	my @t=split(/\t/,$_);

	#next unless ($t[2]>=100); # 100/% identify; from Suiping's blast alignment 
	#$unipro2gene{$t[0]}=$t[1]; # from Suiping's blast alignment

	next if ($t[2] eq ''); # for hg19_kgXref.txt downloaded by AnnoVar
	$unipro2gene{$t[2]}=$t[0]; # from hg19_kgXref.txt downloaded by AnnoVar
}
close IN;

open(IN,$idtxt);
$line=<IN>; print  $line;
$line=<IN>; print  $line;
while(<IN>)
{
	chomp; $line=$_;
	my @t=split(/\;/,$_);
	next if ($t[1] =~ /^co\|/); # rm contaminants
	my ($type,$unipro,$name)=split(/\|/,$t[1]);
	#my $unipro=$t[1];
	if (defined($unipro2gene{$unipro})) 
	{ 
		$t[1]=$unipro2gene{$unipro}; 
		print join(";",@t),"\n";
	}
}
close IN;
