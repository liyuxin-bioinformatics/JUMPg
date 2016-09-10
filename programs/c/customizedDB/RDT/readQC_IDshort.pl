#!/usr/bin/perl -w -I /home/yli4/lib/perl5

use strict;
use FastaqParser;

my $fq=FastaqParser->new();


if (scalar(@ARGV)!=3)
{
        die "Usage: perl readQC.pl in.fq phredSystem['sanger(32 based, i.e., '#')'/'solexa(64 bsaed, i.e. 'B')'] out\n";
}

######################
my $minReadLength=30;
my $phredSystem=$ARGV[1];
#####################

my (@seq);
my $input=$ARGV[0];

open(IN,"$input") or die "cannot open $input\n";
open(OUT,">$ARGV[2]");

while(<IN>)
{
        $seq[0]=$_; chomp($seq[0]);
        $seq[1]=<IN>; chomp($seq[1]);
        $seq[2]=<IN>; chomp($seq[2]);
        $seq[3]=<IN>; chomp($seq[3]);

	#my @t0=$fq->trim3PrimeAdaptor(\@seq);
        my @t=$fq->trim3PrimeLowQuality(\@seq,$phredSystem);

        next if (length($t[1])<$minReadLength); #length control

	#read ID shorten
	#my $id =$seq[0];
	#$id =~  m/^\@\w+\:(.*?)\#\w+\/(\d)$/;
	#$id="$1\/$2";
	#$id =~ s/\:/_/g; $id =~ s/\//-/;

	#@HWI-ST1200:136:C2R9HACXX:7:1101:10000:13117/1
	$seq[0] =~ s/^\@//;
	$seq[0] =~ s/\s/\_/;
	#my $id="$tmp[3]_$tmp[4]_$tmp[5]_$1-$2";
	my $id=$seq[0];

	#print OUT "\@$id\n$t[1]\n$t[2]\n$t[3]\n";
	print OUT "\>$id\n$t[1]\n";
}

close IN;
close OUT;
