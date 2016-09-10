#!/usr/bin/perl -w

use strict;

if (scalar(@ARGV)!=3)
{
	die "Usage: perl seqpad2fas.pl file.seqpad  mut/ref out.fas\n";
}

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[2]");

my $type=$ARGV[1];

while (<IN>)
{
	chomp;
	my @tmp=split(/\t/,$_);

	#check whether AA is changed
	#p.P78P
	if ($tmp[3] =~ m/del/ || $tmp[3] =~ m/ins/) {}
	elsif ($tmp[6] =~ m/^p\.([\w\*])\d+([\w\*])/) 
	{
		next if ( $1 eq $2 );
	}
=head	#p.60_62del
	if ($tmp[6] =~ m/^p\.(\d+)\_(\d+)del/)
	{
		$tmp[6] = "p\.$1d";
		$tmp[6] .= $2-$1+1;
	}
=cut

	# genomic coordinates
	my $ps="$tmp[9]:$tmp[10]-$tmp[11]\.$tmp[12]_$tmp[13]";

	# print out
	#next if ($tmp[8] =~ /\*/);	# ignore if stop codon exist
	#print OUT '>',$tmp[1],'|',$tmp[3],'|',$tmp[6],"\n",$tmp[8],"\n";
	if ( $type eq 'mut' ) {print OUT '>',$tmp[0],'|',$tmp[1],'|',$tmp[6],'|',$ps,"\n",$tmp[8],"\n";}
	elsif ( $type eq 'ref' ) {print OUT '>',$tmp[0],'|',$tmp[1],'|',$tmp[6],'|',$ps,"\n",$tmp[7],"\n";}
}

close IN;
close OUT;
