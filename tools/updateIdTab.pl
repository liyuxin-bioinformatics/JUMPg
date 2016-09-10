#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=2)
{
        die "perl updateIdTab.pl hg19_kgXref.txt uniProt2015toUCSChg19_v1.0.txt\n";
}

#------------------------------------------------------
my $uniprot2ucsc=$ARGV[1];
my $kgXref=$ARGV[0];
#------------------------------------------------------

# parse uniProt2015toUCSChg19_v1.0.txt
# only use 100%
# build %ucsc2uniprot{$ucsc}=$prot
my (%ucsc2uniprot);
open(IN,$uniprot2ucsc) or die "Cannot open file $uniprot2ucsc!!!\n";
while(<IN>) {
	chomp;
	next if (/^UniProt_ID/);
	my ($uniprot,$ucsc,$pct)=split(/\t/,$_);
	next unless ($pct>=100);

	#$ucsc2uniprot{$ucsc}=$uniprot;
	$ucsc2uniprot{$ucsc}=(split(/\|/,$uniprot))[1];
}
close IN;

# parse hg19_kgXref.txt
# next if (!defined($ucsc2uniprot{$ucsc}))
# convert ID
open(IN,$kgXref) or die "Cannot open file $kgXref!!!\n";
open(OUT,">$kgXref\_updated");
while(<IN>) {
	chomp;
	my @t=split(/\t/,$_);
	next if (!defined($ucsc2uniprot{$t[0]}));
	$t[2]=$ucsc2uniprot{$t[0]};
	print OUT join("\t",@t),"\n";
}
close IN;
close OUT;
