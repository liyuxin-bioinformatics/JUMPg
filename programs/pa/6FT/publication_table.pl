#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=3)
{
        die "usage: perl publication_table.pl refFlat BED_file alg.updated2.pep\n";
}

my %ucsc2gene;
open(IN,"$ARGV[0]");
while(<IN>)
{
        chomp;
	my @t=split /\t/,$_;
	if (!defined($ucsc2gene{$t[0]}) or $ucsc2gene{$t[0]} eq '')
	{
		#$ucsc2gene{$t[0]}=$t[7]; # /home/yli4/annotations/knownGenes_uniPro_012314.txt
		$ucsc2gene{$t[0]}=$t[4]; # hg19_kgXref.txt downloaded by AnnoVar
	}
}
close IN;

open(IN,"$ARGV[1]"); # mouseBrain_RDT_peptides.sc.bed.check
my %pep2genome;
my $line=<IN>;
while(<IN>)
{
        my @t=split /\t/,$_;
        $t[3] =~ s/\[\d+\]$//;;
        $pep2genome{$t[3]}{strand}=$t[5];
        $pep2genome{$t[3]}{chr}=$t[0];
        $pep2genome{$t[3]}{start}=$t[1];
        $pep2genome{$t[3]}{end}=$t[2];
}
close IN;

open(IN,"$ARGV[2]"); # mouseBrain_RDT_test1_reads.out3.alg.updated2.pep.check
$line=<IN>; #chomp($line);
print "peptide\tscan_counts\tbest_scan\tsearch_score\tdelta_score\tread_counts\tgene\tlocation_type\tframe\tgenome_chr\tgenome_strand\tgenome_start\tgenome_end\n";
while(<IN>)
{
        chomp;
        my @t=split /\t/,$_;
	my $gene='';
	if (defined($ucsc2gene{$t[7]})) { $gene=$ucsc2gene{$t[7]}; }
        if (defined($pep2genome{$t[3]}))
        {
                print "$t[3]\t$t[4]\t$t[0]\t$t[1]\t$t[2]\t$t[5]\t$gene\t$t[22]\t$t[21]\t$pep2genome{$t[3]}{chr}\t$pep2genome{$t[3]}{strand}\t$pep2genome{$t[3]}{start}\t$pep2genome{$t[3]}{end}\n";
        }
        else { die "$t[3] not found in $ARGV[1]\n"; }
}

close IN;

