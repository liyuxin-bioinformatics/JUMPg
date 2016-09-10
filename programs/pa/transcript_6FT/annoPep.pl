#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/g

use warnings;
use strict;
use Getopt::Long;
use IDtxt_parser;
use PrimarySeq;
use rnaPepDecorator;
#use FindBin;
#use lib $FindBin::Bin;
#use File::Basename;
#use lib dirname (__FILE__);

# module instances
my $idp=IDtxt_parser->new;
my $rpd=rnaPepDecorator->new;
my $ps=PrimarySeq->new;

# input parameters

my @params = (
"-idtxt=s",
"-reference-protein=s",
"-path-to-blat=s",
"-rna=s",
"-reference-genome=s",
"-refFlat=s",
);

my %FLAGS;
GetOptions(\%FLAGS, @params);

my $flag;
$flag="idtxt";
if ($FLAGS{$flag}){
	unless (-e $FLAGS{$flag}) {
		die "Error: $FLAGS{$flag} not exist!!!\n";
	}
} else {
        die "Error: specify -$flag\n";
}

$flag="rna";
if ($FLAGS{$flag}){
	unless (-e $FLAGS{$flag}) {
		die "Error: $FLAGS{$flag} not exist!!!\n";
	}
} else {
        die "Error: specify -$flag\n";
}

$flag="path-to-blat";
if ($FLAGS{$flag}){
	unless (-e $FLAGS{$flag}) {
		die "Error: $FLAGS{$flag} not exist!!!\n";
	}
} else {
        die "Error: specify -$flag\n";
}

$flag="reference-genome";
if ($FLAGS{$flag}){
	unless (-e $FLAGS{$flag}) {
		die "Error: $FLAGS{$flag} not exist!!!\n";
	}
} else {
        die "Error: specify -$flag\n";
}

$flag="refFlat";
if ($FLAGS{$flag}){
	unless (-e $FLAGS{$flag}) {
		die "Error: $FLAGS{$flag} not exist!!!\n";
	}
} else {
        die "Error: specify -$flag\n";
}

# working starts here
# IDtxt_parser
my ($scanhash,$peptidehash,$proteinhash)=$idp->parse_IDtxt($FLAGS{"idtxt"});

# filter peptides if exist on one reference protein file?
$flag="reference-protein";
if ($FLAGS{$flag}){
	unless (-e $FLAGS{$flag}) {
		die "Error: $FLAGS{$flag} not exist!!!\n";
	}

	# read reference protein file
	my $refproseq=$ps->parseFas($FLAGS{"reference-protein"});
	
	# rmKnownPep: filter peptides if exist on one reference protein file
	my $novelscanhash=$rpd->rmKnownPep($scanhash,$refproseq);
	$idp->printIDtxt($novelscanhash,'novelscan.txt');
	$scanhash=$novelscanhash
}

# build %scanhash{$out}{RNA}{$rnaID}{frame/ORFstart/ORFend/seq/pep2ORFstart/pep2ORFend/pep2rnaStart/pep2rnaEnd}
$rpd->buildRNAsubHash($scanhash);

# getRNA: extract RNA sequences for the accepted peptides
my $rnaseqhash=$rpd->getRNA($scanhash,$FLAGS{"rna"});
my $rnaSeqFastaFile='blat_input.fa';
$ps->printFas($rnaseqhash,$rnaSeqFastaFile);

# alignment: blat
my $aligmentFile='blat_output.psl';
$rpd->runBlat($FLAGS{"path-to-blat"},$FLAGS{"reference-genome"},$rnaSeqFastaFile,$aligmentFile);

# parseBlast: RNA-to-genome position
my $rnaalignhash=$rpd->parseBlat($aligmentFile);

# pepGenomePos: peptide-to-genome postion
# = RNA-to-genome start + peptide-to-RNA position 
$rpd->pepGenomePos($scanhash,$rnaalignhash);
my $bedfile='output.bed';
$rpd->printBED($scanhash,$bedfile);

# parseRefFlat: gene-to-genome position
my $transcripthash=$rpd->parseRefFlat($FLAGS{"refFlat"});

# pep2gene: assign gene to a peptide
$rpd->pep2gene($scanhash,$transcripthash);

# pepFunRegion: assign functional region (CDS, intron etc.) to peptide
$rpd->pepFunRegion($scanhash,$transcripthash);

# print results
$rpd->publicTab($scanhash,$peptidehash);
$rpd->publicTabPSM($scanhash);
print "Done\n";
