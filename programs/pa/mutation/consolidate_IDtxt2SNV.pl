#!/usr/bin/perl

use strict;
use warnings;

if ($#ARGV!=3)
{
	die "perl consolidate_IDtxt2SNV.pl high_20.ref.fas high_20.mut.fas SNV_pep_ID2.txt out_prefix\n";
}

#-------------------------------------------------------------------------
my $snv_point_pos=31;
#-------------------------------------------------------------------------

my %refsq;
open(IN,$ARGV[0]) || die "cannot open $ARGV[0]\n";
while(<IN>)
{
	chomp;
	s/>//;
	my $id=$_;

	my $seq=<IN>;
	chomp($seq);
	$refsq{$id}=$seq;
}
close IN;

my %sq;
open(IN,$ARGV[1]) || die "cannot open $ARGV[1]\n";
while(<IN>)
{
	chomp;
	s/>//;
	my $id=$_;

	my $seq=<IN>;
	chomp($seq);
	$sq{$id}=$seq;
}
close IN;

my (%pephash,%snvhash);
open(IN,$ARGV[2]) || die "cannot open $ARGV[2]\n";
open(OUT2,">ID_posCheck.txt");
open(OUT,">$ARGV[3]_unexpectedCases.txt");
while(<IN>)
{
	chomp;
	if (/^Database/) { print OUT2 "$_\n"; next;}
	if (/^Peptide/) { print OUT2 "$_\n"; next;}
	#my @t=split(/\;/,$_);
	my ($pep,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos,$precursor_peak_intensity_percentage)=split(/\;/,$_);

	# starting from mut seq: record any pos that is different from ref (including fragment that not exist in ref) => @snv_point_pos
	$Protein =~ s/^>//;
	if (!defined($sq{$Protein})) { die "$Protein not exist in high_20.mut.fas!!!\n"; }
	if (!defined($refsq{$Protein})) { die "$Protein not exist in high_20.ref.fas!!!\n"; }
	my @snv_point_pos;
	varPos($sq{$Protein},$refsq{$Protein},\@snv_point_pos);

	# check if this peptide contains the SNV point (AA31)
	$pos = correctPos($pep,$sq{$Protein});
	print OUT2 join(";",$pep,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos,$precursor_peak_intensity_percentage),"\n";
	# AA42toAA50
	$pos =~ m/^AA(\d+)toAA(\d+)$/; my ($aa1,$aa2)=($1,$2);

	# overlap: [$aa1,$aa2] intersect @snv_point_pos
	my $snvPoint=checkVar($aa1,$aa2,\@snv_point_pos); 

	next unless ($snvPoint==1);

	#MAP9|uc011cin.2|p.K474R|chr4:156274377-156274377.T_C|SJRHB000026_X2-TB-13-2136
	my ($gene,$ucsc,$psnv,$snv,$sample)=split(/\|/,$Protein);
	if (!defined($sample)) { $sample='sample1'; }

	#my $pep=$t[0]; #K.GEALQAFEK.W
	chop($pep); chop($pep); $pep=reverse($pep); chop($pep); chop($pep); $pep=reverse($pep);

	$pephash{$pep}{outfile}{$Outfile}{XCorr}=$XCorr;
	$pephash{$pep}{outfile}{$Outfile}{dCn}=$dCn;
	$pephash{$pep}{outfile}{$Outfile}{ppm}=$ppm;

	if ( !defined($pephash{$pep}{bestoutfile}) or 
	$pephash{$pep}{outfile}{$pephash{$pep}{bestoutfile}}{XCorr}<$pephash{$pep}{outfile}{$Outfile}{XCorr} )
	{
		$pephash{$pep}{bestoutfile}=$Outfile;
	}

	$pephash{$pep}{SNV}{$snv}{gene}=$gene;
	$pephash{$pep}{SNV}{$snv}{psnv}=$psnv;
	$pephash{$pep}{SNV}{$snv}{sample}{$sample}='';
	
}
close IN;
close OUT;
close OUT2;

open(PEPPUB,">$ARGV[3]_peptides.txt");
print PEPPUB "peptide\tscan_counts\tbest_scan\tsearch_score\tdelta_score\tmutations\n";
open(PEPOUT,">$ARGV[3]_pep.txt");
print PEPOUT "peptide\tSC\tXCorr\tdCn\tppm\tSNV_count\tSNVs\n";
open(SCOUT,">$ARGV[3]_scan.txt");
print SCOUT "outfile\tpeptide\tXCorr\tdCn\tppm\tSNV_count\tSNVs\n";
foreach my $pep (keys %pephash)
{
	my $Outfile=$pephash{$pep}{bestoutfile};
	print PEPOUT "$pep\t",scalar(keys %{$pephash{$pep}{outfile}}),"\t$pephash{$pep}{outfile}{$Outfile}{XCorr}\t$pephash{$pep}{outfile}{$Outfile}{dCn}\t$pephash{$pep}{outfile}{$Outfile}{ppm}\t",scalar(keys %{$pephash{$pep}{SNV}}),"\t";
	my @mutations;
	foreach my $snv (keys %{$pephash{$pep}{SNV}}) { print PEPOUT "$snv,"; push @mutations,$snv;}
	print PEPOUT "\n";

	my @t=split /\//,$Outfile;
	my $bestoutfile=$t[$#t];
	print PEPPUB "$pep\t",scalar(keys %{$pephash{$pep}{outfile}}),"\t$bestoutfile\t$pephash{$pep}{outfile}{$Outfile}{XCorr}\t$pephash{$pep}{outfile}{$Outfile}{dCn}\t",join(",",@mutations),"\n";

	foreach $Outfile (keys %{$pephash{$pep}{outfile}})
	{
		print SCOUT "$Outfile\t$pep\t$pephash{$pep}{outfile}{$Outfile}{XCorr}\t$pephash{$pep}{outfile}{$Outfile}{dCn}\t$pephash{$pep}{outfile}{$Outfile}{ppm}\t",scalar(keys %{$pephash{$pep}{SNV}}),"\t";
		foreach my $snv (keys %{$pephash{$pep}{SNV}}) { print SCOUT "$snv,";}
		print SCOUT "\n";
	}

}
close PEPOUT;
close SCOUT;

# build %snvhash
foreach my $pep (keys %pephash)
{
	#next unless (scalar(keys %{$pephash{$pep}{SNV}})==1); # only use unique-mapped peptide
	my $fistEntry=1; # to deal with one peptide that maps to multiple locus
	foreach my $snv (keys %{$pephash{$pep}{SNV}})
	{
		if ( $fistEntry ) { $snvhash{$snv}{fistEntry}=1; $fistEntry=0; }
		$snvhash{$snv}{peptide}{$pep}='';
		$snvhash{$snv}{gene}=$pephash{$pep}{SNV}{$snv}{gene};
		foreach my $sample (keys %{$pephash{$pep}{SNV}{$snv}{sample}})
		{
			$snvhash{$snv}{sample}{$sample}='';
		}
	}
}

open(OUT1,">$ARGV[3]_mut_all.txt");
open(OUT2,">$ARGV[3]_mut_uniq.txt");
print OUT1 "mutation\tgene\tpeptideN\tscanN\tsampleN\tpeptides\tsamples\tvariantType\tinframe\n";
print OUT2 "mutation\tgene\tpeptideN\tscanN\tsampleN\tpeptides\tsamples\tvariantType\tinframe\n";
open(OUT3,">$ARGV[3]_events.txt");
print OUT3 "mutation\tgene\taa_change\tpeptide_counts\tpeptide_sequences\tscan_counts\tbest_scans\n";
foreach my $snv (keys %snvhash)
{
	my ($variantType,$inframe)=varType($snv);
	my $SC=0;
	my @best_scans;
	my $aa_change='';
	foreach my $pep (keys %{$snvhash{$snv}{peptide}})
	{
		$SC+=scalar(keys %{$pephash{$pep}{outfile}});
		my @t=split /\//,$pephash{$pep}{bestoutfile};
		push @best_scans, $t[$#t];
		$aa_change=$pephash{$pep}{SNV}{$snv}{psnv};
	}
	print OUT1 "$snv\t$snvhash{$snv}{gene}\t",scalar(keys %{$snvhash{$snv}{peptide}}),"\t$SC\t",scalar(keys %{$snvhash{$snv}{sample}}),"\t";
	foreach my $pep (keys %{$snvhash{$snv}{peptide}}) { print OUT1 "$pep,"; }
	print OUT1 "\t";
	foreach my $sample (keys %{$snvhash{$snv}{sample}}) { print OUT1 "$sample,"; }
	print OUT1 "\t$variantType\t$inframe\n";

	print OUT3 "$snv\t$snvhash{$snv}{gene}\t$aa_change\t",scalar(keys %{$snvhash{$snv}{peptide}}),"\t",join(",",(keys %{$snvhash{$snv}{peptide}})),"\t$SC\t",join(",",@best_scans),"\n";

	# print unique mut
	if ( $snvhash{$snv}{fistEntry} )
	{
		print OUT2 "$snv\t$snvhash{$snv}{gene}\t",scalar(keys %{$snvhash{$snv}{peptide}}),"\t$SC\t",scalar(keys %{$snvhash{$snv}{sample}}),"\t";
		foreach my $pep (keys %{$snvhash{$snv}{peptide}}) { print OUT2 "$pep,"; }
		print OUT2 "\t";
		foreach my $sample (keys %{$snvhash{$snv}{sample}}) { print OUT2 "$sample,"; }
		print OUT2 "\t$variantType\t$inframe\n";
	}
}
close OUT;

#------------------------------------------------------------
sub varType
{
	my ($snv)=@_;

	my ($pos,$alleles)=split /\./,$snv;
	my ($ref,$mut)=split /\_/,$alleles;

	my $variantType='snv';
	my $inframe='yes';
	if ($ref eq '-') 
	{ 
		$variantType='ins'; 
		if (length($mut) % 3 == 0) {}
		else { $inframe='no'; }
	}
	elsif ($mut eq '-') 
	{ 
		$variantType='del'; 
		if (length($ref) % 3 == 0) {}
		else { $inframe='no'; }
	}
	return ($variantType,$inframe);
}

sub varPos
{
	my ($mut,$ref,$snv_point_pos)=@_;

	my @m=split //,$mut;
	my @r=split //,$ref;

	for (my $i=0;$i<=$#m; $i++)
	{
		if (defined($r[$i]) and $r[$i] eq $m[$i]) {}
		else { push @{$snv_point_pos},$i+1; }
	}
}

sub checkVar
{
	my ($aa1,$aa2,$snv_point_pos)=@_;

	my $v=0; # contains variant site?
	foreach my $p (@{$snv_point_pos})
	{
		if ( $aa1<=$p and $p<=$aa2 ) { $v=1; last }
	}
	return $v;
}

sub correctPos
{
	my ($pepseq,$proseq)=@_;

	chop($pepseq);chop($pepseq);$pepseq=reverse($pepseq);
	chop($pepseq);chop($pepseq);$pepseq=reverse($pepseq);
	$pepseq =~ s/[\*\#\@\%\&\~\$\^]//g;
	my $fpos = index($proseq, $pepseq) + 1;
	my $lpos = $fpos + length($pepseq) - 1;
	my $position = "AA$fpos"."to"."AA$lpos";
	return $position;
}
