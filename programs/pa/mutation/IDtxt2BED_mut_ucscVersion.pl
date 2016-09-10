#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=3)
{
	die "perl IDtxt2BED.pl ID.txt trackName out\n";
}

my $trackName=$ARGV[1];

# build %prohash
my (%prohash);
#open(IN,'/home/yli4/annotations/refSeq_mRNA2protein_072514.txt') or die "Cannot open annotatiion file!!!\n";
open(IN,'/home/yli4/annotations/hg19_knownGenes_112713.txt') or die "Cannot open annotatiion file!!!\n";
while(<IN>)
{
	next if (/^#/);
	chomp;
	my @t=split(/\t/,$_);
	my ($nm,$np)=($t[0],$t[0]);

	if (!defined($prohash{$np}))
	{
	$prohash{$np}{nm}=$nm;
	$prohash{$np}{strand}=$t[2];
	$prohash{$np}{chr}=$t[1];
	$prohash{$np}{line}=$_;

	my %hash=%{$prohash{$np}};
	#my $ck=checkFrame(\%hash);
	my $ck=checkFrame($prohash{$np});
	$prohash{$np}{check}=$ck;
	}
	#print "$nm\t$ck\n";

=head
	# check wrong frame examples
	if ($ck==1) 
	{ 
		print "$nm\t$np\t"; 
		for (my $i=0;$i<scalar(@{$hash{eb3}}); $i++) {print "$hash{eb5}[$i],";}
		print "\t";
		for (my $i=0;$i<scalar(@{$hash{eb3}}); $i++) {print "$hash{eb3}[$i],";}
		print "\t";	for (my $i=0;$i<scalar(@{$hash{eb3}}); $i++) {print "$hash{elength}[$i],";}
		print "\t";	for (my $i=0;$i<scalar(@{$hash{eb3}}); $i++) {print "$hash{cmllength}[$i],";}
		print "\n";
	}
=cut
}
close IN;

# parse ID.txt
open(IN,"$ARGV[0]");
open(OUT,">$ARGV[2]");
print OUT 'track name=',$trackName,' description="',$trackName,'" visibility=3 itemRgb="On" color="255,0,247"',"\n";
close OUT;
while(<IN>)
{
	if (/^Database/) {  next; }
        if (/^Peptide/) {  next; }

        chomp;
        my @t=split(/\;/,$_);

	my ($mut,$pos,$peptide,$outfile)=($t[1],$t[$#t-1],$t[0],$t[2]);

	# ATP5B|uc001slr.3|p.A170T|chr12:57037720-57037720.C_T
	my ($gene,$np,$pMut,$mutGenPos)=split(/\|/,$mut);
	next unless (defined($prohash{$np}));

	# relative positions
	# AA351toAA358
	if ($pos =~ m/^AA(\d+)toAA(\d+)$/) {} else {die "wrong pos format: $pos!!!\ndifferent ID.txt format?\n";}
        $pos =~ m/^AA(\d+)toAA(\d+)$/;
        my ($aa1,$aa2)=($1,$2);
        my ($b5,$b3)=($aa1*3-3,$aa2*3);

	# 3 classes
	my ($variantType,$inframe,$indelLength)=varType($mutGenPos);
	if ($variantType eq 'snv')
	{
		my ($pepb5_ref,$pepb3_ref)=genomicConvert($np,$b5,$b3,\%prohash);
		if (!defined($pepb5_ref) or !defined($pepb3_ref)) {print $_,"\n";next;}

		printBED($prohash{$np},$pepb5_ref,$pepb3_ref,$peptide,$outfile,$ARGV[2]);
	}
	elsif ($variantType eq 'ins')
	{
		# del position on protein
		my $delPos;
		# p.185_185del
		if ( $pMut =~ m/^p\.(\d+)\_\d+del$/ ) { $delPos=$1; }
		# p.K260R
		elsif ( $pMut =~ m/^p\.[A-Z](\d+)[A-Z]$/ ) { $delPos=$1; }
		else {die "Wrong ins format: $mut,$pMut\n";}
		# del in peptide or before peptide?
		my $delInPep=0;
		if ($aa1<=$delPos and $delPos<=$aa2) {$delInPep=1;} # del in peptide
		else {$delInPep=0;} # del before peptide (that causes frame shift)

		if ($delInPep) {}
		else { $b5=$b5-$indelLength; } # del before peptide (that causes frame shift)

		$b3=$b3-$indelLength;
		my ($pepb5_ref,$pepb3_ref)=genomicConvert($np,$b5,$b3,\%prohash);
		if (!defined($pepb5_ref) or !defined($pepb3_ref)) {print $_,"\n";next;}
		printBED($prohash{$np},$pepb5_ref,$pepb3_ref,$peptide,$outfile,$ARGV[2]);
	}
	elsif ($variantType eq 'del')
	{
		# del position on protein
		my $delPos;
		# p.185_185del
		if ( $pMut =~ m/^p\.(\d+)\_\d+del$/ ) { $delPos=$1; }
		# p.K260R
		elsif ( $pMut =~ m/^p\.[A-Z](\d+)[A-Z]$/ ) { $delPos=$1; }
		else {die "Wrong del format: $mut,$pMut\n";}
		# del in peptide or before peptide?
		my $delInPep=0;
		if ($aa1<=$delPos and $delPos<=$aa2) {$delInPep=1;} # del in peptide
		else {$delInPep=0;} # del before peptide (that causes frame shift)

		if ($delInPep)
		{
			# boundary definitions
			my $frg1_b5=$b5;
			my $frg1_b3=($delPos-1)*3;

			my $frg2_b5=$frg1_b3+$indelLength;
			my $frg2_b3=$b3+$indelLength;

			# BED convert for each fragment 
			my ($frg1_pepb5,$frg1_pepb3)=genomicConvert($np,$frg1_b5,$frg1_b3,\%prohash);
			my ($frg2_pepb5,$frg2_pepb3)=genomicConvert($np,$frg2_b5,$frg2_b3,\%prohash);
			if (!defined($frg1_pepb5) or !defined($frg1_pepb3)) {print $_,"\n";next;}
			if (!defined($frg2_pepb5) or !defined($frg2_pepb3)) {print $_,"\n";next;}

			# merge 2 fragments and print
			my ($pepb5_ref,$pepb3_ref)=mergeFragments($frg1_pepb5,$frg1_pepb3,$frg2_pepb5,$frg2_pepb3,$prohash{$np});
			printBED($prohash{$np},$pepb5_ref,$pepb3_ref,$peptide,$outfile,$ARGV[2]);
		}
		else
		{
			$b5+=$indelLength;
			$b3+=$indelLength;
			my ($pepb5_ref,$pepb3_ref)=genomicConvert($np,$b5,$b3,\%prohash);
			if (!defined($pepb5_ref) or !defined($pepb3_ref)) {print $_,"\n";next;}
			printBED($prohash{$np},$pepb5_ref,$pepb3_ref,$peptide,$outfile,$ARGV[2]);
		}
	}
}
close IN;

#------------------------------------------------------------------------------
sub mergeFragments
{
	my ($frg1_pepb5,$frg1_pepb3,$frg2_pepb5,$frg2_pepb3,$hash)=@_;

	my (@pepb5,@pepb3);
	if ($$hash{strand} eq '+')
	{
		@pepb5=(@$frg1_pepb5,@$frg2_pepb5);
		@pepb3=(@$frg1_pepb3,@$frg2_pepb3);
	}
	else
	{
		@pepb5=(@$frg2_pepb5,@$frg1_pepb5);
		@pepb3=(@$frg2_pepb3,@$frg1_pepb3);
	}

	return (\@pepb5,\@pepb3);
}

sub printBED
{
	my ($prohash,$pepb5_ref,$pepb3_ref,$peptide,$outfile,$output)=@_;

	my @pepb5=@$pepb5_ref;
	my @pepb3=@$pepb3_ref;

	open(OUT,">>$output");
	if (!defined($$prohash{chr})) {die "$peptide,$outfile,$output\n";}
	print OUT "$$prohash{chr}\t$pepb5[0]\t$pepb3[$#pepb3]\t$peptide\t900\t$$prohash{strand}\t$pepb5[0]\t$pepb3[$#pepb3]\t0\t",$#pepb3+1,"\t";
                for (my $i=0; $i<=$#pepb3; $i++) { print OUT $pepb3[$i]-$pepb5[$i],","; }
                print OUT "\t";
                for (my $i=0; $i<=$#pepb3; $i++) { print OUT $pepb5[$i]-$pepb5[0],","; }
                print OUT "\t$outfile";
                print OUT "\n";

	close OUT;
}

sub genomicConvert
{
	my ($np,$b5,$b3,$prohash)=@_;

                if ($$prohash{$np}{strand} eq '-')
                {
                        my $t3=$$prohash{$np}{cdslength}-$b5;
                        my $t5=$$prohash{$np}{cdslength}-$b3;
                        $b5=$t5; $b3=$t3;
                }

	my (@pepb5,@pepb3);
	my $i=0; if (!defined($$prohash{$np}{cmllength}[$i])) {die "$np,$b5,$b3,$i,$$prohash{$np}{cmllength}[$i]\n";}
	if ($b5>$$prohash{$np}{cmllength}[scalar(@{$$prohash{$np}{cmllength}})-1]) {print "peptide b5 is downstream of last CDS exon!!!\n$np,$b5,$$prohash{$np}{cmllength}[scalar(@{$$prohash{$np}{cmllength}})-1]\n"; return;}
	while ( $b5>$$prohash{$np}{cmllength}[$i] ) {$i++;}

	my $j=0; $pepb5[$j]=$$prohash{$np}{eb3}[$i]-($$prohash{$np}{cmllength}[$i]-$b5);
	if ($b3>$$prohash{$np}{cmllength}[scalar(@{$$prohash{$np}{cmllength}})-1]) {print "peptide b3 is downstream of last CDS exon!!!\n$np,$b3,$$prohash{$np}{cmllength}[scalar(@{$$prohash{$np}{cmllength}})-1]\n";return;}
	while ( $b3>$$prohash{$np}{cmllength}[$i] )
                {
                        $pepb3[$j]=$$prohash{$np}{eb3}[$i];
                        $i++;
                        $j++;
                        $pepb5[$j]=$$prohash{$np}{eb5}[$i];
                }
                $pepb3[$j]=$$prohash{$np}{eb3}[$i]-($$prohash{$np}{cmllength}[$i]-$b3);

	# check whether the 1st 'exon' is empty
	if ( $pepb5[0] == $pepb3[0] )
	{
		shift(@pepb5);
		shift(@pepb3);
	}	

	return (\@pepb5,\@pepb3);
}

sub varType
{
        my ($snv)=@_;

        my ($pos,$alleles)=split /\./,$snv;
        my ($ref,$mut)=split /\_/,$alleles;

        my $variantType='snv';
        my $inframe='yes';
	my $indelLength=0;
        if ($ref eq '-')
        {
		# -_CC
                $variantType='ins';
                if (length($mut) % 3 == 0) {}
                else { $inframe='no'; }
		$indelLength=length($mut);
        }
        elsif ($mut eq '-')
        {
                $variantType='del';
                if (length($ref) % 3 == 0) {}
                else { $inframe='no'; }
		$indelLength=length($ref);
        }
        return ($variantType,$inframe,$indelLength);
}

sub checkFrame
{
	my ($hash)=@_;

	my @t=split(/\t/,$$hash{line});
	#if ($t[6]==$t[7]) {return 0;} # non-coding genes
	my $cdsStart=$t[5];
	my $cdsEnd=$t[6];
	if ($cdsStart==$cdsEnd) {return 0;} # non-coding genes

	my @b5=split(/,/,$t[8]); # hg19.knownGene.exonStarts
	my @b3=split(/,/,$t[9]); # hg19.knownGene.exonEnds

	my @exonl;
	for (my $i=0;$i<=$#b5; $i++) { $exonl[$i]=$b3[$i]-$b5[$i]; }
	
	my $i=0;
	while ($cdsStart>$b3[$i]) {$i++;} my $j=0; $$hash{eb5}[$j]=$cdsStart;
	while ($cdsEnd>$b3[$i])
	{
		$$hash{eb3}[$j]=$b3[$i];
		$i++;
		$j++;
		$$hash{eb5}[$j]=$b5[$i];
	}
	$$hash{eb3}[$j]=$cdsEnd;

	for ($i=0;$i<scalar(@{$$hash{eb3}}); $i++) { $$hash{elength}[$i]=$$hash{eb3}[$i]-$$hash{eb5}[$i]; }
	my $cdsL=0;
	for ($i=0;$i<scalar(@{$$hash{eb3}}); $i++) { $cdsL+=$$hash{elength}[$i]; $$hash{cmllength}[$i]=$cdsL; }
	$$hash{cdslength}=$cdsL;
	if ( $cdsL % 3 == 0 )  { return 1;}
	else {return -1;}
}
