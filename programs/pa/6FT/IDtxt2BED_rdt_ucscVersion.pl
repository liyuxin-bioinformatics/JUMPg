#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=3)
{
	die "perl IDtxt2BED_rdt_ucscVersion.pl reads.alg.out  trackName out\n";
}

my $trackName=$ARGV[1];

# build %cdshash
my (%cdshash,%rnahash);
#open(IN,'/home/yli4/annotations/refSeq_mRNA2protein_072514.txt') or die "Cannot open annotatiion file!!!\n";
open(IN,'/home/yli4/annotations/hg19_knownGenes_112713.txt') or die "Cannot open annotatiion file!!!\n";
while(<IN>)
{
	next if (/^#/);
	chomp;
	my @t=split(/\t/,$_);
	my ($nm,$np)=($t[0],$t[0]);

	if (!defined($cdshash{$np}))
	{
		$cdshash{$np}{nm}=$nm;
		$cdshash{$np}{strand}=$t[2];
		$cdshash{$np}{chr}=$t[1];
		$cdshash{$np}{line}=$_;

		my %hash=%{$cdshash{$np}};
		my $ck=buildCDShash($cdshash{$np});
		$cdshash{$np}{check}=$ck;
	}

	if (!defined($rnahash{$np}))
	{
		$rnahash{$np}{nm}=$nm;
		$rnahash{$np}{strand}=$t[2];
		$rnahash{$np}{chr}=$t[1];
		$rnahash{$np}{line}=$_;

		my %hash=%{$rnahash{$np}};
		my $ck=buildRNAhash($rnahash{$np});
		$rnahash{$np}{check}=$ck;
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
open(OUT2,">$ARGV[0].updated");
open(OUT,">$ARGV[2]");
print OUT 'track name=',$trackName,' description="',$trackName,'" visibility=3 itemRgb="On"  color="255,0,247"',"\n";
close OUT;
while(<IN>)
{
        chomp;
	if (/^outfile/) { print OUT2 "$_\tsType\n";  next; }

	#print OUT2 "$_\t"; 
	my $sType='';

        my @t=split(/\t/,$_);

	my ($outfile,$peptide,$level,$np,$b5,$b3,$strand)=($t[0],$t[3],$t[14],$t[7],$t[18],$t[19],$t[20]);
	next unless (defined($b5) and defined($b3) and defined($strand) and $b5 =~ m/^\d+$/ and $b3 =~ m/^\d+$/ and $b5 ne '' and $b3 ne '');

	if ($level==1) # CDS
	{
		$sType='CDS';
		next unless (defined($cdshash{$np}));
		
		my ($pepb5_ref,$pepb3_ref)=genomicConvert($np,$b5,$b3,\%cdshash);
		if (!defined($pepb5_ref) or !defined($pepb3_ref)) {print $_,"\n";next;}

		printBED($cdshash{$np},$pepb5_ref,$pepb3_ref,$peptide,$outfile,$ARGV[2],$strand);
	}
	elsif ($level==2) # RNA
	{
		next unless (defined($rnahash{$np}));
		$sType=$rnahash{$np}{RNAtype};

		my ($pepb5_ref,$pepb3_ref)=genomicConvert($np,$b5,$b3,\%rnahash);
		if (!defined($pepb5_ref) or !defined($pepb3_ref)) {print $_,"\n";next;}

		printBED($rnahash{$np},$pepb5_ref,$pepb3_ref,$peptide,$outfile,$ARGV[2],$strand);
	}
	elsif ($level==5) # genome
	{
		$sType='genome';
		printBEDgenome($outfile,$peptide,$np,$b5,$b3,$strand,$ARGV[2]);
	}
	#else { die "unexpected level: $level\n"; }
	else { next; }

	#print OUT2 "$sType\n";
	print OUT2 join("\t",@t),"\t$sType\n";
}
close IN;
close OUT;
close OUT2;

#------------------------------------------------------------------------------
sub printBEDgenome
{
	my ($outfile,$peptide,$np,$b5,$b3,$strand,$output)=@_;

	open(OUT,">>$output");
	my $fStart=$b5;
	print OUT "$np\t$fStart\t$b3\t$peptide\t900\t$strand\t$fStart\t$b3\t0\t1\t",$b3-$fStart,",\t0,\t$outfile\n";
	close OUT;
}

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
	my ($prohash,$pepb5_ref,$pepb3_ref,$peptide,$outfile,$output,$strand)=@_;

	my @pepb5=@$pepb5_ref;
	my @pepb3=@$pepb3_ref;

	# determine strand
	my $finalStrand=$$prohash{strand};
	if ($strand eq '-')
	{
		$finalStrand=($$prohash{strand} eq '+')?'-':'+';
	}

	open(OUT,">>$output");
	if (!defined($$prohash{chr})) {die "$peptide,$outfile,$output\n";}
	#print OUT "$$prohash{chr}\t$pepb5[0]\t$pepb3[$#pepb3]\t$peptide\t900\t$$prohash{strand}\t$pepb5[0]\t$pepb3[$#pepb3]\t0\t",$#pepb3+1,"\t";
	print OUT "$$prohash{chr}\t$pepb5[0]\t$pepb3[$#pepb3]\t$peptide\t900\t$finalStrand\t$pepb5[0]\t$pepb3[$#pepb3]\t0\t",$#pepb3+1,"\t";
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

sub buildCDShash
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

sub buildRNAhash
{
	my ($hash)=@_;

	my @t=split(/\t/,$$hash{line});
	my $cdsStart=$t[3]; # txStart
	my $cdsEnd=$t[4];  # txEnd

	if ($t[5]==$t[6]) # non-coding genes
	{ $$hash{RNAtype}='ncGene'; }
	else { $$hash{RNAtype}='mRNA'; }

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
	return 1;
	#if ( $cdsL % 3 == 0 )  { return 1;}
	#else {return -1;}
}
