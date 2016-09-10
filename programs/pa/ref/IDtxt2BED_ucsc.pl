#!/usr/bin/perl -I /home/yli4/lib/perl5/

use strict;
use warnings;
use  PrimarySeq;

if (scalar(@ARGV)!=4)
{
	#die "perl IDtxt2BED.pl ID.txt ucscRefFlat ucscProtein.fas trackName \n";
	die "perl IDtxt2BED.pl ID.txt ucscRefFlat trackName output\n";
}

my $idtxt=$ARGV[0];
my $refFlat=$ARGV[1];
#my $proteinSeq=$ARGV[2];
my $trackName=$ARGV[2];
my $output=$ARGV[3];

# build %prohash
my (%prohash,$proseq);
my $ps=PrimarySeq->new();
#$proseq=$ps->parseFas($proteinSeq);
open(IN,$refFlat) or die "Cannot open annotatiion file!!!\n";
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
		my $ck=checkFrame($prohash{$np});
		$prohash{$np}{check}=$ck;
=head
		if (defined($proseq->{$np})) {
			my $proteinLength=length($proseq->{$np});
			if ( $proteinLength*3 + 3 != $prohash{$np}{cdslength} ) {
				$prohash{$np}{check}=0;
				warn "WARNING: protein length ($proteinLength) and annotation (cdslength: $prohash{$np}{cdslength}) unmatched for $np! \n";
			}
		} else {
			$prohash{$np}{check}=0;
		}
=cut
	}
}
close IN;

# parse ID.txt and map to genomic positions
open(IN,$idtxt);
open(OUT,">$output");
print OUT 'track name=',$trackName,' description="',$trackName,'" visibility=3 itemRgb="On" color="0,128,0"',"\n";
while(<IN>)
{
        if (/^Database/) {  next; }
        if (/^Peptide/) {  next; }

        chomp;
        my @t=split(/\;/,$_);

	my ($np,$pos,$peptide,$outfile)=($t[1],$t[$#t-1],$t[0],$t[2]);
	next if ($np =~ /\+/ || /CON_/);
	#$np =~ s/\.(\d+)$//;
#=head
	if (defined($prohash{$np}))
	{
		next unless ($prohash{$np}{check}==1); # only use 'good' coding gene models

		# AA351toAA358
		if ($pos =~ m/^AA(\d+)toAA(\d+)$/) {} else {die "wrong pos format: $pos!!!\ndifferent ID.txt format?\n";}
		$pos =~ m/^AA(\d+)toAA(\d+)$/;
		my ($b5,$b3)=($1*3-3,$2*3);
		if ($prohash{$np}{strand} eq '-')
		{
			my $t3=$prohash{$np}{cdslength}-$b5;
			my $t5=$prohash{$np}{cdslength}-$b3;
			$b5=$t5; $b3=$t3;
		}

		my (@pepb5,@pepb3);
		my $i=0; if (!defined($prohash{$np}{cmllength}[$i])) {die "$np,$pos,$b5,$b3,$i,$prohash{$np}{cmllength}[$i]\n";}
		if ($b5>$prohash{$np}{cmllength}[scalar(@{$prohash{$np}{cmllength}})-1]) {print "peptide b5 is downstream of last CDS exon!!!\n$np,$b5,$prohash{$np}{cmllength}[scalar(@{$prohash{$np}{cmllength}})-1]\n"; next;}
		while ( $b5>$prohash{$np}{cmllength}[$i] ) {$i++;}# print "$np,$pos,$b5,$b3,$i,$prohash{$np}{cmllength}[$i]\n";}

		my $j=0; $pepb5[$j]=$prohash{$np}{eb3}[$i]-($prohash{$np}{cmllength}[$i]-$b5);
		if ($b3>$prohash{$np}{cmllength}[scalar(@{$prohash{$np}{cmllength}})-1]) {print "peptide b3 is downstream of last CDS exon!!!\n$np,$b3,$prohash{$np}{cmllength}[scalar(@{$prohash{$np}{cmllength}})-1]\n"; next;}
		while ( $b3>$prohash{$np}{cmllength}[$i] ) 
		{
			$pepb3[$j]=$prohash{$np}{eb3}[$i];
			$i++;
			$j++;
			$pepb5[$j]=$prohash{$np}{eb5}[$i];
		}
		$pepb3[$j]=$prohash{$np}{eb3}[$i]-($prohash{$np}{cmllength}[$i]-$b3);

		# check whether the 1st 'exon' is empty
		if ( $pepb5[0] == $pepb3[0] )
		{
			#for ($i=0; $i<$#pepb3; $i++) { $pepb5[$i]=$pepb5[$i+1]; $pepb3[$i]=$pepb3[$i+1];  }
			#undef($pepb5[$#pepb3]);
			#undef($pepb3[$#pepb3]);
			shift(@pepb5);
			shift(@pepb3);
		}

		# print BED file
		print OUT "$prohash{$np}{chr}\t$pepb5[0]\t$pepb3[$#pepb3]\t$peptide\t900\t$prohash{$np}{strand}\t$pepb5[0]\t$pepb3[$#pepb3]\t0\t",$#pepb3+1,"\t";
		for ($i=0; $i<=$#pepb3; $i++) { print OUT  $pepb3[$i]-$pepb5[$i],","; } 
		print OUT "\t";
		for ($i=0; $i<=$#pepb3; $i++) { print OUT $pepb5[$i]-$pepb5[0],","; } 
		print OUT "\t$outfile";
		print OUT "\n";
		
	}
	#else {die "Not existed NP id: $np!!!\n";}
	#else {print "Not existed NP id: $np!!!\n";}
#=cut
}
close IN;
close OUT;

#------------------------------------------------------------------------------
sub checkFrame
{
	my ($hash)=@_;

	my @t=split(/\t/,$$hash{line});
	my $cdsStart=$t[5];
	my $cdsEnd=$t[6];
	my $exonStarts=$t[8];
	my $exonEnds=$t[9];

	if ($cdsStart==$cdsEnd) {return 0;} # non-coding genes

	my @b5=split(/,/,$exonStarts);
	my @b3=split(/,/,$exonEnds);

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
