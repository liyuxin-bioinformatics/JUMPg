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
open(IN,'/home/yli4/annotations/refSeq_mRNA2protein_072514.txt') or die "Cannot open annotatiion file!!!\n";
while(<IN>)
{
	next if (/^#/);
	chomp;
	my @t=split(/\t/,$_);
	my ($nm,$np)=($t[1],$t[18]);

	if (!defined($prohash{$np}))
	{
	$prohash{$np}{nm}=$nm;
	$prohash{$np}{strand}=$t[3];
	$prohash{$np}{chr}=$t[2];
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

# parse ID.txt and map to genomic positions
open(IN,"$ARGV[0]");
open(OUT,">$ARGV[2]");
print OUT 'track name=',$trackName,' description="',$trackName,'" visibility=3 itemRgb="On"',"\n";
while(<IN>)
{
        if (/^Database/) {  next; }
        if (/^Peptide/) {  next; }

        chomp;
        my @t=split(/\;/,$_);

	my ($np,$pos,$peptide,$outfile)=($t[1],$t[$#t-1],$t[0],$t[2]);
	next if ($np =~ /\+/ || /CON_/);
	$np =~ s/\.(\d+)$//;
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
		for ($i=0; $i<=$#pepb3; $i++) { print OUT $pepb3[$i]-$pepb5[$i],","; } 
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
	if ($t[6]==$t[7]) {return 0;} # non-coding genes

	my @b5=split(/,/,$t[9]);
	my @b3=split(/,/,$t[10]);

	my @exonl;
	for (my $i=0;$i<=$#b5; $i++) { $exonl[$i]=$b3[$i]-$b5[$i]; }
	
	my $i=0;
	while ($t[6]>$b3[$i]) {$i++;} my $j=0; $$hash{eb5}[$j]=$t[6];
	while ($t[7]>$b3[$i])
	{
		$$hash{eb3}[$j]=$b3[$i];
		$i++;
		$j++;
		$$hash{eb5}[$j]=$b5[$i];
	}
	$$hash{eb3}[$j]=$t[7];

	for ($i=0;$i<scalar(@{$$hash{eb3}}); $i++) { $$hash{elength}[$i]=$$hash{eb3}[$i]-$$hash{eb5}[$i]; }
	my $cdsL=0;
	for ($i=0;$i<scalar(@{$$hash{eb3}}); $i++) { $cdsL+=$$hash{elength}[$i]; $$hash{cmllength}[$i]=$cdsL; }
	$$hash{cdslength}=$cdsL;
	if ( $cdsL % 3 == 0 )  { return 1;}
	else {return -1;}
}
