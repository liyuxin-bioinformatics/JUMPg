#!/usr/bin/perl
use warnings;
use strict;


if (@ARGV != 3)
{
	die "Usage: attachGenomicSeq.pl <input> hg19_raw_genome file.pep.gnSb.fpkm.pos <output> out\n";
}


open(GNM,"$ARGV[0]") or die "Can't open $ARGV[0]\n";
open(PS,"$ARGV[1]") or die "Can't open $ARGV[1]\n";
open(OUT,">$ARGV[2]");


#############################################################################
#my $padding=$ARGV[2];

my $lineNum_unit=15;	#number of lines that one single element of the array contains
my $nt_line=50;				#nucleotide number per line
#############################################################################

my ($line, $i, $j, @est, $est_length, @seq, $est_name, $chrm, $mark, $line_length, $tmp, @genome);
while ( $line=<GNM> )
{
	if ($line =~ m/>/)
	{
		chomp($line);
		$line =~ s/>//;
		$chrm=chr2num($line);

		$i=0;	$j=0; 
		if ( $chrm !~ /^\d+/ ) { die "Error: chrm=$chrm\n$line\n"; }
		$genome[$chrm][$i]='';		
		next;
	}
		
	if ( $j>=$lineNum_unit ) { $i++; $genome[$chrm][$i]=''; $j=0; }
	chomp($line);
	chop($line);
	#$line =~ m"(.*?)\r";
	$tmp=$line;
	$genome[$chrm][$i].=$tmp;
	$j++;
}

my ($plaID, $strand, $mode, $testStart, $testEnd);
$line=<PS>; chomp($line);
print OUT "$line\tgenomicSeq\n";
while ( $line=<PS> )
{
	next if ( $line =~ m/^#/ );
	chomp($line);
	#chr10   100147063       100150355       4       4       +       100147063       100150355       0,128,0
	#my ($chr,$e1,$e2,$id,$spt,$strand)=split(/\t/,$line);
	#ExonJunction    Gene    Annotation_Set  RefGene EnsGene AceView UCSC
	#10:100147064:+,10:100150355:+   PYROXD2 extended        None
	#outfile XCorr   dCn     peptide scanCount       readCount       readID  ucscID  mismatch        gap     readLength      alignLength     qStart  qEnd    levels  sStart  sEndsStrand  b5      b3      strand  frame   sType   geneSymbol      FPKM    genome_chr      genome_strand   genome_start    genome_end      exonCount       exonSizes       exonStarts
	#t027.201205.1.2.spout   79.50   0.230566037735849       GVDEVTIVDILTNR  2       2       HWI-ST1200:136:C2R9HACXX:7:2313:5820:5367/1     uc010bgj.3      1       0       101 101      1       101     1       101     201     -       147     189     +       0       CDS     ANXA2   471.452501      chr15   -       60656681        606745412       41,1,0,17859,
	my @t=split(/\t/,$line);
	#my ($chr,$start,$strand,$exonCount,$exonSizes,$exonStarts)=($t[$#t-6],$t[$#t-4],$t[$#t-5],$t[$#t-2],$t[$#t-1],$t[$#t]);
	my ($chr,$start,$strand,$exonCount)=($t[$#t-6],$t[$#t-4],$t[$#t-5],$t[$#t-2]);
	my @exonSizes=split /\,/,$t[$#t-1];
	my @exonStarts=split /\,/,$t[$#t];

	# chr
	#$chr="chr$chr";
	next if ($chr =~ m/6_ssto_hap7/);
	$chrm=chr2num($chr);
	next if ($chrm<=0);

	# exon boundaries: 0-based
	my (@b5,@b3);
	for (my $i=0; $i<$exonCount;$i++)
	{
		$b5[$i]=$start+$exonStarts[$i];
		$b3[$i]=$b5[$i]+$exonSizes[$i];
	}

	# extract seq
	my $seq='';
	for (my $i=0; $i<$exonCount;$i++)
	{
		my $testStart=$b5[$i]+1;
		my $testEnd=$b3[$i];

		$seq.=genomicSeq($chrm,$testStart,$testEnd);
	}
	print OUT "$line\t$seq\n";

=head	
	if ( $strand eq '+' )	
	{
		genomicSeq($chrm,$testStart,$testEnd);
	}
	else
	{
		genomicSeq($chrm,$testEnd,$testStart);
	}
=cut	
}


#############################################
sub chr2num
{
	my ($chrm) = @_;

	unless ($chrm =~ m"chr")	{ die "Chr unmatchable!!!\n$chrm"; }
	$chrm =~ s/chr//;
	
	if ( $chrm =~ m/random/ )
	{
		#18_gl000207_random
		unless ($chrm =~ m"(\d+)_gl000(\d+)_random") { die "Chr unmatchable!!![line 45]\n$chrm"; }
		$chrm=$2;
	  if ( $chrm !~ /^\d+/ ) { die "Error: 1st if chrm=$chrm\n"; }
	}
	elsif ( $chrm =~ m/Un/ )
	{
		#Un_gl000211
		unless ( $chrm =~ m"Un_gl000(\d+)" ) { die "Chr unmatchable!!![line 52]\n$chrm"; }
		$chrm=$1;			
	}
	else
	{		
		  if ( $chrm =~ m"X" ) { $chrm=23; }
		  elsif ( $chrm =~ m"Y" ) { $chrm=24; }
		  elsif ( $chrm =~ m"M" ) { $chrm=25; }
		  #else { $chrm=-1; }
		  #if ( $chrm !~ /^\d+/ ) { die "Error: 2nd if chrm=$chrm\n"; }
		  if ( $chrm !~ /^\d+/ ) { $chrm=-1;  }
	}

	return $chrm;
	
}

sub complement
{
	my ($b) = @_;
	#my $c;
	if ( $b eq 'A' ) { return ('T'); }
	if ( $b eq 'T' ) { return ('A'); }
	if ( $b eq 'C' ) { return ('G'); }
	if ( $b eq 'G' ) { return ('C'); }
	
	if ( $b eq 'a' ) { return ('t'); }
	if ( $b eq 't' ) { return ('a'); }
	if ( $b eq 'c' ) { return ('g'); }
	if ( $b eq 'g' ) { return ('c'); }
}

#############################################
sub genomicSeq
{
	my ($chrm, $pos_start, $pos_end)=@_;
	my $sequence='';
	#$chrm=7; 
	#my $pos_start=99529632;
	#my $pos_end=99529551;
	
	my ($i, $j, $tmp1, $tmp2, $l, $c);
	
	$i=int( $pos_start/($nt_line*$lineNum_unit) );
	$j=int( $pos_end/($nt_line*$lineNum_unit) );

	$tmp1=$pos_start-$i*$nt_line*$lineNum_unit-1;
	$tmp2=$pos_end-$j*$nt_line*$lineNum_unit-1;

	if ( $pos_start<$pos_end )
	{
		if ( $i==$j ) 
		{
			@seq=split('',$genome[$chrm][$i]);	#	print OUT "@seq\n";
			#for ( $l=$tmp1; $l<=$tmp2; $l++ ) { print OUT $seq[$l]; }
			for ( $l=$tmp1; $l<=$tmp2; $l++ ) { $sequence.=$seq[$l]; }
		}
		else
		{
			@seq=split('',$genome[$chrm][$i]);
			#for ( $l=$tmp1; $l<$nt_line*$lineNum_unit; $l++ ) { print OUT $seq[$l]; }
			for ( $l=$tmp1; $l<$nt_line*$lineNum_unit; $l++ ) { $sequence.=$seq[$l]; }
			for ( my $k=$i+1; $k<=$j-1; $k++ )
		  {
  			@seq=split('',$genome[$chrm][$k]);
	  		#for ( $l=0; $l<$nt_line*$lineNum_unit; $l++ ) { print OUT $seq[$l]; }
	  		for ( $l=0; $l<$nt_line*$lineNum_unit; $l++ ) { $sequence.=$seq[$l]; }
		  }
  		@seq=split('',$genome[$chrm][$j]);
			#for ( $l=0; $l<=$tmp2; $l++ ) { print OUT $seq[$l]; }
			for ( $l=0; $l<=$tmp2; $l++ ) { $sequence.=$seq[$l]; }
		}
	}
	else
	{
		if ( $i==$j ) 
	  {
		  @seq=split('',$genome[$chrm][$i]);	#	print OUT "@seq\n";
	  	#for ( $l=$tmp1; $l>=$tmp2; $l-- ) { $c=complement($seq[$l]);  print OUT $c; }
	  	for ( $l=$tmp1; $l>=$tmp2; $l-- ) { $c=complement($seq[$l]);  $sequence.=$c; }
  	}
	  else
  	{
  		@seq=split('',$genome[$chrm][$i]);
		  #for ( $l=$tmp1; $l>=0; $l-- ) { $c=complement($seq[$l]);  print OUT $c; }
		  for ( $l=$tmp1; $l>=0; $l-- ) { $c=complement($seq[$l]);  $sequence.=$c; }
		  for ( my $k=$i-1; $k>=$j+1; $k-- )
    	{
  	  	@seq=split('',$genome[$chrm][$k]);
	  	  #for ( $l=$nt_line*$lineNum_unit-1; $l>=0; $l-- ) {  $c=complement($seq[$l]);  print OUT $c;  }
	  	  for ( $l=$nt_line*$lineNum_unit-1; $l>=0; $l-- ) {  $c=complement($seq[$l]);  $sequence.=$c;  }
  	  }
    	@seq=split('',$genome[$chrm][$j]);
		  #for ( $l=$nt_line*$lineNum_unit-1; $l>=$tmp2; $l-- ) {  $c=complement($seq[$l]);  print OUT $c;  }
		  for ( $l=$nt_line*$lineNum_unit-1; $l>=$tmp2; $l-- ) {  $c=complement($seq[$l]);  $sequence.=$c;  }
  	}
	}
	return $sequence;
}


close(GNM);
close(PS);
close(OUT);
