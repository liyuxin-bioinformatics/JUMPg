#!/usr/bin/perl
use warnings;
use strict;


if (@ARGV != 4)
{
	die "Usage: extractFlankingGenomicSeq.pl <input> hg19_raw_genome RNAseq_Exon_Junction_novel.txt flank_nt <output> out\n";
}


open(GNM,"$ARGV[0]") or die "Can't open $ARGV[0]\n";
open(PS,"$ARGV[1]") or die "Can't open $ARGV[1]\n";
open(OUT,">$ARGV[3]");


#############################################################################
my $padding=$ARGV[2];
my $blockSize=750;
#my $lineNum_unit=15;	#number of lines that one single element of the array contains
#my $nt_line=50;				#nucleotide number per line
#############################################################################

my ($line, $i, $j, @est, $est_length, @seq, $est_name, $chrm, $mark, $line_length, $tmp, %genome);
my $seq='';
while ( $line=<GNM> )
{
	if ($line =~ m/>/)
	{
		if ($seq ne '') {
			my $maxblock=int(length($seq)/$blockSize);
			for (my $i=0; $i<$maxblock; $i++) {
				$genome{$chrm}[$i]=substr $seq, $i*$blockSize,$blockSize;
			}
			$genome{$chrm}[$maxblock]=substr $seq, $maxblock*$blockSize;
			$seq = '';
		}

		chomp($line);
		$line =~ s/>//;
		#$chrm=chr2num($line);
		$chrm=(split /[\t\s]/,$line)[0];

		#$i=0;	$j=0; 
		#$genome{$chrm}[$i]='';		
		next;
	} else {	
		chomp($line);
		if ($line =~ /\r$/) {
			 chop($line);
		}
		$seq.=$line;
	}
=head	
	#if ( $j>=$lineNum_unit ) { $i++; $genome[$chrm][$i]=''; $j=0; }
	if ( $j>=$lineNum_unit ) { $i++; $genome{$chrm}[$i]=''; $j=0; }
	#$line =~ m"(.*?)\r";
	$tmp=$line;
	#$genome[$chrm][$i].=$tmp;
	$genome{$chrm}[$i].=$tmp;
	$j++;
=cut
}
# last chr
		if ($seq ne '') {
			my $maxblock=int(length($seq)/$blockSize);
			for (my $i=0; $i<$maxblock; $i++) {
				$genome{$chrm}[$i]=substr $seq, $i*$blockSize,$blockSize;
			}
			$genome{$chrm}[$maxblock]=substr $seq, $maxblock*$blockSize;
			$seq = '';
		}

my ($plaID, $strand, $mode, $testStart, $testEnd);
$line=<PS>;
while ( $line=<PS> )
{
	next if ( $line =~ m/^#/ );
	chomp($line);
	#chr10   100147063       100150355       4       4       +       100147063       100150355       0,128,0
	#my ($chr,$e1,$e2,$id,$spt,$strand)=split(/\t/,$line);
	#ExonJunction    Gene    Annotation_Set  RefGene EnsGene AceView UCSC
	#10:100147064:+,10:100150355:+   PYROXD2 extended        None
	#1:883974:+,1:886507:+   2
	my @t=split(/\t/,$line);
	my ($junc,$strand)=split(/\t/,$line);

	#if ($t[0] =~ m/^(\d+)\:(\d+)\:\+\,\d+\:(\d+)\:\+$/) {} else {next;}
	if ($junc =~ m/^(.*?)\:(\d+)\:\+\,.*?\:(\d+)\:\+$/) {} else {next;}
	my ($chr,$e1,$e2)=($1,$2,$3);
	#$chrm="chr$chr";
	my $chrm=($chr =~ m/chr/)?$chr:"chr$chr";
	#$chrm=chr2num($chr);
	#next if ($chrm<=0);

	next unless (defined($genome{$chrm}));

	print OUT "\>$junc strand=$strand\n";
=head
	if (scalar(@t)==2) # simple format: 2 columns: junction & gene
	{
		print OUT "\>$t[0] Gene=$t[1]\n";
	}
	elsif (scalar(@t)>2) # ComBio format
	{
		#sampleN and maxR
		my ($sampleN,$maxR)=(0,0);
		for (my $i=7;$i<=$#t;$i++)
        	{
                	if ($t[$i]) { $sampleN++; }
	                if ( $t[$i]>$maxR ) { $maxR=$t[$i]; }
        	}

		print OUT "\>$t[0] Gene=$t[1]\|Annotation_Set=$t[2]\|sampleN=$sampleN\|maxRead=$maxR\n";
	}
	else { die "Wrong input format for $ARGV[1]!!!\n"; }
=cut
	genomicSeq($chrm,$e1-$padding+1,$e1);
	genomicSeq($chrm,$e2,$e2+$padding-1);
	print OUT "\n";
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
	#$chrm=7; 
	#my $pos_start=99529632;
	#my $pos_end=99529551;
	
	my ($i, $j, $tmp1, $tmp2, $l, $c);
	
	$i=int( $pos_start/($blockSize) ); # start block
	$j=int( $pos_end/($blockSize) ); # end block

	$tmp1=$pos_start-$i*$blockSize-1; # start block offset
	$tmp2=$pos_end-$j*$blockSize-1; # end block offset

	if ( $pos_start<$pos_end ) # positive strand
	{
		if ( $i==$j ) 
		{
			@seq=split('',$genome{$chrm}[$i]);	#	print OUT "@seq\n";
			for ( $l=$tmp1; $l<=$tmp2; $l++ ) { print OUT $seq[$l]; }
		}
		else
		{
			@seq=split('',$genome{$chrm}[$i]);
			for ( $l=$tmp1; $l<$blockSize; $l++ ) { print OUT $seq[$l]; }
			for ( my $k=$i+1; $k<=$j-1; $k++ )
		  {
  			@seq=split('',$genome{$chrm}[$k]);
	  		for ( $l=0; $l<$blockSize; $l++ ) { print OUT $seq[$l]; }
		  }
  		@seq=split('',$genome{$chrm}[$j]);
			for ( $l=0; $l<=$tmp2; $l++ ) { print OUT $seq[$l]; }
		}
	}
	else
	{
		if ( $i==$j ) 
	  {
		  @seq=split('',$genome{$chrm}[$i]);	#	print OUT "@seq\n";
	  	for ( $l=$tmp1; $l>=$tmp2; $l-- ) { $c=complement($seq[$l]);  print OUT $c; }
		  #print OUT "\n"; for ( $l=$tmp2; $l>=$tmp1; $l-- ) {   print OUT $c; }
  	}
	  else
  	{
  		@seq=split('',$genome{$chrm}[$i]);
		  for ( $l=$tmp1; $l>=0; $l-- ) { $c=complement($seq[$l]);  print OUT $c; }
		  for ( my $k=$i-1; $k>=$j+1; $k-- )
    	{
  	  	@seq=split('',$genome{$chrm}[$k]);
	  	  for ( $l=$blockSize-1; $l>=0; $l-- ) {  $c=complement($seq[$l]);  print OUT $c;  }
  	  }
    	@seq=split('',$genome{$chrm}[$j]);
		  for ( $l=$blockSize-1; $l>=$tmp2; $l-- ) {  $c=complement($seq[$l]);  print OUT $c;  }
  	}
	}
}


close(GNM);
close(PS);
close(OUT);
