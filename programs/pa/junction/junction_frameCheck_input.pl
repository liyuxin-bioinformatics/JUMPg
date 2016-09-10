#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=4)
{
	die "perl junction_frameCheck_input.pl junction_peptides_frame.txt novel_AA.fas AApadding out\n";
}

#----------------------------------------------------
my $AApadding=$ARGV[2];
#my $up_padding=67;
#my $down_padding=66;
my $up_padding=$AApadding*3;
my $down_padding=$AApadding*3;
#----------------------------------------------------

my %accpJun;
open(IN,"$ARGV[0]");
my $line=<IN>; chomp($line); 
while(<IN>)
{
	chomp; $line=$_;
	my @t=split(/\t/,$_);
	#outfile peptide XCorr   dCn     ppm     proteins        Junction        Gene    sampleN maxR    BS      BR      h3      h9      h27     juncContained   stopN
	#mmf37.119915.1.2.out    R.AAAIPPIVTK.V  3.2853  0.7924  0.119749060943182       12:124846840:+,12:124848225:+|R1;       12:124846840:+,12:124848225:+   NCOR2   5       212     135     212     50      205     23      1       0
	#mmf05.34029.1.3.out     K.RVDQELNGK.Q   2.3184  0.5411  1.51412543423758        15:59956319:+,15:59960302:+|R2; 15:59956319:+,15:59960302:+     BNIP2   5       25      14      7       8       25      2       1       2
	#my ($pros,$orig_pep,$junc,$pro,$juncContained)=($t[5],$t[1],$t[6],'',$t[15]);
	#my ($pros,$orig_pep,$junc,$pro,$juncContained)=($t[5],$t[1],$t[6],'',$t[$#t-1]);
	my $pro=$t[$#t];
	#next unless ($juncContained>0);

	# 2:190530177:+,2:190532490:+|F1;2:190526336:+,2:190530102:+|F2;
	#my @proteins=split(/\;/,$pros);
	#foreach my $p (@proteins) { $accpJun{$p}=''; }
	$accpJun{$pro}='';
}
close IN;

my (%prohash);
open(DB,"$ARGV[1]");
open(OUT,">$ARGV[3]");
print OUT "junction\tjunction_peptide\tjunction_strand\tjuncWithinPeptidePos1\tjuncWithinPeptidePos2\n";
while(<DB>)
{
	chomp;
	s/^>//;
	my ($pro,$ann)=split(/ /,$_);
	if (!defined($accpJun{$pro})) { my $seq=<DB>; next; }

	$prohash{$pro}{annotation}=$ann;
	#10:100143667:+,10:100144704:+|R0
	my ($junc,$frame)=split(/\|/,$pro);

	my $seq=<DB>;
	chomp($seq);
	# 67 + 1 + 1 + 66 = 135 nt
	# R0: 66/3 = 22 yu 0, 23 => juncPointU / juncPointD
	# R1: 65/3 = 21 yu 2, juncPointU=22, juncPointD=23
	# R2: 64/3 = 21 yu 1, 22 => juncPointU / juncPointD
	# F0: 67/3 = 22 yu 1, 23 => juncPointU / juncPointD
	# F1: 66/3 = 22 yu 0, 23 => juncPointU / juncPointD
	# F2: 65/3 = 21 yu 2, juncPointU=22, juncPointD=23
	$prohash{$pro}{seq}=$seq;
	$prohash{$pro}{frame}=$frame;
	$prohash{$pro}{stopN}=$seq =~ tr/*/*/;
	#if ($frame eq 'R0' || $frame eq 'F0' || $frame eq 'F1') { $prohash{$pro}{juncPointU}=$prohash{$pro}{juncPointD}=23; }
	#elsif ($frame eq 'R1' || $frame eq 'F2') { $prohash{$pro}{juncPointU}=22;$prohash{$pro}{juncPointD}=23; }
	#elsif ($frame eq 'R2') { $prohash{$pro}{juncPointU}=$prohash{$pro}{juncPointD}=22; }
	if ($frame eq 'R0' || $frame eq 'F0') {
		$prohash{$pro}{juncPointU}=$AApadding;
                $prohash{$pro}{juncPointD}=$AApadding+1;
        } else {
                $prohash{$pro}{juncPointU}=$prohash{$pro}{juncPointD}=$AApadding;
        }


	# genomic position and strand
	$frame =~ m/([RF])([012])/;
	my ($strand,$offset)=($1,$2);
	#10:100143667:+,10:100144704:+
	$junc =~ m/^(.*?)\:(\d+)\:\+\,.*?\:(\d+)\:\+$/; #print "$junc\n";
	my ($chr,$junc5,$junc3)=($1,$2,$3);
	$prohash{$pro}{chr}=$chr;
	$prohash{$pro}{junc_up}=$junc5;
	$prohash{$pro}{junc_dn}=$junc3-1; # zero-based
	if ($strand eq 'F') 
	{  
		$prohash{$pro}{strand}='+';
		$prohash{$pro}{b5}=$junc5-$up_padding+$offset-1; # zero-based
		#$prohash{$pro}{b3}=$junc3+$down_padding;
		$prohash{$pro}{b3}=$prohash{$pro}{b5}+length($prohash{$pro}{seq})*3              # b5 + protein seq length
					+( $prohash{$pro}{junc_dn} - $prohash{$pro}{junc_up} );  # junc length
	}
	else
	{
		$prohash{$pro}{strand}='-';
		#$prohash{$pro}{b5}=$junc5-$up_padding-1;  # zero-based
		$prohash{$pro}{b3}=$junc3+$down_padding-$offset;
		$prohash{$pro}{b5}=$prohash{$pro}{b3} - length($prohash{$pro}{seq})*3		# b3 - protein seq length
				- ( $prohash{$pro}{junc_dn} - $prohash{$pro}{junc_up} ); 	# junc length
	}

	# strand
	#my $strand;
	if ( $frame =~ m/^F/ ) { $strand='+'; } else { $strand='-'; }

	# split by stop codon
	my @fr=split(/\*/,$seq);
	my $accuml=0; my $i;
	for ($i=0; $i<=$#fr; $i++)
	{
		$accuml+=length($fr[$i]);
		last if ($accuml>=$prohash{$pro}{juncPointD});
	}
	next if ( 0<$i and $i<$#fr ); # Apply filter: junc-fragment in the middle
	$accuml=$accuml-length($fr[$i]);
	my $juncPosU=$prohash{$pro}{juncPointU}-$accuml;
	my $juncPosD=$prohash{$pro}{juncPointD}-$accuml;

	# select + - 10 AA
	my $spanPep = substr $fr[$i], max(0,$juncPosU-1-10), min($juncPosD-1+10,length($fr[$i])-1)-max(0,$juncPosU-1-10)+1;
	my $juncPos1=$juncPosU-max(0,$juncPosU-1-10);
	my $juncPos2=$juncPosD-max(0,$juncPosU-1-10);

	# print 
	print  OUT "$junc\t$spanPep\t$strand\t$juncPos1\t$juncPos2\n";
}
close DB;
close OUT;

#------------------------------
sub max
{
	my ($a,$b)=@_;
	if ($a>$b) { return $a; } else { return $b; }
}
sub min
{
	my ($a,$b)=@_;
	if ($a<$b) { return $a; } else { return $b; }
}
