#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=5)
{
	die "perl ID2bed.pl valideted_ID.txt novel_AA.fas AApadding trackName out\n";
}

#----------------------------------------------------
my $trackName=$ARGV[3];
my $AApadding=$ARGV[2];
my $up_padding=$AApadding*3;
my $down_padding=$AApadding*3;
#----------------------------------------------------

my (%prohash);
open(DB,"$ARGV[1]");
while(<DB>)
{
	chomp;
	s/^>//;
	my ($pro,$ann)=split(/ /,$_);
	$prohash{$pro}{annotation}=$ann;
	#10:100143667:+,10:100144704:+|R0
	my ($junc,$frame)=split(/\|/,$pro);

	my $seq=<DB>;
	chomp($seq);
	$prohash{$pro}{seq}=$seq;
	$prohash{$pro}{frame}=$frame;
	$prohash{$pro}{stopN}=$seq =~ tr/*/*/;

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
		#$prohash{$pro}{b5}=$junc5-$up_padding+$offset-1; # zero-based
		#$prohash{$pro}{b3}=$prohash{$pro}{b5}+length($prohash{$pro}{seq})*3              # b5 + protein seq length
		#			+( $prohash{$pro}{junc_dn} - $prohash{$pro}{junc_up} );  # junc length
		$prohash{$pro}{b5}=$junc5-$up_padding+$offset; # zero-based
		$prohash{$pro}{b3}=$prohash{$pro}{b5}+length($prohash{$pro}{seq})*3              # b5 + protein seq length
					+( $prohash{$pro}{junc_dn} - $prohash{$pro}{junc_up} );  # junc length
	}
	else
	{
		$prohash{$pro}{strand}='-';
		#$prohash{$pro}{b3}=$junc3+$down_padding-$offset;
		#$prohash{$pro}{b5}=$prohash{$pro}{b3} - length($prohash{$pro}{seq})*3		# b3 - protein seq length
		#		- ( $prohash{$pro}{junc_dn} - $prohash{$pro}{junc_up} ); 	# junc length
		$prohash{$pro}{b3}=$junc3+$down_padding-$offset - 1;
		$prohash{$pro}{b5}=$prohash{$pro}{b3} - length($prohash{$pro}{seq})*3 		# b3 - protein seq length
				- ( $prohash{$pro}{junc_dn} - $prohash{$pro}{junc_up} ); 	# junc length
	}
}
close DB;

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[4]");
my $line;#=<IN>; chomp($line); 
#print OUT 'track name=peptides description="MM non-refSeq junction peptides" visibility=3 itemRgb="On"',"\n";
print OUT 'track name=',$trackName,' description="',$trackName,'" visibility=3 itemRgb="On"  color="255,0,247"',"\n";
while(<IN>)
{
	next unless (/;/);
	chomp; $line=$_;
	my @t=split(/\;/,$_);
	# K.LREPIDNLTPEER.D;20:34313077:+,20:34317384:+|R1;/home/yli4/sudemycin_AS/JUMPg_2015Nov/gnm_test1/intermediate/Pottegrp_090215_TMT_F22/Pottegrp_090215_TMT_F22.1/Pottegrp_090215_TMT_F22.14159.1.3.spout;1810.99076;1810.9859818736;-0.232829059611619;50.53;0.73975855927172;0/0;0;137799;1;1;2;AA29toAA41;56.02%
	my ($orig_pep,$pro,$outfile)=($t[0],$t[1],$t[2]);

	#outfile peptide XCorr   dCn     ppm     proteins        Junction        Gene    sampleN maxR    BS      BR      h3      h9      h27
	#mmf37.119915.1.2.out    R.AAAIPPIVTK.V  3.2853  0.7924  0.119749060943182       12:124846840:+,12:124848225:+|R1;       12:124846840:+,12:124848225:+   NCOR2   5       212     135     212     50      205     23
	#my ($pros,$orig_pep,$junc,$pro)=($t[5],$t[1],$t[6],'');

	# rm mod and side AAs
	my $pep=$orig_pep; $pep =~ s/[\@\#\*\^\~\$\.]//g; chop($pep); $pep=reverse($pep); chop($pep); $pep=reverse($pep);

	# AA pos in protein
	my ($proseq,$juncP);
	if (defined($prohash{$pro})) { $proseq=$prohash{$pro}{seq}; } else { print "Not defnined protein: $pro!!!\n"; next;}
	my $fpos = index($proseq, $pep) + 1;
	my $lpos = $fpos + length($pep) - 1;

	# check juncPoint
        if ( $fpos<=$prohash{$pro}{juncPointU} and $prohash{$pro}{juncPointU}<=$lpos
        and $fpos<=$prohash{$pro}{juncPointD} and $prohash{$pro}{juncPointD}<=$lpos) { $juncP=1; }
        else { $juncP=0; }

	# deal with - strand
	my ($juncPointD_,$juncPointU_);
	if ($prohash{$pro}{strand} eq '-')
	{
		my $lpos_ = length($prohash{$pro}{seq})-$fpos+1;
		my $fpos_ = length($prohash{$pro}{seq})-$lpos+1;
		$fpos = $fpos_;
		$lpos = $lpos_;

		$juncPointD_ = length($prohash{$pro}{seq})-$prohash{$pro}{juncPointU}+1;
		$juncPointU_ = length($prohash{$pro}{seq})-$prohash{$pro}{juncPointD}+1;
		#$prohash{$pro}{juncPointD}=$juncPointD_;
		#$prohash{$pro}{juncPointU}=$juncPointU_;
	}
	else
	{
		$juncPointD_ = $prohash{$pro}{juncPointD};
		$juncPointU_ = $prohash{$pro}{juncPointU};
	}

	#if ( $orig_pep eq 'R.QQQNSNIFFLADR.T' ) {print  "$prohash{$pro}{b5}\t$prohash{$pro}{junc_up}\t$prohash{$pro}{junc_dn}\t$prohash{$pro}{b3}\n";}

	#my $green=0;
	my $color=0;
	# print BED
	if ($juncP)
	{
		# AA to genomic position
		my ($b5,$b3)=(0,0);
		$b5=$prohash{$pro}{b5}+($fpos-1)*3;
		$b3=$prohash{$pro}{b3}-(length($prohash{$pro}{seq})-$lpos)*3;


		if ($prohash{$pro}{junc_up}==$b5)	# 1st exon is pseudo
		{
			$juncP=0;
			$b5=$prohash{$pro}{junc_dn};
			print OUT "chr$prohash{$pro}{chr}\t$b5\t$b3\t$orig_pep\t900\t$prohash{$pro}{strand}\t$b5\t$b3\t$color\t";
			print OUT "1\t",$b3-$b5,"\t0,\t$outfile\n";
		}
		elsif ($b3==$prohash{$pro}{junc_dn}) # 2nd exon is pseudo
		{
			$juncP=0;
			$b3=$prohash{$pro}{junc_up};
			print OUT "chr$prohash{$pro}{chr}\t$b5\t$b3\t$orig_pep\t900\t$prohash{$pro}{strand}\t$b5\t$b3\t$color\t";
			print OUT "1\t",$b3-$b5,"\t0,\t$outfile\n";
		}
		else  # real junc. pep.
		{
			print OUT "chr$prohash{$pro}{chr}\t$b5\t$b3\t$orig_pep\t900\t$prohash{$pro}{strand}\t$b5\t$b3\t$color\t";
			print OUT "2\t",$prohash{$pro}{junc_up}-$b5,',',$b3-$prohash{$pro}{junc_dn},"\t0,",$prohash{$pro}{junc_dn}-$b5,"\t$outfile\n";
		}
	}
	else
	{
		#print "$fpos  $lpos  $prohash{$pro}{juncPointU}  $prohash{$pro}{juncPointD}\n";
		# AA to genomic position
		my ($b5,$b3)=(0,0);
		#if ( $lpos<= $prohash{$pro}{juncPointU})
		if ( $lpos<= $juncPointU_)
		{
			$b5=$prohash{$pro}{b5}+($fpos-1)*3;
			$b3=$prohash{$pro}{b5}+$lpos*3;
		}
		else
		{
			$b5=$prohash{$pro}{b3}-(length($prohash{$pro}{seq})-$fpos+1)*3;
			$b3=$prohash{$pro}{b3}-(length($prohash{$pro}{seq})-$lpos)*3;
		}

		print OUT "chr$prohash{$pro}{chr}\t$b5\t$b3\t$orig_pep\t900\t$prohash{$pro}{strand}\t$b5\t$b3\t$color\t";
		print OUT "1\t",$b3-$b5,"\t0,\t$outfile\n";
	}
}
close IN;
close OUT;
