#!/usr/bin/perl -w -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/g

use strict;
use PrimarySeq;

my $prsq=PrimarySeq->new();


if (scalar(@ARGV)!=2)
{
        die "Usage: perl junction_seq_translation.pl RNAseq_Exon_Junction_novel.fas  out\n";
}

my $output=$ARGV[1];

open(FAS,"$ARGV[0]") or die "cannot open $ARGV[0].\n";
#open(STR,"$ARGV[1]") or die "cannot open $ARGV[1].\n";
open(OUT,">$output"); close OUT;



=head
my (%strand);
while(<STR>)
{
        chomp;
        next if (/^#/);
        my @tmp=split(/\t/,$_);
        #uc001aaa.3      chr1    +       11873   14409   DDX11L1
        $strand{$tmp[5]}=$tmp[2];
}
my (%junc);
while(<STR>)
{
	chomp;
	#JUNC222989      chrX    153043549       153043851       4       PLXNB3  +
	my @tmp=split(/\t/,$_);
	$junc{$tmp[0]}{'strand'}=$tmp[6];
	$junc{$tmp[0]}{'gene'}=$tmp[5];
}
=cut
my $title;
while($title=<FAS>)
{
	chomp($title);
	#>10:100147064:+,10:100150355:+ Gene=PYROXD2|Annotation_Set=extended|sampleN=5|maxRead=5
	#>1:928737:+,1:986633:+ strand=1
	my ($id,$annotation)=split(/\s/,$title);
	$id =~ s/>//;

	if ($annotation =~ m/^strand\=(.*?)$/) {}
	#elsif ($annotation =~ m/^Gene\=(.*?)$/) {}
	my $strand=$1;
	my $st='';
	if ($strand == 1) {
		$st='+';
	} elsif ($strand == 2) {
		$st='-';
	}

	#if (defined($strand{$gene})) { $st=$strand{$gene}; }

	my $seq=<FAS>; 
	chomp($seq);
	my @tr;
	#my $seqID="$id\|$junc{$id}{gene}";
	if ( $st eq '+' ) { @tr=$prsq->translate_3frames($seq); printSeq(\@tr,$id,'F',$output,$annotation); }
	elsif ( $st eq '-' ) { $seq=$prsq->revcom($seq); @tr=$prsq->translate_3frames($seq); printSeq(\@tr,$id,'R',$output,$annotation);}
	else { @tr=$prsq->translate_6frames($seq); printSeq(\@tr,$id,'B',$output,$annotation);}

}

close FAS;
#close STR;

#---------------------------------------------------------------------
sub printSeq
{
	my ($tr,$id,$s,$out,$annotation)=@_;

	open(OUT,">>$out");

	if ($s eq 'F' or $s eq 'R')
	{
		for (my $i=0; $i<scalar(@$tr); $i++)
		{
			print OUT ">$id\|$s$i $annotation\n";
			print OUT "$$tr[$i]\n";
		}
	}
	elsif ($s eq 'B')
	{
		for (my $i=0; $i<3; $i++)
		{
			print OUT ">$id\|F$i $annotation\n";
			print OUT "$$tr[$i]\n";
		}
		for (my $i=3; $i<6; $i++)
		{
			print OUT ">$id\|R",$i-3," $annotation\n";
			print OUT "$$tr[$i]\n";
		}
	}

	close OUT;
}
