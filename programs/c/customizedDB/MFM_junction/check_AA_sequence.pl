#!/usr/bin/perl -w -I /usr/local/lib/perl5

use strict;

if (scalar(@ARGV)!=2)
{
        die "Usage: perl check_AA_sequence.pl translated_peptides.fas out\n";
}

######################
#$ARGV[2] ||= 'KR';
#$ARGV[3] ||= 'P';

my $minFragmentLength=7;
my $minOverhangAA=3;
#my @sites = split(/\s*/,$ARGV[2]);
#my $uncleavage = $ARGV[3];
######################

open(IN,"$ARGV[0]") or  die "Cannot open $ARGV[0]\n";
open(OUT,">$ARGV[1]");

while(<IN>)
{
	chomp;
	my $id=$_;
	my $seq=<IN>;
        chomp($seq);
	my $oriSeq=$seq;

	# in silico digest
        $seq =~ s/\*/\* /g;
        $seq =~ s/ \Z//;
        my @a = split(/ /,$seq);

        #next if (scalar(@a)<=1); # no tryptic sites

	#count length to anchor the central peptide
	my $cml=0;
	for my $i (0..$#a)
	{
		$cml+=length($a[$i]);
		#if ($cml>length($oriSeq)/2)
		if ($cml>length($oriSeq)/2 + $minOverhangAA)
		{
			if (length($a[$i])>=$minFragmentLength # min peptide length
				&& $i==0)	# ensure no stop codon before ss
			{
				print OUT "$id\n";
				#print OUT "$id $cml ",length($oriSeq)/2,"\n";
				$a[$i] =~ s/\*$//;
				print OUT "$a[$i]\n";
			}
			last;
		}
	}
}

close IN;
close OUT;

