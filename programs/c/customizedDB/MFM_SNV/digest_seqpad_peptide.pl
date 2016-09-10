#!/usr/bin/perl -w -I /usr/local/lib/perl5

use strict;

if (scalar(@ARGV)<3)
{
        die "Usage: perl digest_seqpad_peptide.pl seqpad.fas seqpad_length[int, e.g., 30] out enzymeSite[e.g. KR] uncleavage[e.g. P]\n";
}

######################
my $seqpad_length=$ARGV[1];

$ARGV[3] ||= 'KR';
$ARGV[4] ||= 'P';

my $minFragmentLength=7;
my @sites = split(/\s*/,$ARGV[3]);
my $uncleavage = $ARGV[4];
######################

open(IN,"$ARGV[0]") or  die "Cannot open $ARGV[0]\n";
open(OUT,">$ARGV[2]");

while(<IN>)
{
	chomp;
	my $id=$_;
	my $seq=<IN>;
        chomp($seq);
	my $oriSeq=$seq;

	# in silico digest
        for my $i (0..$#sites)
        {
                $seq =~ s/($sites[$i])(?!$uncleavage)/$1 /g;
        }
        $seq =~ s/ \Z//;
        my @a = split(/ /,$seq);

        next if (scalar(@a)<=1); # no tryptic sites

	# check if non-sense mutation
	if (length($oriSeq)<2*$seqpad_length+1 && $id =~ /\*/)
	{
		$a[$#a] =~ s/\*//;
		if ( length($a[$#a])>=$minFragmentLength ) { print OUT "$id\n$a[$#a]\n"; }
		next;
	}

	#
	my $fl=substr($oriSeq,0,1); my $ll=substr($oriSeq,length($oriSeq)-1,1);
	if (length($oriSeq)==2*$seqpad_length+1) { print_posSpecific_pep(\@a,$id,$seqpad_length+1); next;}
	elsif ( $ll eq '*' ) { $a[$#a] =~ s/\*//; print_posSpecific_pep(\@a,$id,$seqpad_length+1); next; }
	#elsif ( $fl eq 'M' ) { print_posSpecific_pep(\@a,$id,length($oriSeq)-$seqpad_length); next;}
	else { print_posSpecific_pep(\@a,$id,length($oriSeq)-$seqpad_length); next;}
	#else {die "Special case? $id\n";}
}

close IN;
close OUT;

#----------------------------------------------
sub print_posSpecific_pep
{
	my ($a,$id,$pos)=@_;

	#count length to anchor the central peptide
	my $cml=0;
	for my $i (0..(scalar(@$a)-1))
	{
		$cml+=length($$a[$i]);
		if ($cml>=$pos) # check if central peptide
		{
			if (length($$a[$i])>=$minFragmentLength) # additional conditions
			{
				print OUT "$id\n";
				print OUT "$$a[$i]\n";
			}
			last;
		}
	}
}
