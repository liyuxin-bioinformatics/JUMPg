#!/usr/bin/perl -w 

use strict;


if (scalar(@ARGV)!=4)
{
        die "Usage: perl tryptic_digest.pl in.fas enzymeSite[e.g. KR] uncleavage[e.g. P] out\n";
}



######################
my $minFragmentLength=7;
my @sites = split(/\s*/,$ARGV[1]);
my $uncleavage = $ARGV[2];
my $input=$ARGV[0];
my $output=$ARGV[3];
#####################

open(IN,"$input") or die "cannot open $input\n";
open(OUT,">$output");

#my ($idLine);
while(<IN>)
{
	my $idLine=$_; chomp($idLine);
	my $seq=<IN>; chomp($seq);

	# in silico digest
 	for my $i (0..$#sites)
	{
                $seq =~ s/($sites[$i])(?!$uncleavage)/$1 /g;
        }
        $seq =~ s/ \Z//;
        my @a = split(/ /,$seq);

	next if (scalar(@a)<=1); # no tryptic sites

	# parse ID
	my $idInfor=parseID($idLine);

	# calculate relative position in RNAseq read
	my %ps;
	for my $i (0..$#a) { $ps{$i}{length}=length($a[$i]); }
	$ps{0}{end}=$ps{0}{length};
	for my $i (1..$#a) 
	{
		$ps{$i}{start}=$ps{$i-1}{end}+1; 
		$ps{$i}{end}=$ps{$i-1}{end}+$ps{$i}{length};
	}
	for my $i (1..$#a)
	{
		$ps{$i}{start}=($ps{$i}{start}-1)*3+$idInfor->{offset};
		$ps{$i}{end}=($ps{$i}{end})*3+$idInfor->{offset};
	}
	
	# output
	for my $i (1..($#a-1))
	{
		if ( length($a[$i])>=$minFragmentLength and $a[$i] !~ /\*/ )
		{
			print OUT '>',$idInfor->{readID},'|',$idInfor->{strand},'|',$ps{$i}{start},'-',$ps{$i}{end},"\n";
			my @b=split(//,$a[$i-1]); my @c=split(//,$a[$i+1]);
			print OUT "$a[$i]\n";
		}
	}

	#for last fragment
	if (length($a[$#a])>=$minFragmentLength && $a[$#a] =~ m/\*$/)
	{
		$a[$#a] =~ s/\*$//;
#		my @c=split(//,$a[$#a]); my $mark=0;
#		for my $i (0..$#sites) #loop to see if this final fragment is fully trptic
#		{
#			if ($c[$#c] eq $sites[$i]) { $mark=1; }
#		}

#		if ($mark)
#		{
			my $i=$#a;
			print OUT '>stop|',$idInfor->{readID},'|',$idInfor->{strand},'|',$ps{$i}{start},'-',$ps{$i}{end},"\n";
			my @b=split(//,$a[$i-1]);
			print OUT "$a[$i]\n";
#		}
	}
}



close IN;
close OUT;

#---------------------------------------------------------------------------------------------
sub parseID
{
	my ($id)=@_;

	$id =~ s/^\>//;
	my @a=split(/\|/,$id);

	my $idInfor;
	$idInfor->{readID}=$a[0];
#	my @b=split(//,$a[1]);
	$idInfor->{strand}=$a[1];
	$idInfor->{offset}=0;

	return $idInfor;
}
