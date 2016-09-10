#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/c

use strict;
use PrimarySeq;

my $prsq=PrimarySeq->new();


if (scalar(@ARGV)!=2)
{
        die "Usage: perl fas_tr6_noQC.pl in.fas minORFlength[AA]\n";
}

######################
my $minORFLength=$ARGV[1];
#####################

#my (%seq);
my $input=$ARGV[0];
my $seqhash=$prsq->parseFas($input);

foreach my $id (keys %$seqhash) {
	my @tr=$prsq->translate_6frames($seqhash->{$id});
	for (my $i=0; $i<3; $i++) {
		my @fr=split /\*/,$tr[$i];
		foreach my $f (@fr) {
			next if (length($f)<$minORFLength); # min ORF length
			my $aaStart=index $tr[$i],$f;
			my $aaEnd=$aaStart+length($f);
			print ">$id|F$i|",$aaStart*3+$i,'-',$aaEnd*3+$i,"\n";
			print "$f\n";
		}
	}
	for (my $i=3; $i<6; $i++) {
		my $tmp=$i-3;
		my @fr=split /\*/,$tr[$i];
		foreach my $f (@fr) {
			next if (length($f)<$minORFLength); # min ORF length
			my $aaStart=index $tr[$i],$f;
			my $aaEnd=$aaStart+length($f);
			print ">$id|R$tmp|",$aaStart*3+$tmp,'-',$aaEnd*3+$tmp,"\n";
			print "$f\n";
		}
	}
}
=head
sub ORFfragment {
	my $AAseq=@_;
	my @t=split /\*/,$AAseq;
	# %frg{$orf}=''
}

my @seq;
while(<IN>)
{




	$seq[0]=$_; chomp($seq[0]);
        $seq[1]=<IN>; chomp($seq[1]);
#        $seq[2]=<IN>; chomp($seq[2]);
#        $seq[3]=<IN>; chomp($seq[3]);

	my $id =$seq[0];
	$id =~ s/^>//;
#	$id =~ s/@//;
#	$id =~  m/^\>\w+\:(.*?)\#\w+\/(\d)$/;
#	$id="$1\/$2";
#	$id =~ s/\:/_/g; $id =~ s/\//-/;

	my @t=@seq;

	next if (length($t[1])<$minReadLength); #length control

	#for (my $i=0; $i<4; $i++) {print OUT "$t[$i]\n";}

	my @tr=$prsq->translate_6frames($t[1]);

	my $mark=0; # check if read through
	my $marks1=0; # check if exists a frame with a single stop codon
	for (my $i=0; $i<3; $i++)
	{
		if ($tr[$i] =~ /\*/)
		{
			my @tmp=($tr[$i] =~ /\*/);
			if (scalar(@tmp)==1) {$marks1=1;}
		}
		else
		{
			my $pptid="$id|F$i";
			print PPT ">$pptid\n",$tr[$i],"\n";
			print R2P "$id\t$pptid\n";
			$mark=1;
		}
	}
	for (my $i=3; $i<6; $i++)
        {
                if ($tr[$i] =~ /\*/)
                {
                        my @tmp=($tr[$i] =~ /\*/);
                        if (scalar(@tmp)==1) {$marks1=1;}
                }
                else 
                {
			my $tmp=$i-3;
                        my $pptid="$id|R$tmp";
                        print PPT ">$pptid\n",$tr[$i],"\n";
                        print R2P "$id\t$pptid\n";
			$mark=1;
                }
        }
	unless ($mark) {print LG "$seq[0]\n$seq[1]\n"; if ($marks1) {$singleStopCount++;}}
#	if ($marks1 and $mark) {$singleStopCount++;}
	#print OUT "$t[0]\n$t[1]\n";
	#foreach (@tr) {print OUT "$_\n";}
}
print LG "$singleStopCount\n";

close IN;
close PPT;
close R2P;
close LG;
=cut
