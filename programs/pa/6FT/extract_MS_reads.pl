#!/usr/bin/perl
use strict;
use warnings;
use Storable;

# purposes:
# 1) extract reads that give rise to RDT peptides
# 2) peptide to read position (preparing for visualization): peptide, readID, pos
# 3) print out2read.hash

if (scalar(@ARGV)!=4)
{
	die "usage: perl extract_MS_reads.pl ID.txt seqIDs reads.fas out_prefix\n";
}

# parse ID.txt, build %rdtpep{$intpep}{$readID}=$pos
#print "reading ID.txt\n";
my (%tmppep,%rdtpep,%reads,%out2read);
open(IN,"$ARGV[0]") or die "cannot open $ARGV[0]\n";
print "reading $ARGV[0]\n";
while(<IN>)
{
	next if (/^Database/);
	next if (/^Peptide/);

	my @t=split /\;/,$_;
	my $pep=$t[0];

	chop($pep);chop($pep);$pep=reverse($pep);
	chop($pep);chop($pep);$pep=reverse($pep);
	my $nomod=$pep;
	$nomod =~ s/[\@\*\%]//g;

	# build %out2read
	my @s=split /\//,$t[2];
	my $out=$s[$#s];
	if (!defined($out2read{$out}))
	{
		#$tmppep{$pep}{$out}='';
		$tmppep{$nomod}{$out}='';

		$out2read{$out}{XCorr}=$t[6];
		$out2read{$out}{dCn}=$t[7];
		$out2read{$out}{peptide}=$pep;
		$out2read{$out}{scanCount}=1;
	}
	else
	{
	}
}
close IN;

open(IN,"$ARGV[1]") or die "cannot open $ARGV[1]\n";
print "reading $ARGV[1]\n";
while(<IN>)
{
	chomp;
	my ($id,$pep,$ids)=split /\t/,$_;

	chop($pep);$pep=reverse($pep);
	chop($pep);$pep=reverse($pep);

	if (defined($tmppep{$pep}))
	{
		my @readids=split /\,/,$ids;
		foreach my $r (@readids)
		{
			next if ($r eq '');
			my ($readID,$frame,$pos)=split /\|/,$r;
			$rdtpep{$pep}{$readID}{pos}=$pos;
			$rdtpep{$pep}{$readID}{frame}=$frame;
			$reads{$readID}='';

			foreach my $out (keys %{$tmppep{$pep}})
			{
				$out2read{$out}{readID}{$readID}='';
			}
		}
	}
}
close IN;
#store \%out2read ,'out2read.hash';
store \%out2read ,"$ARGV[3]_out2read.hash";

# print _pep2readPos.txt
open(OUT,">$ARGV[3]_pep2readPos.txt");
print OUT "peptide\treadID\tpos\tframe\n";
foreach my $pep (keys %rdtpep)
{
	foreach my $readID (keys %{$rdtpep{$pep}})
	{
		print OUT "$pep\t$readID\t$rdtpep{$pep}{$readID}{pos}\t$rdtpep{$pep}{$readID}{frame}\n";
	}
}
close OUT;

#=head
open(IN,"$ARGV[2]");
open(OUT,">$ARGV[3]_reads.fas");

while(<IN>)
{
	chomp;
	s/>//;
	my $id=$_;

	my $seq=<IN>;
	chomp($seq);

	if (defined($reads{$id}))
	{
		print OUT ">$id\n$seq\n";
	}
}

close IN;
close OUT;
#=cut
