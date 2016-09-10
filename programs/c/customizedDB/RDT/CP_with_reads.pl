#!/usr/bin/perl -w
use strict;

if (scalar(@ARGV)!=3)
{
	die "Usage: perl CP_with_reads.pl gene_model.fas (ref. used in mapping, for gene length esimation) file.map out\n";
}

my %cp=establish_CP_hash($ARGV[0]);

parse_map_file($ARGV[1],\%cp);

open(OUT,">$ARGV[2]");
foreach my $g (keys %cp) 
{
	print OUT "$g\t",$cp{$g}{length},"\t",$cp{$g}{readSpt},"\t",$cp{$g}{mut},"\n";
}
close OUT;

#classify_CP($ARGV[2],\%cp);
#-----------------------------
sub establish_CP_hash
{
	my ($in)=@_;
	open(IN,"$in") or die "Cannot open $in\n";

	my (%cp,$id);
	while(<IN>)
	{
		chomp;
		if (/^>/)
		{
			my @tmp=split(/\s/,$_);
			$tmp[0] =~ s/>hg19_knownGene_//;
			$id=$tmp[0];
			$cp{$id}{'length'}=0;
			$cp{$id}{'readSpt'}=0;
			$cp{$id}{'mut'}=0;
		}
		else
		{
			$cp{$id}{'length'}+=length($_);
		}
	}

	return %cp;
}

sub parse_map_file
{
	my ($in,$cp)=@_;

	open(IN,"$in") or die "Cannot open $in\n";

	while(<IN>)
	{
		my @tmp=split(/\t/,$_);
		$tmp[2] =~ s/hg19_knownGene_//;
		my $id=$tmp[2];
		if (defined($$cp{$id}))
		{
			$$cp{$id}{'readSpt'}++;
			if (defined($tmp[7]) and $tmp[7] =~ m/>/) { $$cp{$id}{'mut'}++; }
		}
		else {die "not existed gene: $id\n";}
		
	}

	close IN;
}

sub classify_CP
{
	my ($in,$cp)=@_;

	open(IN,"$in") or die "Cannot open $in\n";
	open(RSP,">RNAseq_spt_CP.fas");
	open(NSP,">No_RNAseq_spt_CP.fas");

	my (%rspt,%nspt,$line);
	foreach my $g (keys %{$cp})
	{
		if ($$cp{$g}{readSpt}) { $rspt{$g}=''; }
		else { $nspt{$g}=''; }
	}

	my $spt=-1;
	while($line=<IN>)
	{
		chomp($line);
		if ($line =~ m/^>/)
		{
			my @tmp=split(/\s/,$line);
                        $tmp[0] =~ s/>//;
                        my $id=$tmp[0];
			if (defined($rspt{$id})) { print RSP "$line\n"; $spt=1; }
			else { print NSP "$line\n"; $spt=0; }
		}
		elsif ($spt) { print RSP "$line\n";}
		else { print NSP "$line\n";}
	}

	close RSP;
	close NSP;
}
