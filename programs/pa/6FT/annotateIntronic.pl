#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=3)
{
        die "perl annotateIntronic.pl reads.alg.out.updated reads_vs_geneLocus.tab pep2readPos.txt\n";
}

# parse reads_vs_geneLocus.tab and build %geneLocusReads{$readID}{strand/ucscID}
my %geneLocusReads;
open(IN,"$ARGV[1]");
while(<IN>)
{
        next if (/^#/);
        my @t=split(/\t/,$_);
        my ($readID,$start,$end,$ucscID)=($t[0],$t[8],$t[9],$t[1]);

        # hg19_knownGene_uc010qoe.2
        $ucscID =~ s/hg19_knownGene_//;
        my $strand;
        if ( $start<$end ) { $strand='+'; }
        else { $strand='-'; }

        $geneLocusReads{$readID}{strand}=$strand;
        $geneLocusReads{$readID}{ucscID}=$ucscID;
}
close IN;

# read pep2readPos.txt; build %%pep2readPos{$pep}{$readID}{frame/pos/strand}
my (%pep2readPos);
open(IN,"$ARGV[2]");
while(<IN>)
{
        chomp;
        next if (/^peptide/);
        my ($pep,$readID,$pos,$frame)=split /\t/,$_;
        my ($pStrand,$offset)=split //,$frame;
        $pep2readPos{$pep}{$readID}{frame}=$frame;
        $pep2readPos{$pep}{$readID}{pos}=$pos;
        $pep2readPos{$pep}{$readID}{pStrand}=$pStrand;#print "$pStrand,";
}
close IN;

# parse reads.alg.out.updated row by row
open(IN,"$ARGV[0]");
while(<IN>)
{
        chomp;
	if (/^outfile/) { print  "$_\n";  next; }

	my $sType='';

        my @t=split(/\t/,$_);

        my ($outfile,$pep,$level,$np,$b5,$b3,$strand,$readID)=($t[0],$t[3],$t[14],$t[7],$t[18],$t[19],$t[20],$t[6]);
        next unless (defined($b5) and defined($b3) and defined($strand) and $b5 =~ m/^\d+$/ and $b3 =~ m/^\d+$/);
	$pep =~ s/[\*\@\$\^]//g;

        if ($level==5) # genome
        {
		if (defined($geneLocusReads{$readID}))
		{
			$sType='intronic';
                        $t[7]=$geneLocusReads{$readID}{ucscID};

			# determine strand
			my $pStrand=$pep2readPos{$pep}{$readID}{pStrand};
			my $sStrand=$geneLocusReads{$readID}{strand};
			my $strand='';
			if ($sStrand eq '+' and $pStrand eq 'F' ) { $strand='+'; }
			elsif ($sStrand eq '+' and $pStrand eq 'R' ) { $strand='-'; }
			elsif ($sStrand eq '-' and $pStrand eq 'F' ) { $strand='-'; }
			elsif ($sStrand eq '-' and $pStrand eq 'R' ) { $strand='+'; }
			$t[20]=$strand;
		}
		else { $sType='intergenic'; }
		$t[22]=$sType;
	}

	print  join("\t",@t),"\n";
}
close IN;
