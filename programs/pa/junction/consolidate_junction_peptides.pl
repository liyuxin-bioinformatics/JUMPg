#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/g

use IDtxt_parser;
use File::Basename;
use strict;
use warnings;

if (scalar(@ARGV)!=4)
{
        die "perl consolidate_junction_peptides.pl junc_ID.txt novel_AA.fas AApadding out_prefix\n";
}

my $AApadding=$ARGV[2];
my $ouput_prefix=$ARGV[3];


my (%prohash);
open(DB,"$ARGV[1]") || die "cannot open novel_AA.fas: $ARGV[1]!!!\n";
while(<DB>)
{
        chomp;
        s/^>//;
        my ($pro,$ann)=split(/ /,$_);
        $prohash{$pro}{annotation}=$ann;
        my ($junc,$frame)=split(/\|/,$pro);
	# Gene=DTNBP1|Annotation_Set=extended
	my $gene=(split /\|/,$ann)[0];
	$gene =~ s/^Gene=//;
	$prohash{$pro}{gene}=$gene;

        my $seq=<DB>;
        chomp($seq);
        $prohash{$pro}{seq}=$seq;
        $prohash{$pro}{stopN}=$seq =~ tr/*/*/;

        #$prohash{$pro}{juncPointU}=int(length($seq)/2+0.5);
        #if ($prohash{$pro}{juncPointU} == length($seq)/2) { $prohash{$pro}{juncPointD}=$prohash{$pro}{juncPointU}+1; }

        #if ($frame eq 'R0' || $frame eq 'F0' || $frame eq 'F1') { $prohash{$pro}{juncPointU}=$prohash{$pro}{juncPointD}=23; }
        #elsif ($frame eq 'R1' || $frame eq 'F2') { $prohash{$pro}{juncPointU}=22;$prohash{$pro}{juncPointD}=23; }
        #elsif ($frame eq 'R2') { $prohash{$pro}{juncPointU}=$prohash{$pro}{juncPointD}=22; }
        if ($frame eq 'R0' || $frame eq 'F0') {
		$prohash{$pro}{juncPointU}=$AApadding;
		$prohash{$pro}{juncPointD}=$AApadding+1;
	} else {
		$prohash{$pro}{juncPointU}=$prohash{$pro}{juncPointD}=$AApadding;
	}
}
close DB;

my $idp=IDtxt_parser->new();
my ($scanhash,$peptidehash)=$idp->parse_IDtxt($ARGV[0]);

#print scalar(keys %$scanhash),"\n";
#print scalar(keys %$peptidehash);

my $junctionhash;
foreach my $pep (keys %$peptidehash)
{
	my $nomod=$idp->nomod_pep($pep);
	foreach my $pro (keys %{$peptidehash->{$pep}->{proteins}})
	{
		# 10:104167023:+,10:104168652:+|R1
		my ($junc,$frame)=split /\|/,$pro;
		$peptidehash->{$pep}->{junctions}->{$junc}->{frame}=$frame;

		# AA pos in protein
		my ($proseq,$juncP);
		if (defined($prohash{$pro})) { $proseq=$prohash{$pro}{seq}; }
		else { print "Not defnined protein: $pro!!!\n"; next;}

		my $fpos = index($proseq, $pep) + 1;
		my $lpos = $fpos + length($pep) - 1;

		# check juncPoint
		if ( $fpos<=$prohash{$pro}{juncPointU} and $prohash{$pro}{juncPointU}<=$lpos
	        and $fpos<=$prohash{$pro}{juncPointD} and $prohash{$pro}{juncPointD}<=$lpos) { $juncP=1; }
        	else { $juncP=0; }
		$peptidehash->{$pep}->{junctions}->{$junc}->{containSS}=$juncP;

		# currently only consider peptides that contain SS
		# But keep in mind that if the novel AS generates frame shifted peptides, even peptide not containing SS could be great evidence
		next unless ($juncP); 

		# junctionhash
		$junctionhash->{$junc}->{peptides}->{$pep}->{containSS}=$juncP;
		if (!defined($junctionhash->{$junc}->{SC}))
		{
			$junctionhash->{$junc}->{SC}=0;
			#$junctionhash->{$junc}->{pepN}=0;
		}
		$junctionhash->{$junc}->{peptides}->{$pep}=scalar(keys %{$peptidehash->{$pep}->{outfiles}}); # SC for this peptide
		$junctionhash->{$junc}->{SC}+=$junctionhash->{$junc}->{peptides}->{$pep};
		if (!defined($junctionhash->{$junc}->{XCorr}) || $junctionhash->{$junc}->{XCorr}<$peptidehash->{$pep}->{XCorr})
		{
			$junctionhash->{$junc}->{XCorr}=$peptidehash->{$pep}->{XCorr};
			$junctionhash->{$junc}->{top_score_outfile}=$peptidehash->{$pep}->{top_score_outfile};
		}
		$junctionhash->{$junc}->{gene}=$prohash{$pro}{gene};
	}

}

# output
# currently only print peptides that contain SS
# But keep in mind that if the novel AS generates frame shifted peptides, even peptide not containing SS could be great evidence
open(JUN,">$ouput_prefix\_events.txt");
print JUN "junction\tgene\tpeptide_counts\tpeptide_sequences\tscan_counts\tbest_scan\n";
foreach my $junc (keys %$junctionhash)
{
	print JUN "$junc\t$junctionhash->{$junc}->{gene}\t",scalar(keys %{$junctionhash->{$junc}->{peptides}}),"\t";
	my @peptides;
	foreach my $pep (keys %{$junctionhash->{$junc}->{peptides}})
	{
		push @peptides, $pep;
	}
	my $out=basename($junctionhash->{$junc}->{top_score_outfile});
	print JUN join(",",@peptides),"\t$junctionhash->{$junc}->{SC}\t$out\n";
}
close JUN;

# currently only print peptides that contain SS
open(PEP,">$ouput_prefix\_peptides.txt");
print PEP "peptide\tscan_counts\tbest_scan\tsearch_score\tdelta_score\tjunction\n";
foreach my $pep (keys %$peptidehash)
{
	foreach my $junc (keys %{$peptidehash->{$pep}->{junctions}})
	{
		if ($peptidehash->{$pep}->{junctions}->{$junc}->{containSS})
		{
			my $out=basename($peptidehash->{$pep}->{top_score_outfile});
			# peptides could be redundant if one peptide map to multiple junctions
			print PEP "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$out\t$peptidehash->{$pep}->{XCorr}\t$peptidehash->{$pep}->{dCn}\t$junc\n";
		}
	}
}
close PEP;

# accepted ID.txt
# valideted_ID.txt
open(OUT,">valideted_ID.txt");
foreach my $pep (keys %$peptidehash)
{
	foreach my $junc (keys %{$peptidehash->{$pep}->{junctions}})
	{
		if ($peptidehash->{$pep}->{junctions}->{$junc}->{containSS})
		{
			foreach my $Outfile (keys %{$peptidehash->{$pep}->{outfiles}})
			{
				print OUT "$scanhash->{$Outfile}->{IDtxt_line}\n";
			}
		}
	}
}
close OUT;

# output another file as junction_frameCheck_input.pl input that countains 'frame'
# currently only print peptides that contain SS
open(PEP,">$ouput_prefix\_peptides_frame.txt");
print PEP "peptide\tscan_counts\tbest_scan\tsearch_score\tdelta_score\tjunction\tjuncFrame\n";
foreach my $pep (keys %$peptidehash)
{
	foreach my $junc (keys %{$peptidehash->{$pep}->{junctions}})
	{
		if ($peptidehash->{$pep}->{junctions}->{$junc}->{containSS})
		{
			my $out=basename($peptidehash->{$pep}->{top_score_outfile});
			# peptides could be redundant if one peptide map to multiple junctions
			my $juncFrame=$junc;
			$juncFrame.="\|$peptidehash->{$pep}->{junctions}->{$junc}->{frame}";
			print PEP "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$out\t$peptidehash->{$pep}->{XCorr}\t$peptidehash->{$pep}->{dCn}\t$junc\t$juncFrame\n";
		}
	}
}
close PEP;
