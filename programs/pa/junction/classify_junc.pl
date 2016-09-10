#!/bin/env perl

use strict;
use warnings;

if ($#ARGV != 2)
{
        die "perl classify_junc.pl refflat.file frame.tab out\n";
}

#-------------------------------------------------
	my %annotations=(
		'coding'=>1, 
		'CDS'=>1,
		'5UTR'=>2, 
		'5utr'=>2,
		'3UTR'=>3, 
		'3utr'=>3,
		'intron'=>4,
		'intergenic'=>5
	);
my $floatingSplicing_dist=30; # nt
#-------------------------------------------------

# parse annotation file
my (%prohash,%geneStartEnd);
#open(IN,'/home/yli4/annotations/refSeq_mRNA2protein_072514.txt') or die "Cannot open annotatiion file!!!\n";
open(IN,$ARGV[0]) || die "cannot open $ARGV[0]\n";
while(<IN>)
{
        next if (/^#/);
        chomp;
        my @t=split(/\t/,$_);
        my ($nm,$chr,$txStart,$txEnd)=($t[1],$t[2],$t[4],$t[5]);
	next if ($chr =~ /\_/);

        if (defined($prohash{$nm}{$txStart}{$txEnd})) { die "duplicated record: {$nm}{$txStart}{$txEnd}!!!\n"; }
	else
        {
	        $prohash{$nm}{$txStart}{$txEnd}{gene}=$t[0];
        	$prohash{$nm}{$txStart}{$txEnd}{strand}=$t[3];
	        $prohash{$nm}{$txStart}{$txEnd}{chr}=$t[2];
        	$prohash{$nm}{$txStart}{$txEnd}{CDSstart}=$t[6];
	        $prohash{$nm}{$txStart}{$txEnd}{CDSend}=$t[7];
        	$prohash{$nm}{$txStart}{$txEnd}{line}=$_;

        	my @b5=split(/,/,$t[9]);
	        my @b3=split(/,/,$t[10]);
		for (my $i=0; $i<=$#b5; $i++) { $prohash{$nm}{$txStart}{$txEnd}{b5}[$i+1]=$b5[$i]+1; } # convert from 0-based to 1-based
		for (my $i=0; $i<=$#b3; $i++) { $prohash{$nm}{$txStart}{$txEnd}{b3}[$i+1]=$b3[$i]; }

		my $g=$t[0];  # gene name
		next if ($t[6]==$t[7]); # non-coding gene?
		if (!defined($geneStartEnd{$g}{CDSstart}) or $geneStartEnd{$g}{CDSstart}>$t[6]) { $geneStartEnd{$g}{CDSstart}=$t[6]; }
		if (!defined($geneStartEnd{$g}{CDSend}) or $geneStartEnd{$g}{CDSend}<$t[7]) { $geneStartEnd{$g}{CDSend}=$t[7]; }
		if (!defined($geneStartEnd{$g}{strand})) { $geneStartEnd{$g}{strand}=$t[3]; }
        }
}
close IN;


# parse frame.tab
open(IN,$ARGV[1]) || die "cannot open $ARGV[1]\n";
open(OUT,">$ARGV[2]");
open(OUT2,">isoformSpecific_$ARGV[2]");
#print OUT "junction\tstatusA\tstatusB\tannA\tannB\tframeA\tframeB\tstrand\tflankAA\tjuncPepClass\n";
print OUT "junction\tstrand\tgene\tgeneNameType\tjuncType\tgeneA\tgeneB\tbestAnnA\tbestAnnB\tdetailAnnA\tdetailAnnB\tbestDistA\tbestDistB\trelativePosTypeA\trelativePosTypeB\tframeA\tframeB\tpanCdsA\tpanCdsB\n";
print OUT2 "junction\tstrand\ttranscript\tannA\tannB\texonA\texonB\texonSideA\texonSideB\tdistA\tdistB\n";
while (<IN>)
{
        chomp;
	next if (/^junction/);
	my @t=split /\t/, $_;
	if ( !defined($t[31]) or $t[31] eq '' ) { $t[31]= '-'; }
        if ( !defined($t[32]) or $t[32] eq '' ) { $t[32]= '-'; }
	my ($junc,$strand,$frameA,$frameB,$flankAA,$annA,$annB)=($t[0],$t[2],$t[31],$t[32],$t[1],$t[33],$t[34]);

	# junc side position
	my ($posA,$posB);
	# 10:104140443:+,10:104141813:+
	$junc =~ m/(.*?)\:(\d+)\:\+\,.*?\:(\d+)\:\+/;
	if ($strand eq '+') { ($posA,$posB)=($2,$3); }
	else { ($posA,$posB)=($3,$2); }

	# build %transcripts
	# %transcripts{$NM_accession}{annA/annB}=coding/5utr/3utr/intron/intergenic
	my %transcripts;
	transcriptAnno($annA,\%transcripts,'A',$posA);
	transcriptAnno($annB,\%transcripts,'B',$posB);


	# Distance to known exon boundary: 
	$junc =~ m/(.*?)\:(\d+)\:\+\,.*?\:(\d+)\:\+/;
	my ($leftPos,$rightPos)=($2,$3);
	foreach my $nm (keys %transcripts)
	{
		my %tmphash;
		next if ($nm eq '' or $nm eq 'no_transcripts_found');
		if (!defined($prohash{$nm})) { die "$nm not defined in annotation file!!!\n"; }
		my ($txStart,$txEnd)=(0,0);
		foreach my $txStart1 (keys %{$prohash{$nm}})
		{
			foreach my $txEnd1 (keys %{$prohash{$nm}{$txStart1}})
			{
				if ($txStart1<=$leftPos and $leftPos<=$txEnd1
				or $txStart1<=$rightPos and $rightPos<=$txEnd1)
				{
					($txStart,$txEnd)=($txStart1,$txEnd1);
					last;
				}
			}
		}
		if ($txStart==0) { die "transcript $nm not in range for junction ($leftPos,$rightPos)\n"; }
		else
		{  #print "$nm,$txStart,$txEnd,$leftPos,$rightPos\n";
			annotateLeftJunSide($leftPos,$prohash{$nm}{$txStart}{$txEnd},$transcripts{$nm},\%tmphash);
			annotateRightJunSide($rightPos,$prohash{$nm}{$txStart}{$txEnd},$transcripts{$nm},\%tmphash);
			checkStrand($strand,$prohash{$nm}{$txStart}{$txEnd},$transcripts{$nm},\%tmphash);
		}
	}

	# pick best annotation
	my ($bestDistA,$bestAnnA,$detailAnnA)=bestAnnotation('geneA',\%transcripts);
	my ($bestDistB,$bestAnnB,$detailAnnB)=bestAnnotation('geneB',\%transcripts);

	foreach my $t (keys %transcripts)
	{ 
		if ($t eq '' or $t eq 'no_transcripts_found') 
		{ print OUT2 "$junc\t$strand\t$t\t$transcripts{$t}{geneA}{ann}\t$transcripts{$t}{geneB}{ann}\t-\t-\t-\t-\t-\t-\n"; }
		else
		{
			print OUT2 "$junc\t$strand\t$t\t$transcripts{$t}{geneA}{ann}\t$transcripts{$t}{geneB}{ann}\t$transcripts{$t}{geneA}{exon}\t$transcripts{$t}{geneB}{exon}\t$transcripts{$t}{geneA}{exonSide}\t$transcripts{$t}{geneB}{exonSide}\t$transcripts{$t}{geneA}{dist}\t$transcripts{$t}{geneB}{dist}\n"; 
		}
	}


	# geneNameAnnotation
	my ($geneNameType,$geneA,$geneB,$gene)=geneNameAnnotation($detailAnnA,$detailAnnB);

	# annotate junctions
	my $juncType=annotateJunctions($bestDistA,$bestDistB);

	# relativePosType
	my $relativePosTypeA=checkRelativePosType($geneA,$posA);
	my $relativePosTypeB=checkRelativePosType($geneB,$posB);

	# panCDS
	my $panCdsA=checkPanCDS($bestAnnA,$relativePosTypeA,$bestDistA,$posA,$geneA);
	my $panCdsB=checkPanCDS($bestAnnB,$relativePosTypeB,$bestDistB,$posB,$geneB);

	#print "$junc,$strand,$bestDistA,$bestAnnA,$detailAnnA,$bestDistB,$bestAnnB,$detailAnnB,$geneNameType,$geneA,$geneB,$gene,$juncType,$relativePosTypeA,$relativePosTypeB,$frameA,$frameB\n";
	print OUT "$junc\t$strand\t$gene\t$geneNameType\t$juncType\t$geneA\t$geneB\t$bestAnnA\t$bestAnnB\t$detailAnnA\t$detailAnnB\t$bestDistA\t$bestDistB\t$relativePosTypeA\t$relativePosTypeB\t$frameA\t$frameB\t$panCdsA\t$panCdsB\n";

}
close IN;
close OUT;
close OUT2;

#------------------------------------------
sub checkPanCDS
{
	my ($bestAnn,$relativePosType,$bestDist,$pos,$gene)=@_;

	my $panCds='no';

	if (defined($annotations{$bestAnn}))
	{
		if ($annotations{$bestAnn} == 1) { $panCds='yes'; } # coding
		elsif ( $annotations{$bestAnn} == 4) # intron
		{
			if ( $relativePosType eq 'internal' and abs($bestDist)<=$floatingSplicing_dist ) { $panCds='yes'; }
			elsif ( $relativePosType eq 'beforeStart' and abs($pos-$geneStartEnd{$gene}{CDSstart})<=$floatingSplicing_dist ) { $panCds='yes'; }
			elsif ( $relativePosType eq 'afterStop' and abs($pos-$geneStartEnd{$gene}{CDSend})<=$floatingSplicing_dist ) { $panCds='yes'; }
		}
		elsif ( $annotations{$bestAnn} == 2) # 5'UTR
		{
			if ( abs($pos-$geneStartEnd{$gene}{CDSstart})<=$floatingSplicing_dist ) { $panCds='yes'; }
		}
		elsif ( $annotations{$bestAnn} == 3) # 3'UTR
		{
			if ( abs($pos-$geneStartEnd{$gene}{CDSend})<=$floatingSplicing_dist ) { $panCds='yes'; }
		}
	}

	return $panCds;
}

sub checkRelativePosType
{
	my ($gene,$pos)=@_;

	my $type;
	# relativePosType: beforeStart/internal/afterStop
	if ( $gene eq '' or $gene eq '-' ) { $type='-'; }
	elsif ( $geneStartEnd{$gene}{CDSstart}<=$pos and $pos<=$geneStartEnd{$gene}{CDSend} ) { $type='internal'; }
	elsif ( $geneStartEnd{$gene}{strand} eq '+' and $geneStartEnd{$gene}{CDSstart}>$pos
	  or  $geneStartEnd{$gene}{strand} eq '-' and $geneStartEnd{$gene}{CDSend}<$pos) 	{ $type='beforeStart'; }
	elsif ( $geneStartEnd{$gene}{strand} eq '-' and $geneStartEnd{$gene}{CDSstart}>$pos
	  or  $geneStartEnd{$gene}{strand} eq '+' and $geneStartEnd{$gene}{CDSend}<$pos)	{ $type='afterStop'; }
	else { die "unknown relativePosType: $geneStartEnd{$gene}{CDSstart},$pos,$geneStartEnd{$gene}{CDSend}\n"; }
#if ($gene eq 'GMPR2') {print "$geneStartEnd{$gene}{strand} ,$geneStartEnd{$gene}{CDSstart},$pos,$geneStartEnd{$gene}{CDSend}\n";}
	return $type;
}

sub annotateJunctions
{
	my ($bestDistA,$bestDistB)=@_;

	my $type;
	if ($bestDistA eq '0' and $bestDistB eq '0') {$type='bothSidesKnown';}
	elsif ($bestDistA eq '0' or $bestDistB eq '0') {$type='oneSideNovel';}
	else {$type='bothSidesNovel';}

	return $type;
}

#sub classifyJunction
sub geneNameAnnotation
{
	my ($detailAnnA,$detailAnnB)=@_;

	my (%genesA,%genesB,%genes);
	getGeneNames($detailAnnA,\%genesA);
	getGeneNames($detailAnnB,\%genesB);

	foreach my $a (keys %genesA)
	{
		$genes{$a}{genenA}=1;
		$genes{$a}{genenB}=0;
	}

	foreach my $b (keys %genesB)
	{
		$genes{$b}{genenB}=1;
		if (!defined($genes{$b}{genenA})) {$genes{$b}{genenA}=0;}
	}

	my %targets;
	foreach my $g (keys %genes)
	{
		if ( $genes{$g}{genenA} +  $genes{$g}{genenB} == 2 ) { push @{$targets{both}},$g; }
		elsif ( $genes{$g}{genenA} ) { push @{$targets{A}},$g; }
		elsif ( $genes{$g}{genenB} ) { push @{$targets{B}},$g; }
	}

	my ($type,$geneA,$geneB,$gene);
	if (defined($targets{both}) and scalar(@{$targets{both}})) {$type='bothSidesInOneGene';$geneA=$geneB=$gene=shift(@{$targets{both}});}
	elsif (defined($targets{A}) and scalar(@{$targets{A}}) and defined($targets{B}) and scalar(@{$targets{B}})) {$type='acrossTwoGenes';$geneA=shift(@{$targets{A}});$geneB=shift(@{$targets{B}});$gene="$geneA,$geneB";}
	elsif (defined($targets{A}) and scalar(@{$targets{A}})) {$type='OneSideGene';$geneA=shift(@{$targets{A}});$geneB="-";$gene="$geneA,$geneB";}
	elsif (defined($targets{B}) and scalar(@{$targets{B}})) {$type='OneSideGene';$geneB=shift(@{$targets{B}});$geneA="-";$gene="$geneA,$geneB";}
	else {$type='intergenic';$geneA=$geneB=$gene="-";}

	return ($type,$geneA,$geneB,$gene);
}

sub getGeneNames
{
	my ($detailAnn,$genes)=@_;

	my @t=split(/\;/,$detailAnn);
	foreach my $a (@t)
	{
		next if ($a eq '-');
		my ($nm,$exon)=split(/\|/,$a);
		foreach my $a (keys %{$prohash{$nm}})
		{
			foreach my $b (keys %{$prohash{$nm}{$a}})
			{
				$$genes{$prohash{$nm}{$a}{$b}{gene}}='';
			}
		}
	}

}

sub bestAnnotation
{
	my ($side,$transcripts)=@_;

	unless ($side =~ /gene[AB]/) {die "Wrong side paramster in bestAnnotation: $side\n\n";}

	# 1st level: min dist
	my $minDist=10000000;my $bestDist='';
	foreach my $nm (keys %{$transcripts})
	{
		if (!defined($$transcripts{$nm}{$side}{ann}) or $$transcripts{$nm}{$side}{ann} eq '') {$$transcripts{$nm}{$side}{ann}='intergenic';}
		next if ($nm eq 'no_transcripts_found');
		next if ($$transcripts{$nm}{$side}{dist} eq '-');
		if ($minDist>abs($$transcripts{$nm}{$side}{dist})) {$minDist=abs($$transcripts{$nm}{$side}{dist});$bestDist=$$transcripts{$nm}{$side}{dist};}
	}
	if ($minDist==10000000) {return ('-','intergenic','-');} # intergenic
	my @targets;
	foreach my $nm (keys %{$transcripts})
	{
		next if ($nm eq 'no_transcripts_found');
		next if ($$transcripts{$nm}{$side}{dist} eq '-');
		if ($minDist==abs($$transcripts{$nm}{$side}{dist})) {push @targets,$nm;}
		#if ($nm eq 'NM_006720') {print "NM_006720:$minDist\n";}
	}

	# 2nd level: min annotation score
	my $minScore=100;my $bestAnn;
	foreach my $nm (@targets)
	{
		if (!defined($$transcripts{$nm}{$side}{ann}) or $$transcripts{$nm}{$side}{ann} eq '') {$$transcripts{$nm}{$side}{ann}='intergenic';}
		if (!defined($annotations{$$transcripts{$nm}{$side}{ann}})){die "not defined $$transcripts{$nm}{$side}{ann} in 2nd level: min annotation score\n\n";}
		if ($minScore>$annotations{$$transcripts{$nm}{$side}{ann}}) { $minScore=$annotations{$$transcripts{$nm}{$side}{ann}}; $bestAnn=$$transcripts{$nm}{$side}{ann};}
		#if ($nm eq 'NM_006720') {print "NM_006720:$minDist,$minScore,$bestAnn\n";}
	}
	my @targets2;
	foreach my $nm (@targets)
	{
		if ($minScore==$annotations{$$transcripts{$nm}{$side}{ann}})  {push @targets2,$nm;} 
	}

	# output
	my $detailAnn='';
	foreach my $nm (@targets2)
	{
		my $exonType=($side eq 'geneA')?3:5;
		my $s=($bestDist>0)?'+':'';
		$detailAnn.="$nm|exon".$$transcripts{$nm}{$side}{exon}."|b$exonType|$s$bestDist;";
		#$detailAnn.="$nm|exon".$$transcripts{$nm}{$side}{exon}."|b$exonType|$s$bestDist|$$transcripts{$nm}{$side}{relativePosType};";
	}
	return ($bestDist,$bestAnn,$detailAnn);
}

sub annotateLeftJunSide
{
	my ($pos,$hash,$transcripts,$tmphash)=@_;

	#my %tmphash;
	my $i=1;
	while ($pos>$$hash{b3}[$i]) {$i++; if ($i>=scalar(@{$$hash{b3}})) {die "Unexpected out of boundary for left side junction: $pos,$$hash{b3}[$i]";}}
	if  ($pos==$$hash{b3}[$i]) { $$tmphash{geneA}{exon}=$i; $$tmphash{geneA}{exonSide}=3; $$tmphash{geneA}{dist}=0; return; } # known ss
	#elsif ($i==scalar(@{$$hash{b3}})-1) { die "Unexpected out of boundary for left side junction: $pos,$$hash{b3}[$i]" } # more than last exon, should not be because it's left junction site
	elsif ($$hash{b5}[$i] <= $pos and $pos <= $$hash{b3}[$i]) { $$tmphash{geneA}{exon}=$i; $$tmphash{geneA}{exonSide}=3; $$tmphash{geneA}{dist}=$pos-$$hash{b3}[$i]; }
	elsif ($$hash{b5}[$i] > $pos)
	{
		if ($i==1) { $$tmphash{geneA}{exon}=0; $$tmphash{geneA}{exonSide}='-'; $$tmphash{geneA}{dist}='-'; return; } # beyond first exon: intergenic
		else { $$tmphash{geneA}{exon}=$i-1; $$tmphash{geneA}{exonSide}=3; $$tmphash{geneA}{dist}=$pos-$$hash{b3}[$i-1];; return; }
	}

}

sub annotateRightJunSide
{
	my ($pos,$hash,$transcripts,$tmphash)=@_;

	my $i=1;
	while ($i<scalar(@{$$hash{b3}}) and $pos>$$hash{b5}[$i]) {$i++;}

	if ($i==scalar(@{$$hash{b3}}))  # think about it
	{
		if ( $pos>$$hash{b3}[$i-1] ) { $$tmphash{geneB}{exon}=0; $$tmphash{geneB}{exonSide}='-'; $$tmphash{geneB}{dist}='-'; return; } # beyond last exon: intergenic 
		else { $$tmphash{geneB}{exon}=$i-1; $$tmphash{geneB}{exonSide}=5; $$tmphash{geneB}{dist}=$pos-$$hash{b5}[$i-1]; return; } # within last exon
	}

	if ( $pos==$$hash{b5}[$i] ) { $$tmphash{geneB}{exon}=$i; $$tmphash{geneB}{exonSide}=5; $$tmphash{geneB}{dist}=0; return; } # known ss
	elsif ($i==1) { die "Unexpected out of boundary for right side junction: $pos,$i,$$hash{b5}[$i]"; $$tmphash{geneB}{exon}=0; $$tmphash{geneB}{exonSide}='-'; $$tmphash{geneB}{dist}='-'; return; } # intergenic
	elsif ( $$hash{b3}[$i-1]<=$pos and $pos<$$hash{b5}[$i]) { $$tmphash{geneB}{exon}=$i; $$tmphash{geneB}{exonSide}=5; $$tmphash{geneB}{dist}=$pos-$$hash{b5}[$i]; return; } # upstream intron
	elsif ( $$hash{b3}[$i-1]>$pos ) { $$tmphash{geneB}{exon}=$i-1; $$tmphash{geneB}{exonSide}=5; $$tmphash{geneB}{dist}=$pos-$$hash{b5}[$i-1]; return; } # upstream exon
}

sub checkStrand
{
	my ($strand,$hash,$transcripts,$tmphash)=@_;

	if ($strand eq '+')
	{
		foreach my $a (keys %{$tmphash})
		{
			foreach my $b (keys %{$$tmphash{$a}}) { $$transcripts{$a}{$b}=$$tmphash{$a}{$b}; }
		}
	}
	else # negative strand
	{
		foreach my $a (keys %{$tmphash})
		{
			#$a='geneA'; 
			my $b=($a eq 'geneA')?'geneB':'geneA';
			if ($$tmphash{$a}{exon}>0)
			{
				$$transcripts{$b}{exon}=scalar(@{$$hash{b3}})-$$tmphash{$a}{exon}+1;
				$$transcripts{$b}{exonSide}=($$tmphash{$a}{exonSide}==5)?3:5;
				$$transcripts{$b}{dist}=0-$$tmphash{$a}{dist};
			}
			else
			{ foreach my $k (keys %{$$tmphash{$a}}) { $$transcripts{$b}{$k}=$$tmphash{$a}{$k}; } }
		}
	}
}

sub classifyJun
{
	my ($statusA,$statusB,$annA,$annB)=@_;

	my $juncClass='';

	# both CDS
	if ( $annotations{$annA}==1 and $annotations{$annB}==1 )
	{
		if ( $statusA eq 'inframe' and $statusB eq 'inframe') { $juncClass='bothinframe'; }
		elsif ( $statusA eq 'inframe' or $statusB eq 'inframe' ) { $juncClass='newORF'; }
		else { $juncClass='drop'; }
	}
	# one CDS, one 5utr
	elsif ( $annotations{$annA}==1 and $annotations{$annB}==2 
	or $annotations{$annA}==2 and $annotations{$annB}==1 )
	{
		if ( $statusA eq 'inframe' or $statusB eq 'inframe' ) { $juncClass='altStart'; }
		else { $juncClass='drop'; }
	}
	# one CDS, one 3utr
	elsif ( $annotations{$annA}==1 and $annotations{$annB}==3
	or $annotations{$annA}==3 and $annotations{$annB}==1 )
	{
		if ( $statusA eq 'inframe' or $statusB eq 'inframe' ) { $juncClass='altStop'; }
		else { $juncClass='drop'; }
	}
	# one CDS, one intron
	elsif ( $annotations{$annA}==1 and $annotations{$annB}==4
	or $annotations{$annA}==4 and $annotations{$annB}==1 )
	{
		if ( $statusA eq 'inframe' or $statusB eq 'inframe' ) { $juncClass='newCDSexon|intron'; }
		else { $juncClass='drop'; }
	}
	# one CDS, one intergenic
	elsif ( $annotations{$annA}==1 and $annotations{$annB}==5
	or $annotations{$annA}==5 and $annotations{$annB}==1 )
	{
		if ( $statusA eq 'inframe' or $statusB eq 'inframe' ) { $juncClass='newCDSexon|intergenic'; }
		else { $juncClass='drop'; }
	}

	return $juncClass;
}

sub frameStatus
{
        my ($frame,$ann)=@_;

        my $status;
        if ( $frame == 1 )
        {
                if ($annotations{$ann} == 1) { $status='inframe'; }
                elsif ( $ann =~ m/UTR/ or  $ann =~ m/utr/) { $status='UTRinframe';}
                else { $status='error'; }
        }
        else
        {
                if ( $annotations{$ann}!=1) { $status='unknown'; }
                else {  $status='notinframe';}
        }

        return $status;
}



sub transcriptAnno
{
	my ($ann,$transcripts,$AB,$pos)=@_;

	unless ($AB eq 'A' or $AB eq 'B') { die "wrong parameter for \$AB in transcriptAnno: $AB\n"; }
	my $side=($AB eq 'A')?'geneA':'geneB';

	#intron=NM_001270965,NM_001270966,NM_002779      
	#coding=NM_001270965,NM_001270966,NM_002779
	#coding=NM_001282280,NM_001282281,NM_006907,NM_153824;intron=NM_001282279
	my @types=split(/\;/,$ann);
	#my ($bestType,$bestCode)=('',10);
	foreach my $t (@types)
	{
		my ($a,$b)=split(/\=/,$t);
		if (!defined($annotations{$a})) {die "unexpected annotation types: $a\n";}
		#if ( $bestCode>$annotations{$a} ) { $bestCode=$annotations{$a}; $bestType=$a; }
		my @nm=split(/\,/,$b);
		foreach my $n (@nm)
		{
			$$transcripts{$n}{$side}{ann}=$a;

			# relativePosType: beforeStart/internal/afterStop
=head
			if ( $n eq 'no_transcripts_found' ) 
				{ $$transcripts{$n}{$side}{relativePosType}='-'; }
			elsif ( $prohash{$n}{CDSstart}<=$pos and $pos<=$prohash{$n}{CDSend} ) 
				{ $$transcripts{$n}{$side}{relativePosType}='internal'; }
			elsif ( $prohash{$n}{strand} eq '+' and $prohash{$n}{CDSstart}>$pos
			  or  $prohash{$n}{strand} eq '-' and $prohash{$n}{CDSend}<$pos)
				{ $$transcripts{$n}{$side}{relativePosType}='beforeStart'; }
			elsif ( $prohash{$n}{strand} eq '-' and $prohash{$n}{CDSstart}>$pos
			  or  $prohash{$n}{strand} eq '+' and $prohash{$n}{CDSend}<$pos)
				{ $$transcripts{$n}{$side}{relativePosType}='afterStop'; }
			else { die "unknown relativePosType: $prohash{$n}{CDSstart},$pos,$prohash{$n}{CDSend}\n"; }
=cut
		}
	}
}

sub checkFrame
{
        # check
        # 1) coding gene?
	# 2) return CDS exon start and end pos (relative to + strand genome)
        # 3) check 3 dividable
        my ($hash)=@_;

        my @t=split(/\t/,$$hash{line});
        if ($t[6]==$t[7]) {return 0;} # non-coding genes

        my @b5=split(/,/,$t[9]);
        my @b3=split(/,/,$t[10]);

        my @exonl; # exon lengths
        for (my $i=0;$i<=$#b5; $i++) { $exonl[$i]=$b3[$i]-$b5[$i]; }

        my $i=0;
        while ($t[6]>$b3[$i]) {$i++;} my $j=0; $$hash{eb5}[$j]=$t[6]; # 1st CDS exon
        while ($t[7]>$b3[$i])
        {
                $$hash{eb3}[$j]=$b3[$i];
                $i++;
                $j++;
                $$hash{eb5}[$j]=$b5[$i];
        }
        $$hash{eb3}[$j]=$t[7]; # last CDS exon

        for ($i=0;$i<scalar(@{$$hash{eb3}}); $i++) { $$hash{elength}[$i]=$$hash{eb3}[$i]-$$hash{eb5}[$i]; }
        my $cdsL=0;
        for ($i=0;$i<scalar(@{$$hash{eb3}}); $i++) { $cdsL+=$$hash{elength}[$i]; $$hash{cmllength}[$i]=$cdsL; }
        $$hash{cdslength}=$cdsL;
        if ( $cdsL % 3 == 0 )  { return 1;}
        else {return -1;}
}

