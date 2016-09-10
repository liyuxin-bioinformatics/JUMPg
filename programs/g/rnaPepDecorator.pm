#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.2.6/programs/g
package rnaPepDecorator;
use strict;
use PrimarySeq;
sub new
{
        my($class) = @_;
        my $self = {};
        bless ($self,$class);
        return $self;
}

#---------------------------------------------------------------------------------------------------------

1;

sub publicTab2 {
	shift @_;
	my ($scanhash,$peptidehash, $output)=@_;
	$output||='id_peptide.txt';
	open(OUT,">",$output);
	#print OUT "peptide\tPSMs\tbest_PSM\tJscore\tdJn\t\tread_counts\tgene\tlocation_type\tframe\tgenome_chr\tgenome_strand\tgenome_start\tgenome_end\n";
	print OUT "peptide\tPSMs\tbest_PSM\tJscore\tdJn\tread_counts\tgene\tlocation_type\tgenome_chr\tgenome_strand\tgenome_start\tgenome_end\tGenomicSeq\tref_peptide\tAAdiff\n";
	foreach my $out (keys %$scanhash) {
		my $pep=$scanhash->{$out}->{peptide};
		next unless ($out eq $peptidehash->{$pep}->{top_score_outfile});

		my $rnaID=$scanhash->{$out}->{rnaID};
		if (defined($scanhash->{$out}->{RNA}->{$rnaID}->{b5})) {
			my $blkN=scalar(@{$scanhash->{$out}->{RNA}->{$rnaID}->{b5}});

			print OUT "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$out\t$scanhash->{$out}->{XCorr}\t$scanhash->{$out}->{dCn}\tnot_available\t$scanhash->{$out}->{transcriptID}\t$scanhash->{$out}->{functionalType}\t$scanhash->{$out}->{RNA}->{$rnaID}->{chr}\t$scanhash->{$out}->{RNA}->{$rnaID}->{genomicStrand}\t$scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[0]\t$scanhash->{$out}->{RNA}->{$rnaID}->{b3}->[$blkN-1]\t";
			if (defined($scanhash->{$out}->{ref_peptide})) {
				print OUT "$scanhash->{$out}->{GenomicSeq}\t$scanhash->{$out}->{ref_peptide}\t$scanhash->{$out}->{AAdiff}";
			}
			print OUT "\n";
		} else {
			print OUT "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$out\t$scanhash->{$out}->{XCorr}\t$scanhash->{$out}->{dCn}\tnot_available\t$scanhash->{$out}->{transcriptID}\t$scanhash->{$out}->{functionalType}\tNA\tNA\tNA\tNA\n";
		}
	}
	close OUT;
}

sub annCDSpep {
	#shift @_;
	my ($self,$scanhash,$peptidehash,$genomeSeq)=@_;

	foreach my $out (keys %$scanhash) {
		my $pep=$scanhash->{$out}->{peptide};
		next unless ($out eq $peptidehash->{$pep}->{top_score_outfile});
		# only for CDS peptides
		next unless ($scanhash->{$out}->{functionalType} eq 'CDS');

		my $rnaID=$scanhash->{$out}->{rnaID};
		my $chr=$scanhash->{$out}->{RNA}->{$rnaID}->{chr};
		my $strand=$scanhash->{$out}->{RNA}->{$rnaID}->{genomicStrand};
		if (defined($scanhash->{$out}->{RNA}->{$rnaID}->{b5})) {
			# get genomic sequnce via peptide genomic positions
			my $pepGenomicSeq=$self->getGenomic($genomeSeq->{$chr},$scanhash->{$out}->{RNA}->{$rnaID});
			# translate genomic sequnce to AA sequence
			my $ps=PrimarySeq->new;
			my $refRNAseq=($strand eq '+')?$pepGenomicSeq:$ps->revcom($pepGenomicSeq);
			my $AAseq=$ps->translate($refRNAseq);
			# count mismatches by comparing identified peptide sequence with ref AA sequence
			my $AAdiff=$ps->sequenceDiff($scanhash->{$out}->{nomod},$AAseq);
			# store information in %scanhash
			$scanhash->{$out}->{GenomicSeq}=$pepGenomicSeq;
			$scanhash->{$out}->{ref_peptide}=$AAseq;
			$scanhash->{$out}->{AAdiff}=$AAdiff;
			
		} else {
		}
	}
	
}

sub getGenomic {
	shift @_;
	my ($chrSeq,$locationhash)=@_;

	my $blkN=scalar(@{$locationhash->{b5}});
	my $resultSeq='';

	for (my $i=0; $i<$blkN; $i++) {
		my $blkSize=$locationhash->{b3}->[$i] - $locationhash->{b5}->[$i];
		$resultSeq.= substr $chrSeq, $locationhash->{b5}->[$i], $blkSize;
	}

	return $resultSeq;
}

sub publicTabPSM {
	shift @_;
	my ($scanhash,$output)=@_;
	$output||='id_PSM.txt';
	open(OUT,">",$output);
	print OUT "outfile\tpeptide\tJscore\tdJn\tppm\tgene\tlocation_type\n";
	foreach my $out (keys %$scanhash) {
		my $pep=$scanhash->{$out}->{peptide};
		
		print OUT "$out\t$pep\t$scanhash->{$out}->{XCorr}\t$scanhash->{$out}->{dCn}\t$scanhash->{$out}->{ppm}\t$scanhash->{$out}->{transcriptID}\t$scanhash->{$out}->{functionalType}\n";
	}
	close OUT;
}

sub publicTab {
	shift @_;
	my ($scanhash,$peptidehash, $output)=@_;
	$output||='id_peptide.txt';
	open(OUT,">",$output);
	#print OUT "peptide\tPSMs\tbest_PSM\tJscore\tdJn\t\tread_counts\tgene\tlocation_type\tframe\tgenome_chr\tgenome_strand\tgenome_start\tgenome_end\n";
	print OUT "peptide\tPSMs\tbest_PSM\tJscore\tdJn\tread_counts\tgene\tlocation_type\tgenome_chr\tgenome_strand\tgenome_start\tgenome_end\n";
	foreach my $out (keys %$scanhash) {
		my $pep=$scanhash->{$out}->{peptide};
		next unless ($out eq $peptidehash->{$pep}->{top_score_outfile});

		my $rnaID=$scanhash->{$out}->{rnaID};
		if (defined($scanhash->{$out}->{RNA}->{$rnaID}->{b5})) {
			my $blkN=scalar(@{$scanhash->{$out}->{RNA}->{$rnaID}->{b5}});

			print OUT "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$out\t$scanhash->{$out}->{XCorr}\t$scanhash->{$out}->{dCn}\tnot_available\t$scanhash->{$out}->{transcriptID}\t$scanhash->{$out}->{functionalType}\t$scanhash->{$out}->{RNA}->{$rnaID}->{chr}\t$scanhash->{$out}->{RNA}->{$rnaID}->{genomicStrand}\t$scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[0]\t$scanhash->{$out}->{RNA}->{$rnaID}->{b3}->[$blkN-1]\n";
		} else {
			print OUT "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$out\t$scanhash->{$out}->{XCorr}\t$scanhash->{$out}->{dCn}\tnot_available\t$scanhash->{$out}->{transcriptID}\t$scanhash->{$out}->{functionalType}\tNA\tNA\tNA\tNA\n";
		}
	}
	close OUT;
}

sub pepFunRegion {
	shift @_;
	my ($scanhash,$transcripthash)=@_;

	foreach my $hash (values %$scanhash) {
		#foreach my $hash2 (values %{$hash->{RNA}}) {
		foreach my $rnaID (keys %{$hash->{RNA}}) {
			my $hash2=$hash->{RNA}->{$rnaID};
			next unless (defined($hash2->{transcripts}));

			my $chr=$hash2->{chr};
			my $blkN=scalar(@{$hash2->{b5}});
			my $pepStart=$hash2->{b5}->[0];
			my $pepEnd=$hash2->{b3}->[$blkN-1];
			my $b5=$hash2->{b5};
			my $b3=$hash2->{b3};

			foreach my $id (keys %{$hash2->{transcripts}}) {
				my $thash=$transcripthash->{$chr}->{$id};
				my $frEnds=$thash->{frEnds};
				my $frTypes=$thash->{frTypes};
				$hash2->{transcripts}->{$id}->{type}='';

				# near gene?
				if ($pepStart < $thash->{trStart} or $pepEnd > $thash->{trEnd}) {
					$hash2->{transcripts}->{$id}->{type}='nearGene';
					next;
				}

				# ncGene?
				if ($thash->{cdsStart} == $thash->{cdsEnd}) {
					$hash2->{transcripts}->{$id}->{type}='ncGene';
					next;
				}

				# check each block: any non-CDS?
				my $i=0; # index for transcript fragments
				for (my $k=0; $k<$blkN; $k++) {
					# anchor block start
					while ($b5->[$k] >= $frEnds->[$i]) {
						$i++;
					}

					# not in CDS?
					if ($frTypes->[$i] ne 'CDS') {
						$hash2->{transcripts}->{$id}->{type}=$frTypes->[$i];
						last;
					} else { # block starts in CDS fragment
						if ($b3->[$k] <= $frEnds->[$i]) { # block ends in CDS fragment
							next; # check next block
						} else {
							$hash2->{transcripts}->{$id}->{type}=$frTypes->[$i+1]; 
						}
					}
				}

				# all blocks checked: all are CDS
				if ($hash2->{transcripts}->{$id}->{type} eq '') {
					$hash2->{transcripts}->{$id}->{type}='CDS';
				}
			}
		}
	}

	# select representative transcript according to type: CDS > utr > intron > ncGene > nearGene > n/a (intergenic)
	foreach my $hash (values %$scanhash) {
		my %typehash;
		my %rnaIDhash;
		#foreach my $hash2 (values %{$hash->{RNA}}) {
		foreach my $rnaID (keys %{$hash->{RNA}}) {
			my $hash2=$hash->{RNA}->{$rnaID};
			if (!defined($hash2->{b5})) {
				$typehash{'unmapped'}='not_available';
				$rnaIDhash{'unmapped'}='not_available';
				next;
			}
			if (!defined($hash2->{transcripts})) {
				$typehash{'intergenic'}='not_available';
				$rnaIDhash{'intergenic'}=$rnaID;
				next;
			}
			foreach my $id (keys %{$hash2->{transcripts}}) {
				$typehash{$hash2->{transcripts}->{$id}->{type}}=$id; # could use gene name in future?
				$rnaIDhash{$hash2->{transcripts}->{$id}->{type}}=$rnaID; # Trnity ID
			}
		}
		# select representative transcript according to type
		if (defined($typehash{'CDS'})) {
			$hash->{functionalType}='CDS';
			$hash->{transcriptID}=$typehash{'CDS'};
			$hash->{rnaID}=$rnaIDhash{'CDS'};
		} elsif (defined($typehash{'UTR'})) {
			$hash->{functionalType}='UTR';
			$hash->{transcriptID}=$typehash{'UTR'};
			$hash->{rnaID}=$rnaIDhash{'UTR'};
		} elsif (defined($typehash{'intron'})) {
			$hash->{functionalType}='intron';
			$hash->{transcriptID}=$typehash{'intron'};
			$hash->{rnaID}=$rnaIDhash{'intron'};
		} elsif (defined($typehash{'ncGene'})) {
			$hash->{functionalType}='ncGene';
			$hash->{transcriptID}=$typehash{'ncGene'};
			$hash->{rnaID}=$rnaIDhash{'ncGene'};
		} elsif (defined($typehash{'nearGene'})) {
			$hash->{functionalType}='nearGene';
			$hash->{transcriptID}=$typehash{'nearGene'};
			$hash->{rnaID}=$rnaIDhash{'nearGene'};
		} elsif (defined($typehash{'intergenic'})) {
			$hash->{functionalType}='intergenic';
			$hash->{transcriptID}=$typehash{'intergenic'};
			$hash->{rnaID}=$rnaIDhash{'intergenic'};
		} elsif (defined($typehash{'unmapped'})) {
			$hash->{functionalType}='unmapped';
			$hash->{transcriptID}=$typehash{'unmapped'};
			$hash->{rnaID}=$rnaIDhash{'unmapped'};
		}
	}
}

sub pep2gene {
	shift @_;
	my ($scanhash,$transcripthash)=@_;

	# %scanhash{$out}{RNA}{$rnaID}{frame/ORFstart/ORFend/seq/pep2ORFstart/pep2ORFend/pep2rnaStart/pep2rnaEnd}
	# @{$scanhash{$out}{RNA}{$rnaID}{b5/b3}}
	foreach my $hash (values %$scanhash) {
		foreach my $hash2 (values %{$hash->{RNA}}) {
			next unless (defined($hash2->{b5}));
			my $blkN=scalar(@{$hash2->{b5}});
			my $pepStart=$hash2->{b5}->[0];
			my $pepEnd=$hash2->{b3}->[$blkN-1];
			my $chr=$hash2->{chr};

			foreach my $id (keys %{$transcripthash->{$chr}}) {
				if ($transcripthash->{$chr}->{$id}->{trStart}<=$pepStart and $pepStart<=$transcripthash->{$chr}->{$id}->{trEnd}
				or $transcripthash->{$chr}->{$id}->{trStart}<=$pepEnd and $pepEnd<=$transcripthash->{$chr}->{$id}->{trEnd} ) {
					$hash2->{transcripts}->{$id}->{trStart}=$transcripthash->{$chr}->{$id}->{trStart};
					$hash2->{transcripts}->{$id}->{trEnd}=$transcripthash->{$chr}->{$id}->{trEnd};
				}
			}
		}
	}
}

sub parseRefFlat {
	shift @_;
	my ($refFlat)=@_;
	my $transcripthash;

	open(IN,$refFlat) || die "cannot open ID.txt: $refFlat\n";
	while (<IN>) {
		chomp;
		my ($id,$chr,$strand,$trStart,$trEnd,$cdsStart,$cdsEnd,$exonN,$exonStarts,$exonEnds,$proID,$id2)=split /\t/,$_;
		$transcripthash->{$chr}->{$id}->{trStart}=$trStart;
		$transcripthash->{$chr}->{$id}->{trEnd}=$trEnd;
		$transcripthash->{$chr}->{$id}->{strand}=$strand;
		$transcripthash->{$chr}->{$id}->{cdsStart}=$cdsStart;
		$transcripthash->{$chr}->{$id}->{cdsEnd}=$cdsEnd;

		# ncGene?
		my $ncGene=($cdsStart==$cdsEnd)?1:0;

		# fill in @frEnds / @frTypes
		my ($frEnds,$frTypes);
		my $k=0; # index of @frEnds / @frTypes
		my @b5=split /,/,$exonStarts;
		my @b3=split /,/,$exonEnds;
		# 1) exon that contains CDS start / end: split into two fragments
		# 2) check if last exon:
		#  a) no: one exon, two fragments (except start/stop codon exons): exon itself + intron on the 3'side
		#  b) yes: only one fragment (exon itself; no upcoming intron)
		my $frtype=($ncGene)?'ncGene':'UTR';
		for (my $i=0; $i<$exonN; $i++) {
			if ($b5[$i]<$cdsStart and $cdsStart<$b3[$i]) { # exon with CDS start?
				$frEnds->[$k]=$cdsStart;
				$frTypes->[$k]=$frtype;
				$k++;
				$frtype='CDS';
			} elsif ($b5[$i]==$cdsStart and !$ncGene) { # start codon happens to be the 1st nucleotide
				$frtype='CDS';
			}

			if ($b5[$i]<$cdsEnd and $cdsEnd<=$b3[$i]) { # exon with CDS end?
				$frEnds->[$k]=$cdsEnd;
				$frTypes->[$k]=$frtype;
				$k++;
				$frtype='UTR';
			}

			# normal exons
			unless ($cdsEnd==$b3[$i]) {
				$frEnds->[$k]=$b3[$i];
				$frTypes->[$k]=$frtype;
				$k++;
			}

			# unless last exon, assign intron
			unless ($i==$#b5) {
				$frEnds->[$k]=$b5[$i+1];
				$frTypes->[$k]='intron';
				$k++;
			}
		}
		$transcripthash->{$chr}->{$id}->{frEnds}=$frEnds;
		$transcripthash->{$chr}->{$id}->{frTypes}=$frTypes;
	}
	close IN;

	return $transcripthash;
}

sub printBED {
	shift @_;
	my ($scanhash,$output)=@_;

	open(OUT,">",$output);
	foreach my $out (keys %$scanhash) {
		foreach my $rnaID (keys %{$scanhash->{$out}->{RNA}}) {
			next unless (defined($scanhash->{$out}->{RNA}->{$rnaID}->{b5}));
			my $blkN=scalar(@{$scanhash->{$out}->{RNA}->{$rnaID}->{b5}});
			print OUT $scanhash->{$out}->{RNA}->{$rnaID}->{chr},"\t",
				$scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[0],"\t",
				$scanhash->{$out}->{RNA}->{$rnaID}->{b3}->[$blkN-1],"\t",
				$scanhash->{$out}->{peptide},"\t900\t",
				$scanhash->{$out}->{RNA}->{$rnaID}->{genomicStrand},"\t",
				$scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[0],"\t",
				$scanhash->{$out}->{RNA}->{$rnaID}->{b3}->[$blkN-1],"\t0\t",
				$blkN,"\t";
			# print block size
			for (my $i=0; $i<$blkN; $i++) {
				print OUT $scanhash->{$out}->{RNA}->{$rnaID}->{b3}->[$i] - $scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[$i],",";
			}
			print OUT "\t";
			# print block starts
			for (my $i=0; $i<$blkN; $i++) {
				print OUT $scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[$i] - $scanhash->{$out}->{RNA}->{$rnaID}->{b5}->[0],",";
			}
			print OUT "\t$out\n";
		}
	}
	close OUT;
}

sub pepGenomePos {
	shift @_;
	my ($scanhash,$rnaalignhash)=@_;

	# %scanhash{$out}{RNA}{$rnaID}{frame/ORFstart/ORFend/seq/pep2ORFstart/pep2ORFend/pep2rnaStart/pep2rnaEnd}
	# @{$scanhash{$out}{RNA}{$rnaID}{b5/b3}}
	foreach my $hash (values %$scanhash) {
		foreach my $rnaID (keys %{$hash->{RNA}}) {
			if (defined($rnaalignhash->{$rnaID})) {

				# defined ($b5,$b3): relative to genome of postive strand 
				my ($b5,$b3);
				if ($rnaalignhash->{$rnaID}->{strand} eq '+') {
					$b5=$hash->{RNA}->{$rnaID}->{pep2rnaStartAbs};
					$b3=$hash->{RNA}->{$rnaID}->{pep2rnaEndAbs};
				} else {
					$b5= $hash->{RNA}->{$rnaID}->{rnalength} - $hash->{RNA}->{$rnaID}->{pep2rnaEndAbs};
					$b3= $hash->{RNA}->{$rnaID}->{rnalength} - $hash->{RNA}->{$rnaID}->{pep2rnaStartAbs};
				}

				# check if peptide is within 'unaligned' part
				next if ($b5 < $rnaalignhash->{$rnaID}->{cmllength}->[0] # peptide in 5' unaligned part
				or $b3 > $rnaalignhash->{$rnaID}->{cmllength}->[$rnaalignhash->{$rnaID}->{blockCount}]); # peptide in 3' unaligned part

				# chr
				$hash->{RNA}->{$rnaID}->{chr}=$rnaalignhash->{$rnaID}->{chr}; 

				# peptide  strand on genome
				if ( $hash->{RNA}->{$rnaID}->{strand} eq $rnaalignhash->{$rnaID}->{strand} ) {
					$hash->{RNA}->{$rnaID}->{genomicStrand}='+';
				} else {
					$hash->{RNA}->{$rnaID}->{genomicStrand}='-';
				}

				# locate the exon in which pep2rnaStartAbs is in
				my $i=1;
				while ( $b5 >= $rnaalignhash->{$rnaID}->{cmllength}->[$i]) {
					$i++;
				}

				# genomic position of pep2rnaStartAbs
				my $j=0;
				$hash->{RNA}->{$rnaID}->{b5}->[$j] = $rnaalignhash->{$rnaID}->{eb3}->[$i] - ($rnaalignhash->{$rnaID}->{cmllength}->[$i]-$b5);

				# fill in @{$scanhash{$out}{RNA}{$rnaID}{b5/b3}} [note: may cross multiple introns]
				while ( $b3>$rnaalignhash->{$rnaID}->{cmllength}->[$i] ) {

					$hash->{RNA}->{$rnaID}->{b3}->[$j] = $rnaalignhash->{$rnaID}->{eb3}->[$i];
					$i++; $j++;
					$hash->{RNA}->{$rnaID}->{b5}->[$j] = $rnaalignhash->{$rnaID}->{eb5}->[$i];
				}

				# genomic position of pep2rnaEndAbs
				$hash->{RNA}->{$rnaID}->{b3}->[$j] = $rnaalignhash->{$rnaID}->{eb3}->[$i] - ($rnaalignhash->{$rnaID}->{cmllength}->[$i]-$b3);
			}
		}
	}
}

sub parseBlat {
	shift @_;
	my ($aligmentFile)=@_;

	# %rnaalignhash{$rnaID}{strand/chr/line/@eb5/@eb3/@elength/@cmllength}
	my $rnaalignhash;

	# load alignment
	open(IN,$aligmentFile) || die "cannot open ID.txt: $aligmentFile\n";
	while (<IN>) {
		chomp; 
		my $line=$_;
		my ($matches,$misMatches ,$repMatches ,$nCount , $qNumInsert ,$qBaseInsert ,$tNumInsert ,$tBaseInsert ,$strand ,$qName,$qSize ,$qStart ,$qEnd  ,$tName,$tSize ,$tStart ,$tEnd ,$blockCount ,$blockSizes,$qStarts  ,$tStarts )=split /\t/,$_;

		next if ($qNumInsert>0); # if gap in query? discard in this version!

		if (!defined($rnaalignhash->{$qName}) or $rnaalignhash->{$qName}->{matches} < $matches) {
			$rnaalignhash->{$qName}->{"strand"}=$strand;
			$rnaalignhash->{$qName}->{"matches"}=$matches;
			$rnaalignhash->{$qName}->{"tieScore"}=1;
			$rnaalignhash->{$qName}->{"chr"}=$tName;
			$rnaalignhash->{$qName}->{"line"}=$line;

			$rnaalignhash->{$qName}->{"misMatches"}=$misMatches;
			$rnaalignhash->{$qName}->{"repMatches"}=$repMatches;
			$rnaalignhash->{$qName}->{"nCount"}=$nCount;
			$rnaalignhash->{$qName}->{"qNumInsert"}=$qNumInsert;
			$rnaalignhash->{$qName}->{"qBaseInsert"}=$qBaseInsert;
			$rnaalignhash->{$qName}->{"tNumInsert"}=$tNumInsert;
			$rnaalignhash->{$qName}->{"blockCount"}=$blockCount;
			$rnaalignhash->{$qName}->{"blockSizes"}=$blockSizes;
			$rnaalignhash->{$qName}->{"qStarts"}=$qStarts;
			$rnaalignhash->{$qName}->{"tStarts"}=$tStarts;
		} elsif ($rnaalignhash->{$qName}->{matches} == $matches) {
			#print "$qName: tie score!\n";
			$rnaalignhash->{$qName}->{"tieScore"}++;
		}
	}
	close IN;

	# build @eb5/@eb3/@elength/@cmllength
	foreach my $hash (values %$rnaalignhash) {
		my @b5=split /,/,$hash->{"tStarts"};
		my @bsize=split /,/,$hash->{"blockSizes"};
		my @b3;
		for (my $i=0; $i<=$#b5; $i++) {
			$b3[$i]=$b5[$i]+$bsize[$i];
		}

		# @eb5 / @eb3: genomic positions of exon boundaries
		# index: 1 .. n
		$hash->{'eb5'}=\@b5;
		unshift @{$hash->{'eb5'}},-1;
		push @{$hash->{'eb5'}},-1;
		$hash->{'eb3'}=\@b3;
		unshift @{$hash->{'eb3'}},-1;
		push @{$hash->{'eb3'}},-1;

		# @elength: length of each exon (nt)
		# index: 0 .. (n+1)
		# 0 / n+1: unmapped portion
		$hash->{'elength'}=\@bsize;
		if ($hash->{'strand'} eq '+') {
			unshift @{$hash->{'elength'}},$hash->{'qStart'}; # unaligned nucleotides at 5'end
			push @{$hash->{'elength'}},$hash->{'qSize'}-$hash->{'qEnd'}; # unaligned nucleotides at 3'end
		} else {
			unshift @{$hash->{'elength'}},$hash->{'qSize'}-$hash->{'qEnd'};# unaligned nucleotides at 5'end (relative to genome positive strand)
			push @{$hash->{'elength'}},$hash->{'qStart'}; # unaligned nucleotides at 3'end (relative to genome positive strand)
		}

		# @cmllength: cumulative exon length
		# index: 0 .. (n+1)
		# 0 / n+1: unmapped portion
		$hash->{'cmllength'}->[0]=$hash->{'elength'}->[0];
		for (my $i=1; $i<scalar(@{$hash->{'elength'}}); $i++) {
			$hash->{'cmllength'}->[$i]=$hash->{'cmllength'}->[$i-1] + $hash->{'elength'}->[$i];
		}
	}

	return $rnaalignhash;
}

sub runBlat {
	shift @_;
	my ($blat,$database,$query,$output)=@_;
	$output||='output.psl';
	system(qq($blat $database $query -noHead $output >/dev/null 2>&1));
}

=head
sub pep2rnaPos {
	shift @_;
	my ($scanhash)=@_;

	
}
=cut
sub getRNA {
	shift @_;
	my ($scanhash,$rna_fas)=@_;
	my $rnaseqhash;

	# build %rna2out
	my %rna2out;
	foreach my $out (keys %$scanhash) {
		foreach my $rna (keys %{$scanhash->{$out}->{RNA}}) {
			$rna2out{$rna}{$out}='';
		}
	}

	# assumming the rna_fas is huge, the file will be read one row by another, rather than put into RAM at once
	open(IN,"$rna_fas") || die "cannot open ID.txt: $rna_fas\n";
	while (<IN>) {
		chomp;
		s/>//;
		#my $rna=$_;
		my $rna = (split /[\s\t]+/,$_)[0];

		my $seq=<IN>;
		chomp($seq);	

		if (defined($rna2out{$rna})) {
			foreach my $out (keys %{$rna2out{$rna}}) {
				# scanhash
				# RNA sequence
				$scanhash->{$out}->{RNA}->{$rna}->{seq}=$seq;
				# pep2rnaStartAbs
				my $rnalength=length($seq);
				$scanhash->{$out}->{RNA}->{$rna}->{rnalength}=$rnalength;
				if ($scanhash->{$out}->{RNA}->{$rna}->{strand} eq '+') {
					$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaStartAbs}=$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaStart};
					$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaEndAbs}=$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaEnd};
				} else {
					$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaStartAbs}=$rnalength-$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaEnd};
					$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaEndAbs}=$rnalength-$scanhash->{$out}->{RNA}->{$rna}->{pep2rnaStart};
				}
				# rnaseqhash
				$rnaseqhash->{$rna}=$seq;
			}
		}
	}
	close IN;

	return $rnaseqhash;
}

sub buildRNAsubHash {
	shift @_;
	my ($scanhash)=@_;

	# %scanhash{$out}{RNA}{$rnaID}{frame/ORFstart/ORFend/seq/pep2ORFstart/pep2ORFend/pep2rnaStart/pep2rnaEnd}
	foreach my $hash (values %$scanhash) {
		# AA2toAA13
		$hash->{pos} =~ /AA(\d+)toAA(\d+)$/;
		my $pep2ORFstart=$1 - 1; # 0-based
		my $pep2ORFend=$2;
		foreach my $item (keys %{$hash->{proteins}}) {
			# HWI-ST1200:152:C34NRACXX:5:1201:15343:14910/1|F0|21-78
			my ($rna,$frame,$ORFrange)=split /\|/,$item;
			my ($ORFstart,$ORFend)=split /\-/,$ORFrange;
			$hash->{RNA}->{$rna}->{frame}=$frame;
			$hash->{RNA}->{$rna}->{ORFstart}=$ORFstart;
			$hash->{RNA}->{$rna}->{ORFend}=$ORFend;
			$hash->{RNA}->{$rna}->{pep2ORFstart}=$pep2ORFstart;
			$hash->{RNA}->{$rna}->{pep2ORFend}=$pep2ORFend;

			# peptide to RNA
			$hash->{RNA}->{$rna}->{pep2rnaStart}=$ORFstart + 3 * $pep2ORFstart;
			$hash->{RNA}->{$rna}->{pep2rnaEnd}= $ORFstart + 3 * $pep2ORFend;

			# strand
			$hash->{RNA}->{$rna}->{strand}= ($frame =~ m/^F/)?'+':'-';
		}
	}
}

sub rmKnownPep {
	shift @_;
	my ($scanhash,$seqhash)=@_;
	my $novelScanHash;

	# for each scan in %scanhash: mark if ref peptides
	foreach my $outfile (keys %$scanhash) {
		my $nomod=$scanhash->{$outfile}->{nomod};
		foreach my $pro (keys %$seqhash) {
			if ($seqhash->{$pro} =~ m/$nomod/) {
				#$novelScanHash->{$outfile}=$scanhash->{$outfile};
				$scanhash->{$outfile}->{ref_protein}=1;
				last;
			}
		}
	}

	# assign terms to novelScanHash
	foreach my $outfile (keys %$scanhash) {
		if (defined($scanhash->{$outfile}->{ref_protein})
	                and $scanhash->{$outfile}->{ref_protein} == 1) {
		} else {
			$novelScanHash->{$outfile}=$scanhash->{$outfile};
		}
	}

	return $novelScanHash;
}
