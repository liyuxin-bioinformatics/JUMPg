#!/usr/bin/perl

## Release date: 05/01/2015
## Release version: version 12.0.0

package Utils::DatabaseUtils;
use strict;
use warnings;
use File::Basename;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new {
	my ($class, %arg)=@_;
    my $self = {};
    bless $self, $class;
	return $self;
}

sub generateFastaAndPit {
	## In this subroutine, $params is from builddb.params
	my ($self, $params, $dbName, $log) = @_;
	
	## Check input databases (fasta files)
	my %inputFastas;
	foreach my $key (keys %{$params}) {
		if ($key =~ /input_database/) {
			$key =~ /input_database(\d+)/;
			if (defined $1) {
				$inputFastas{$1} = $$params{$key};
			} else {
				print "Please specify the number of input database(s) from 1\n";
				print "e.g. input_database1 = ...\n";
				print "     input_database2 = ...\n";
				print $log "Please specify the number of input database(s) from 1\n";
				print $log "e.g. input_database1 = ...\n";
				print $log "     input_database2 = ...\n";
				exit;
			}
		}
	}
	
	## Generate output file names
	$dbName = basename($dbName);
	my $pitFile = $dbName.".pit";
	my $fastaFile = $dbName.".fasta";
	
	## Read protein abundances and annotations for PIT (according to input parameters)
	my @proteinAbundances;
	my %proteinAnnotations;
	foreach my $key (keys %{$params}) {
		if ($key =~ /list_protein_abundance/) {
			push (@proteinAbundances, $$params{$key})
		} elsif ($key =~ /list_/) {
			$proteinAnnotations{$key} = $$params{$key};
		}
	}
	my %abundances = readAbundance(\@proteinAbundances, $log);
	print "\n";
	print $log "\n";
	my %annotations = readAnnotation(\%proteinAnnotations, $log);
	print "\n";
	print $log "\n";

	## Generate @FASTAs array containing the names (and paths) of SOURCE .fasta files
	my @FASTAs;
	for (my $i = 1; $i <= scalar (keys %inputFastas); $i++) {
		push (@FASTAs, $inputFastas{$i});
	}
	
	## Whether or not include contaminants
	if ($$params{'include_contaminants'} == 1) {
		unshift (@FASTAs, $$params{'input_contaminants'});
	}
	
	## Read .fasta files (including contaminants if necessary) and
	## generate hashes	
	my @priorities;
	my (%protInfo, %protGroups, %pep2Group, %pits);
	for (my $i = 0; $i < scalar(@FASTAs); $i++) {
		my $fasta = $FASTAs[$i];
		my $priority = basename($fasta); #my $priority = $fasta;
		$priority =~ s/\.[^.]+$//; # $priority =~ s/\.\S+//;
		push (@priorities, $priority);
		($protInfo{$priority}, $protGroups{$priority}, $pep2Group{$priority}) =	pitHashGeneration($log, $fasta, \%abundances);
		($protInfo{$priority}, $protGroups{$priority}, $pep2Group{$priority}) = pitProteinGrouping($protInfo{$priority}, $protGroups{$priority}, $pep2Group{$priority}, \%abundances);
		@{$pits{$priority}} = pitGeneration($protGroups{$priority}); 
	}
	
	## Open an output .pit file and make a header line	
	open (PIT, ">", $pitFile) or die "Cannot open $pitFile\n";
	my $annotHeader = join("\t", sort (keys (%annotations)));
	print PIT "UniprotAC\tSJPGnumber\tGroupName\tAbundance\tResidueLength\tProteinName\tFullDescription\t".
					$annotHeader."\n";

	## Open an output .fasta file
	open (FASTA, ">", $fastaFile) or die "Cannot open $fastaFile\n";
	my %protHash;	# Protein hash for .fasta file (necessary for generating decoys)
	
	## Write protein (and contaminants, if necessary) information to the output files
	my $nGroups = 0;
	for (my $i = 0; $i < scalar(@priorities); $i++) {
		## Write to the output .pit file
		if ($$params{'include_contaminants'} eq 1) {
			if ($i == 0) {
				$nGroups = 0;
			} elsif ($i == 1) {
				$nGroups = 100000;
			}
		} elsif ($$params{'include_contaminants'} eq 0 && $i == 0) {
			$nGroups = 100000;
		}
		my @pit = @{$pits{$priorities[$i]}};
		for (my $j = 0; $j < scalar(@pit); $j++) {
			for (my $k = 0; $k < scalar(@{$pit[$j]{'PROTEINS'}}); $k++) {
				my $protein = $pit[$j]{'PROTEINS'}[$k];
				my $groupNumber = sprintf("SJPG%06d", ($j + 1 + $nGroups)).".".sprintf("%03d", ($k + 1));
				my $groupName = $pit[$j]{'NAME'};
				my $abundance = $protInfo{$priorities[$i]}{$protein}{'ABUN'};
				my $length = $protInfo{$priorities[$i]}{$protein}{'LEN'};
				my @annotValues;
				foreach my $key (sort (keys %annotations)) {
					if (defined $annotations{$key}{$protein}) {
						push (@annotValues, "O");
					} else {
						push (@annotValues, "-");
					}
				}
				my $proteinAnnot = join("\t", @annotValues);
				my $header = $protInfo{$priorities[$i]}{$protein}{'HEAD'};
				my ($proteinName, $fullDescription) = ($header =~ /^\>([a-z]{2}\|\S+\|\S+)\s(.+)$/);
				if ($priorities[$i] eq "contaminants") {
					$proteinName =~ s/Contaminant/Con/;
				}
				print PIT "$protein\t$groupNumber\t$groupName\t$abundance\t$length\t$proteinName\t$fullDescription\t".
							 	$proteinAnnot."\n";
			}
		}
		$nGroups += scalar(@pit);
		
		## Write to the output .fasta file
		foreach my $key (sort keys %{$protInfo{$priorities[$i]}}) {
			my $header = $protInfo{$priorities[$i]}{$key}{'HEAD'};
			my $seq = $protInfo{$priorities[$i]}{$key}{'SEQ'};
			$protHash{$priorities[$i]}{$header} = $seq;
			print FASTA "$header\n$seq\n";
		}
	}
	
	## Write decoy information to the output files, if necessary
	print "\n";
	print $log "\n";
	if ($$params{'decoy_generation'} == 1) {
		my $nDecoys = 0;
		for (my $i = 0; $i < scalar(@priorities); $i++) {

			## Write to the output .pit file
			my @pit = @{$pits{$priorities[$i]}};
			for (my $j = 0; $j < scalar(@pit); $j++) {
				for (my $k = 0; $k < scalar(@{$pit[$j]{'PROTEINS'}}); $k++) {
					my $protein = $pit[$j]{'PROTEINS'}[$k];
					my $groupNumber = sprintf("SJPG%06d", ($j + 1 + $nGroups)).".".sprintf("%03d", ($k + 1));
					my ($groupName, $abundance) = ("-", "-");
					my $length = $protInfo{$priorities[$i]}{$protein}{'LEN'};
					my @decoyAnnotValues;
					foreach my $key (keys %annotations) {
						push (@decoyAnnotValues, "-");
					}
					my $decoyAnnot = join("\t", @decoyAnnotValues);
					my $header = $protInfo{$priorities[$i]}{$protein}{'HEAD'};
					my ($proteinName, $fullDescription) = ($header =~ /^\>([a-z]{2}\|\S+\|\S+)\s(.+)$/);
					if ($priorities[$i] eq "contaminants") {
						$proteinName =~ s/Contaminant/Con/;
					}
					$protein = "##Decoy__".$protein;
					$proteinName = "##Decoy__".$proteinName;
					print PIT "$protein\t$groupNumber\t$groupName\t$abundance\t$length\t$proteinName\t$fullDescription\t".
									$decoyAnnot."\n";
				}
			}
			$nGroups += scalar(@pit);
			
			## Write to the output .fasta file
			my %decoys = generateDecoys($protHash{$priorities[$i]}, $$params{'decoy_generation_method'}, $log);
			foreach my $key (sort keys %decoys) {
				print FASTA "$key\n$decoys{$key}\n";
				$nDecoys++;
				print "Generating $nDecoys decoys\r";
			}
		}
		print $log "Generating $nDecoys decoys\n";
	} else {
		print "Decoys are not considered\n";
		print $log "Decoys are not considered\n";
	}
	close (FASTA);
	close (PIT);
	print "\n";
	print "$fastaFile has been successfully generated\n";
	print "$pitFile has been successfully generated\n";
	print $log "\n";
	print $log "$fastaFile has been successfully generated\n";
	print $log "$pitFile has been successfully generated\n";	
	
}

sub pitHashGeneration {
	my ($log, $fasta, $uniprot2abundance) = @_;
	my %uniprot2abundance = %{$uniprot2abundance};
	
	print "Loading $fasta (may take a while)\n";
	print $log "Loading $fasta (may take a while)\n";
	
	####################
	## Initialization ##
	####################
	my ($dbInfo, $protAC, $geneName, $isoNo, $seq);
	my %protInfo;
	## Keys used in %protInfo
	#### First key: UniprotAC
	#### Second keys
	###### HEAD: header information as in .fasta file
	###### DB: database information - the first two characters in .fasta file (e.g. sp, tr, etc)
	###### GN: Gene Name in the header of fasta entry
	###### ISO: isoform number (e.g. for A12345-5, isoform number is 5)
	###### SEQ: amino acid sequence in the fasta entry
	###### LEN: length of amino acid sequence
	###### PEP: fully tryptic-digested peptides of the sequence
	###### ABUN: abundance level as in the abundance file derived from our in-house AD proteome data with emPAI method
	###### GROUP: group name (i.e. gene name) to which a protein belong (not grouped, then undef)
	
	my %pep2Group;
	## Keys used in %pep2Group
	#### key: fully tryptic peptide generated by 'digest' subroutine
	#### value: Gene name (GN) from which the key peptide is derived
	
	my %protGroups;
	## Keys used in %protGroups
	#### First key: Gene name (GN) equivalent to a group name
	#### Second keys
	###### ABUN: group abundance derived from protein abundances in the group
	###### PROTEINS: proteins (UniprotACs) included in the group
	
	my $nProteinsTot = 0;
	my $nProteinsGn = 0;
	my $nProteinsGnAbundance = 0;
	my $nProteinsNoGnAbundance = 0;
	
	############################################
	## FASTA parsing and generation of hashes ##
	############################################
	open (FASTA, "<", $fasta) or die "Cannot open $fasta file\n";
	while (<FASTA>) {
		chomp;
		## Header information
		if ($_ =~ /^\>/) {
			$nProteinsTot++;
			## Generation of fully tryptic peptides
			if (defined $protAC && defined $protInfo{$protAC}{'SEQ'}) {
				$protInfo{$protAC}{'LEN'} = length($protInfo{$protAC}{'SEQ'});
				@{$protInfo{$protAC}{'PEP'}} = digest($protInfo{$protAC}{'SEQ'});
				if (defined $protInfo{$protAC}{'GROUP'}) {
					foreach my $peptide (@{$protInfo{$protAC}{'PEP'}}) {
						$pep2Group{$peptide} = $geneName;
						if (defined $pep2Group{$peptide} && $pep2Group{$peptide} ne $geneName) {
							die "\nPeptide $peptide is mapped to multiple protein groups\n";
						}
					}
				}
			}
			## Must use original .fasta file(s) including protein entries ONLY
			## Neither contaminants nor decoys are allowed 
			## $dbInfo and $protAC should be defined
			($dbInfo, $protAC) = $_ =~ /^\>([a-z]{2})\|(\S+)\|/;
			if (!defined $dbInfo || !defined $protAC) {
				die "\nCheck $fasta file for the entry $_\n";
			}
			## $geneName may be "undef" (when there's no GN field in the header)
			if ($_ =~ /PE=/) {
				($geneName) = $_ =~ /\sGN=(.+)(?=\sPE)/;
			} else {
				($geneName) = $_ =~ /\sGN=(.+)(?=$)/;
			}
			if (defined $geneName) {
				$nProteinsGn++;
				## Pre-grouping of proteins based on gene name
				push (@{$protGroups{$geneName}{'PROTEINS'}}, $protAC);
				$protInfo{$protAC}{'GROUP'} = $geneName;
			} else {
				$protInfo{$protAC}{'GROUP'} = undef;
			}
			## $isoNo may be "undef" (when there's no isoform associated with $protAC)
			($isoNo) = $protAC =~ /-(\d)/;
			## Put the information of a protein to %protInfo
			my $header = $_;
			if (($header =~ /^\>co/) && ($header =~ /\|Contaminant/)) {
				$header =~ s/\|Contaminant/\|Con/;
				$protInfo{$protAC}{'HEAD'} = $header;
			} else {
				$protInfo{$protAC}{'HEAD'} = $header;
			}
			$protInfo{$protAC}{'DB'} = $dbInfo;
			$protInfo{$protAC}{'GN'} = $geneName;
			$protInfo{$protAC}{'ISO'} = $isoNo;
			$protInfo{$protAC}{'SEQ'} = "";
			## Obtain abundance information from our in-house AD proteome data with emPAI method
			if (defined $uniprot2abundance{$protAC}) {
				$protInfo{$protAC}{'ABUN'} = $uniprot2abundance{$protAC};
				if (defined $geneName) {
					$nProteinsGnAbundance++;
				} else {
					$nProteinsNoGnAbundance++;
				}
			} else {
				$protInfo{$protAC}{'ABUN'} = "-";
			}
		## Sequence information
		} elsif (defined $protAC) {
			## Concatenate multi-line sequences to one sequence for a protein
			$protInfo{$protAC}{'SEQ'} = $protInfo{$protAC}{'SEQ'}.$_;
		}
	}
	## Dealing with the last entry
	if (defined $protAC && defined $protInfo{$protAC}{'SEQ'}) {
		$protInfo{$protAC}{'LEN'} = length($protInfo{$protAC}{'SEQ'});
		@{$protInfo{$protAC}{'PEP'}} = digest($protInfo{$protAC}{'SEQ'});
		if (defined $protInfo{$protAC}{'GROUP'}) {
			foreach my $peptide (@{$protInfo{$protAC}{'PEP'}}) {
				$pep2Group{$peptide} = $geneName;
				if (defined $pep2Group{$peptide} && $pep2Group{$peptide} ne $geneName) {
					die "Peptide $peptide is mapped to multiple protein groups\n";
				}
			}
		}
	}
	close (FASTA);
	
	## Stats for proteins
	my $nProteinsGnNoAbundance = $nProteinsGn - $nProteinsGnAbundance;
	my $nProteinsNoGn = $nProteinsTot - $nProteinsGn;
	my $nProteinsNoGnNoAbundance = $nProteinsNoGn - $nProteinsNoGnAbundance;	
	print "  Total $nProteinsTot proteins\n";
	print "    $nProteinsGn proteins with gene names (GN) used for protein grouping\n";
	print "      $nProteinsGnAbundance proteins with abundance information\n";
	print "      $nProteinsGnNoAbundance proteins without abundance information\n";
	print "    $nProteinsNoGn proteins without GNs\n";
	print "      $nProteinsNoGnAbundance proteins with abundance information\n";
	print "      $nProteinsNoGnNoAbundance proteins without abundance information\n";
	print $log "  Total $nProteinsTot proteins\n";
	print $log "    $nProteinsGn proteins with gene names (GN) used for protein grouping\n";
	print $log "      $nProteinsGnAbundance proteins with abundance information\n";
	print $log "      $nProteinsGnNoAbundance proteins without abundance information\n";
	print $log "    $nProteinsNoGn proteins without GNs\n";
	print $log "      $nProteinsNoGnAbundance proteins with abundance information\n";
	print $log "      $nProteinsNoGnNoAbundance proteins without abundance information\n";
	return (\%protInfo, \%protGroups, \%pep2Group);
}

sub pitProteinGrouping {
	my ($protInfo, $protGroups, $pep2Group, $uniprot2abundance) = @_;
	my %protInfo = %{$protInfo};
	my %protGroups = %{$protGroups};
	my %pep2Group = %{$pep2Group};
	my %uniprot2abundance = %{$uniprot2abundance};
	
	#########################################################
	## Form protein groups using GN and update information ##
	#########################################################
	## For each protein group formed based on GN
	my @dbPriority = ("sp", "tr", "co", "cu");
	my @groupAbundances;
	my @groupKeys;
	my $nMeasuredGroups = 0;
	my $nMeasuredProteins = 0;
	foreach my $GN (keys %protGroups) {
		push (@groupKeys, $GN);
		my @proteins = @{$protGroups{$GN}{'PROTEINS'}};	## proteins in a group
		my $nProteins = scalar(@proteins);
		my @abundances;	
		
		## Generate hash %proteinsWithPriority
		## Keys: database names (e.g. sp, tr, etc.)
		## Values: proteins (in the group) belong to the database
		my %proteinsWithPriority;
		my $count = 0;
		for (my $i = 0; $i < scalar(@dbPriority); $i++) {
			foreach my $protein (@proteins) {
				if ($protInfo{$protein}{'DB'} eq $dbPriority[$i]) {
					$count++;
					push (@{$proteinsWithPriority{$dbPriority[$i]}}, $protein);
				}			
				if (defined $uniprot2abundance{$protein}) {
					push (@abundances, $uniprot2abundance{$protein});
					$nMeasuredProteins++;
				}
			}
		}
		
		## Define the group abundance
		## Assign the largest abundance of proteins in the group to the group abundance
		@abundances = uniq(@abundances);
		if (scalar(@abundances) > 0) {
			@abundances = sort {$b <=> $a} @abundances;
			$protGroups{$GN}{'ABUN'} = shift (@abundances);
			$nMeasuredGroups++;	## number of groups which have representiative abundances
		} else {
			$protGroups{$GN}{'ABUN'} = 0;
		}
		push (@groupAbundances, $protGroups{$GN}{'ABUN'});
		
		## Check
		if ($nProteins ne $count) {
			die "Some proteins in group $GN are not considered\n";
		}
		
		## For the proteins in each prioritized database, for example, proteins in SWISSPROT (sp),
		#### 1. Count the number of isoforms in a protein family
		#### 2. Put a representative (canonical or less isoform numbered) protein on top of the list
		#### 2.1. If there are multiple proteins which have the same number of isoforms,
		####      the one with longer length of residues comes first
		#### 3. Sort isoforms of the representative protein in alphanumerical order and
		####    put them after the representative one
		####
		#### Example
		#### 1. If there are proteins in a group as following,
		####	A12345, A12345-2, A12345-4, B45678, B45678-1, C67890
		####	then, put A12345 on the top of list and its isoforms after it
		####	The rationale is that a protein with more isoforms is likely to be a major component in the group
		#### 2. If there are proteins as following,
		####	A12345, A12345-1, B45678, B45678-2
		####	(A12345 and B45678 have the same number of isoforms)
		####	then, compare the length of residues between A12345 and B45678, and put a longer one
		####	on the top as a representative group protein 
		my @sortedProteins;
		for (my $i = 0; $i < scalar(@dbPriority); $i++) {
			my $db = $dbPriority[$i];
			if (!defined $proteinsWithPriority{$db}) {
				next;
			}
			my %isoHash;
			foreach my $protein (@{$proteinsWithPriority{$db}}) {
				my $key = $protein;
				$key =~ s/\-\d+//;
				$isoHash{$key}{'FREQ'}++;
				push(@{$isoHash{$key}{'MEM'}}, $protein);
			}
			foreach my $key (keys %isoHash) {
				## Sort member proteins of a hash in alphanumerical ascending order
				#### Subroutine "isoformSort" is used to avoid the case that
				#### P12345-10 comes before P12345-9 (when "cmp" is used) 
				@{$isoHash{$key}{'MEM'}} = sort isoformSort @{$isoHash{$key}{'MEM'}};
				## For each key of %isoHash, find the most canonical-like protein and its length
				$isoHash{$key}{'CANONICAL'} = $isoHash{$key}{'MEM'}[0];
				$isoHash{$key}{'LEN'} = $protInfo{$isoHash{$key}{'CANONICAL'}}{'LEN'}; 
			}
			
			## If there are multiple protein families,
			## 1. sort according to the number of isoforms
			## 2. sort according to the length of residues
			my @sortedIsoHashKey;
			if (scalar(keys %isoHash) > 1) {
				@sortedIsoHashKey = sort {$isoHash{$b}{'FREQ'} <=> $isoHash{$a}{'FREQ'} ||
					$isoHash{$b}{'LEN'} <=> $isoHash{$a}{'LEN'} ||
					$isoHash{$a} cmp $isoHash{$b}} keys %isoHash;
				for (my $j = 0; $j < scalar(@sortedIsoHashKey); $j++) {
					push(@sortedProteins, @{$isoHash{$sortedIsoHashKey[$j]}{'MEM'}});
				}
			## If there's only one protein (family)
			} else {
				my $key = (keys %isoHash)[0];
				push (@sortedProteins, @{$isoHash{$key}{'MEM'}});
			}
		}
		@{$protGroups{$GN}{'PROTEINS'}} = @sortedProteins;
	}	
	
	###############################################################
	## Sort protein groups according to either abundance or name ##
	###############################################################
	
	#### Sort them according to "group abundance" in descending order
	#### For the groups without "group abundance", sort them in alphanumerical ascending order
	
	#### 1. Obtain the sorted order (@groupOrder) based on @groupAbundances 
	#### 2. Sort @groupKeys according to the @groupOrder
	####    The first element of @groupKeys is the key (name) of protein group with the largest abundance
	#### 3. Divide @groupKeys into @measuredGroupKeys (groups with abundance) and 
	####    @unmeasuredGroupKeys (groups without abundance)
	####    Since @groupAbundances have real abundance values as well as nominal abundance values (=0)
	####    for the groups without abundance information, 
	####    some keys (names) of protein groups actually sorted by abundance and others should be separated
	#### 4. Sort @unmeasuredGroupKeys in alphanumerical ascending order  
	my @groupOrder = sort {$groupAbundances[$b] <=> $groupAbundances[$a]} 0..$#groupAbundances;
	@groupKeys = @groupKeys[@groupOrder];
	my @unmeasuredGroupKeys = @groupKeys[$nMeasuredGroups..(scalar(@groupKeys) - 1)];
	@unmeasuredGroupKeys = sort {$a cmp $b} @unmeasuredGroupKeys;
	my @measuredGroupKeys = @groupKeys[0..($nMeasuredGroups - 1)];
	@groupKeys = (@measuredGroupKeys, @unmeasuredGroupKeys);
	my $nGroups = scalar(@groupKeys);
	
	## Assign orders to protein groups based on the sorted @groupKeys
	my @sortedProtGroups;
	for (my $i = 0; $i < $nGroups; $i++) {
		$protGroups{$groupKeys[$i]}{'ORDER'} = $i;
		$sortedProtGroups[$i]{'GN'} = $groupKeys[$i];
		$sortedProtGroups[$i]{'ABUN'} = $protGroups{$groupKeys[$i]}{'ABUN'};
		$sortedProtGroups[$i]{'PROTEINS'} = $protGroups{$groupKeys[$i]}{'PROTEINS'};
	}
	
	############################################
	## Handling of proteins not grouped by GN ##
	############################################
	## Assign the 'ungrouped' proteins to each prioritized database (i.e. sp, tr, etc.)
	my %ungroupedProteins;
	foreach my $protein (keys %protInfo) {
		if (!defined $protInfo{$protein}{'GROUP'}) {	## i.e. not grouped by GN yet
			for (my $i = 0; $i < scalar(@dbPriority); $i++) {
				if ($protInfo{$protein}{'DB'} eq $dbPriority[$i]) {
					push (@{$ungroupedProteins{$dbPriority[$i]}}, $protein);
				}
			}
		}
	}
	
	## For the ungrouped proteins in each prioritized database
	## 1. Sort them according to abundance (if possible) in descending order
	## 2. Sort the remainder (i.e. ungrouped proteins without abundance information)
	##    in alphanumerical ascending order
	my @sortedUngroupedProteins;
	for (my $i = 0; $i < scalar(@dbPriority); $i++) {
		my @measured;
		my @abundances;
		my @unmeasured;
		foreach my $protein (@{$ungroupedProteins{$dbPriority[$i]}}) {
			if ($protInfo{$protein}{'ABUN'} ne "-") {
				push (@measured, $protein);
				push (@abundances, $protInfo{$protein}{'ABUN'});
			} else {
				push (@unmeasured, $protein);
			}
		}
		## Sort @measured according to @abundances in descending order
		if (scalar(@measured) > 0) {
			my @order = sort {$abundances[$b] <=> $abundances[$a]} 0..$#abundances;
			@measured = @measured[@order];
		}		
		## Sort @unmeasured in alphanumeric order
		## or sort @unmeasured according to the length of residues ??? ##
		if (scalar(@unmeasured) > 0) {
			@unmeasured = sort {$a cmp $b} @unmeasured;
		}
		push (@sortedUngroupedProteins, @measured, @unmeasured);
	}
	
	## For each ungrouped protein (in @sortedUngroupedProteins) in each prioritized database,
	## 1. Compare tryptic peptides of the protein with protein in each group
	## 2. If there are shared peptides, 
	##    sort the groups according to the number of shared peptides and then according to the group order
	##    Add the protein into the bottom of the group and assign a subgroup number
	## 3. If not, generate a new group composed of the protein
	for (my $i = 0; $i < scalar(@sortedUngroupedProteins); $i++) {
		my $protein = $sortedUngroupedProteins[$i];
		my %sharedGroups;
		foreach my $peptide (@{$protInfo{$protein}{'PEP'}}) {
			if (defined $pep2Group{$peptide}) {
				my $GN = $pep2Group{$peptide};
				$sharedGroups{$GN}{'numShared'}++;
				$sharedGroups{$GN}{'groupOrder'} = $protGroups{$GN}{'ORDER'};
			}
		}
		## If there are shared peptides with existing groups,
		## assign the protein to one of the groups
		if (scalar(keys %sharedGroups) > 0) {
			my $numShared = 0;
			my @sharedGroupNames;
			my @sharedGroupOrders;
			foreach my $GN (keys %sharedGroups) {
				if ($sharedGroups{$GN}{'numShared'} >= $numShared) {
					push (@sharedGroupNames, $GN);
					push (@sharedGroupOrders, $sharedGroups{$GN}{'groupOrder'});
					$numShared = $sharedGroups{$GN}{'numShared'};
				}
			}
			## If there are multiple groups sharing the same number of peptides with the protein,
			## choose a group with a higher order (rank)
			my $sharedGroupName;
			if (scalar(@sharedGroupNames) > 1) {
				my @order = sort {$sharedGroupOrders[$a] <=> $sharedGroupOrders[$b]} 0..$#sharedGroupOrders;
				$sharedGroupName = $sharedGroupNames[shift (@order)];
			} else {
				$sharedGroupName = shift (@sharedGroupNames);
			}
			## Add the protein to the group corresponding to $sharedGroupName
			push (@{$protGroups{$sharedGroupName}{'PROTEINS'}}, $protein);
	
		## If there's no shared peptides with existing groups,
		## generate a new protein group composed of the protein
		} else {
			## Generation of a new protein group
			$nGroups++;
			my $newGroupName;
			if (defined $protInfo{$protein}{'GN'}) {
				$newGroupName = $protInfo{$protein}{'GN'};
			} else {
				$newGroupName = $protein;
			}
			$protGroups{$newGroupName}{'ABUN'} = $protInfo{$protein}{'ABUN'};
			$protGroups{$newGroupName}{'ORDER'} = $nGroups - 1;
			push (@{$protGroups{$newGroupName}{'PROTEINS'}}, $protein);
			
			## Update pep2Group
			foreach my $peptide (@{$protInfo{$protein}{'PEP'}}) {
				$pep2Group{$peptide} = $newGroupName;
			}
		}
	}
	return (\%protInfo, \%protGroups, \%pep2Group);
}

sub readAbundance {
	## Protein abundance: Col1 = UniprotAC (e.g. P12345), Col2 = Abundance, Col3 = additional information, ...
	## The file should be a tab-delimited text format	
	my ($proteinAbundances, $log) = @_;
	my %abundances;
	if (scalar (@{$proteinAbundances}) > 0) {
		print "Loading protein abundance information\n";
		print $log "Loading protein abundance information\n";
		foreach my $proteinAbundanceFile (@{$proteinAbundances}) {
			my $nAbundances = 0;
			open (INFILE, "<", $proteinAbundanceFile) or die "Cannot open $proteinAbundanceFile\n";
			while (<INFILE>) {
				chomp;
				my @elems = split(/\t/, $_);
				next if ($elems[0] !~ /[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/);
				$abundances{$elems[0]} = $elems[1];
				$nAbundances++;
			}
			close (INFILE);
			print "  $nAbundances protein abundance information from $proteinAbundanceFile\n";
			print $log "  $nAbundances protein abundance information from $proteinAbundanceFile\n";
		}	
	} else {
		print "There is no protein abundance information\n";
		print $log "There is no protein abundance information\n";
	}
	return (%abundances);
}

sub readAnnotation {
	## Protein annotation: Col1 = UniprotAC (e.g. P12345), Col2 = name of annotation (e.g. TF motif name, kinase name, etc.)
	## The file should be a tab-delimited text format
	my ($proteinAnnotations, $log) = @_;
	my %annotations;
	if (scalar (keys %{$proteinAnnotations}) > 0) {
		print "Loading protein annotation information\n";
		print $log "Loading protein annotation information\n";
		foreach my $key (sort (keys %{$proteinAnnotations})) {
			my $annotationName = $key;
			$annotationName =~ s/list\_//;
			my $nAnnotations = 0;			
			open (INFILE, "<", $$proteinAnnotations{$key}) or die "Cannot open $$proteinAnnotations{$key}\n";
			while (<INFILE>) {
				chomp;
				my @elems = split(/\t/, $_);
				next if ($elems[0] !~ /[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/);
				if (!defined $elems[1]) {
					print "  No annotation information for protein $elems[0] in $$proteinAnnotations{$key}\n";
					print "  Please check the input annotation file\n";
					print $log "  No annotation information for protein $elems[0] in $$proteinAnnotations{$key}\n";
					print $log "  Please check the input annotation file\n";
					exit;
				} else {
					$annotations{$annotationName}{$elems[0]} = $elems[1];
					$nAnnotations++;
				}				
			}
			close (INFILE);
			print "  $nAnnotations $annotationName annotation information from $$proteinAnnotations{$key}\n";
			print $log "  $nAnnotations $annotationName annotation information from $$proteinAnnotations{$key}\n";
		}		
	} else {
		print "There is no protein annotation information\n";
		print $log "There is no protein annotation information\n";
	}
	return (%annotations);
}

sub pitGeneration {
	#################################################
	## Generation of protein inference table (PIT) ##
	#################################################
	my ($protGroups) = @_;
	my @pit;
	foreach my $key (keys %{$protGroups}) {
		my $order = $$protGroups{$key}{'ORDER'};
		$pit[$order]{'NAME'} = $key;
		@{$pit[$order]{'PROTEINS'}} = @{$$protGroups{$key}{'PROTEINS'}}; 
	}
	return (@pit);
}

sub generateDecoys {
	my ($proteins, $option, $log) = @_;
	my %decoys;	
	foreach my $key (keys %{$proteins}) {
		my $decoyKey = $key;
		$decoyKey =~ s/>(\w+)/>##Decoy__$1/;
		if (!defined $decoyKey) {
			print "There's a problem in the header of a protein entry, $key\n";
			print $log "There's a problem in the header of a protein entry, $key\n";
			exit;
		} else {
			my $decoySequence;
			if ($option == 1) {
				$decoySequence = reverse($$proteins{$key});
			} elsif ($option == 2) {
				$decoySequence = reverse($$proteins{$key});
				$decoySequence =~ s/([^KR]{1})([KR]+)/$2$1/g;
			} else {
				print "Please set \"decoy_generation_method\" parameter correctly\n";
				print "It should be either 1 (reverse) or 2 (reverse and then swap every K/R with the preceding AA)\n";
				print $log "Please set \"decoy_generation_method\" parameter correctly\n";
				print $log "It should be either 1 (reverse) or 2 (reverse and then swap every K/R with the preceding AA)\n";
				exit;
			}
			if (defined $decoySequence) {
				$decoys{$decoyKey} = $decoySequence;
			} else {
				print "There's a problem in generating the decoy of a protein entry, $key\n";
				print $log "There's a problem in generating the decoy of a protein entry, $key\n";
				exit;
			}
		}
	}
	return (%decoys);
}

sub digest {
	##################################################################
	## Generation of fully tryptic peptides (>= 7AAs and <= 25AAs)	##
	##################################################################
	my $seq = shift (@_);
	my @allPeptides = split (/(?<=[KR])/, $seq);
	my @selPeptides;
	foreach my $peptide (@allPeptides) {
		if (length($peptide) >= 7 && length($peptide) <= 25) {
			push (@selPeptides, $peptide);
		}
	}
	return (@selPeptides);
}

sub isoformSort {
	##############################################################
	## When 'sort' function is used for P12345-10 and P12345-8, ##
	## P12345-10 comes before P12345-8                          ##
	## This subroutine prevents such cases                      ##
	##############################################################
	my ($numberA) = $a =~ /-(\d+)/;
	if (!defined $numberA) {
		$numberA = 0;
	}
	my ($numberB) = $b =~ /-(\d+)/;
	if (!defined $numberB) {
		$numberB = 0;
	}
	return $numberA <=> $numberB;
}

sub uniq (@) {
    my %seen = ();
    grep { not $seen{$_}++ } @_;
}

sub generateDbName {
	## In this subroutine, $params is from jump.params	
	my ($self, $params, $log) = @_;
	print "Conditions for a new database\n";
	print $log "Conditions for a new database\n";
	
	##################################################
	## Automatic generation of a new database name	##
	##################################################
	
	my @dbSuffix;
	## Enzyme information
	if ($$params{'enzyme_info'} =~ /^Tryptic KR P/ && $$params{'digestion'} eq "full") {
		push (@dbSuffix, "ft");
		print "  Enzyme: fully tryptic\n";
		print $log "  Enzyme: fully tryptic\n";
	} elsif ($$params{'enzyme_info'} =~ /^Tryptic KR P/ && $$params{'digestion'} eq "partial") {
		push (@dbSuffix, "pt");
		print "  Enzyme: partial tryptic\n";
		print $log "  Enzyme: partial tryptic\n";
	} elsif ($$params{'enzyme_info'} =~ /^Chymotrypsin/) {
		push (@dbSuffix, "chymo");
		print "  Enzyme: chymotrypsin\n";
		print $log "  Enzyme: chymotrypsin\n";
	} elsif ($$params{'enzyme_info'} =~ /^No_Enzyme/) {
		push (@dbSuffix, "none");
		print "  Enzyme: none\n";
		print $log "  Enzyme: none\n";
	} elsif ($$params{'enzyme_info'} =~ /^Lys_C/) {
		push (@dbSuffix, "lycs");
		print "  Enzyme: Lys_C\n";
		print $log "  Enzyme: Lys_C\n";
	} elsif ($$params{'enzyme_info'} =~ /^Arg_C/) {
		push (@dbSuffix, "argc");
		print "  Enzyme: Arg_C\n";
		print $log "  Enzyme: Arg_C\n";
	} elsif ($$params{'enzyme_info'} =~ /^Glu_C/)	{
		push (@dbSuffix, "gluc");
		print "  Enzyme: Glu_C\n";
		print $log "  Enzyme: Glu_C\n";
	} elsif ($$params{'enzyme_info'} =~ /^Asp_N/) {
		push (@dbSuffix, "aspn");
		print "  Enzyme: Asp_N\n";
		print $log "  Enzyme: Asp_N\n";
	} elsif ($$params{'enzyme_info'} =~ /^Cyanogen_Bromide/) {
		push (@dbSuffix, "cnbr");
		print "  Enzyme: cyanogen_bromide\n";
		print $log "  Enzyme: cyanogen_bromide\n";
	} else {
		print "Please set a correct enzyme information in a JUMP parameter file\n";
		print $log "Please set a correct enzyme information in a JUMP parameter file\n";
		exit;
	}
	## Miscleavage
	if (defined $$params{'max_mis_cleavage'}) {
		my $mc = "mc".$$params{'max_mis_cleavage'};
		push (@dbSuffix, $mc);
		print "  Maxmimum miscleavage: $$params{'max_mis_cleavage'}\n";
		print $log "  Maxmimum miscleavage: $$params{'max_mis_cleavage'}\n";
	} else {
		print "Please set \"max_mis_cleavage\" parameter in a JUMP parameter file\n";
		print $log "Please set \"max_mis_cleavage\" parameter in a JUMP parameter file\n";
		exit;
	}	
	## Cysteine modification
	if ($$params{'add_C_Cysteine'} == 0) {
		push (@dbSuffix, "c0");
		print "  Static modification: Cysteine 0\n";
		print $log "  Static modification: Cysteine 0\n";
	} elsif (int($$params{'add_C_Cysteine'}) == 57) {
		push (@dbSuffix, "c57");
		print "  Static modification: C_Cysteine $$params{'add_C_Cysteine'}\n";
		print $log "  Static modification: C_Cysteine $$params{'add_C_Cysteine'}\n";
	}
	## Other static modifications
	foreach my $key (keys %{$params}) {
		if ($key =~ /^add_([A-BD-JL-Z])_([a-zA-Z0-9\_]+)/ && int($$params{$key}) > 0) {	## Do not consider "K" here
			my $mod = lc($1).int($$params{$key});
			push (@dbSuffix, $mod);
			print "  Static modification: $1\_$2 $$params{$key}\n";
			print $log "  Static modification: $1\_$2 $$params{$key}\n";
		}
	}
	## Phosphorylation
	if (lc($$params{'search_engine'}) eq "jump") {
		my $phosphoAAs;
		foreach my $key (sort (keys %{$params})) {
			if ($key =~ /dynamic_([A-Z])/) {
				my $dynamicAA = $1;
				if ($dynamicAA =~ /[STY]/) {
					$phosphoAAs = $phosphoAAs.$dynamicAA;
				}
				print "  Dynamic modification: $dynamicAA $$params{$key}\n";
				print $log "  Dynamic modification: $dynamicAA $$params{$key}\n";
			}
		}
		if (defined $phosphoAAs) {
			if (length($phosphoAAs) == 3 && int($$params{'dynamic_S'}) == 79 &&
				$$params{'dynamic_S'} == $$params{'dynamic_T'} && $$params{'dynamic_T'} == $$params{'dynamic_Y'}) {
				push (@dbSuffix, "pho");					
			} else {
				print "For phosphorylation, dynamic modifications for STY should be defined and they should be identical\n";
				print $log "For phosphorylation, dynamic modifications for STY should be defined and they should be identical\n";
				exit;				
			}
		}
	}
	## TMT modification
	if (int($$params{'add_Nterm_peptide'}) == 229 && int($$params{'add_K_Lysine'}) == 229) {
		push (@dbSuffix, "TMT_K229");
		print "  N-terminal and Lysine modifications for TMT\n";
		print "  N-terminal peptide modification: $$params{'add_Nterm_peptide'}\n";
		print "  K_Lysine modification: $$params{'add_K_Lysine'}\n";
		print $log "  N-terminal and Lysine modifications for TMT\n";
		print $log "  N-terminal peptide modification: $$params{'add_Nterm_peptide'}\n";
		print $log "  K_Lysine modification: $$params{'add_K_Lysine'}\n";
	} elsif (int($$params{'add_Nterm_peptide'}) == 0 && int($$params{'add_K_Lysine'}) == 0) {
		print "  No TMT modification\n";
		print $log "  No TMT modification\n";
=head
	} elsif (int($$params{'add_Nterm_peptide'}) != int($$params{'add_K_Lysine'})) {
		print "Please set \"add_Nterm_peptide\" and \"add_K_Lysine\" paremeters correctly in a JUMP parameter file for TMT\n";
		print "Two modifications should be the same\n";
		print $log "Please set \"add_Nterm_peptide\" and \"add_K_Lysine\" paremeters correctly in a JUMP parameter file for TMT\n";
		print $log "Two modifications should be the same\n";
		exit;
=cut
	}
	my $suffix = join("_", @dbSuffix);
	return ($suffix);
}

1;
