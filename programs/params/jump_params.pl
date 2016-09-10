#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use File::Basename;

## Directory information
my $currDir = getcwd()."/ParameterFiles";
if (!-e $currDir) {
	system ("mkdir $currDir");
}

## Version information
my $version = "12.0.0";
my $releaseDate = "05/01/2015";

## Search engines/Conditions
my @searchEngines = ("JUMP", "SEQUEST");
my @conditions  = ("HH", "HL", "TMThh", "TMThhpho");

##################################################
## Search parameter file generation (jump -s)	##
##################################################
foreach my $searchEngine (@searchEngines) {
	foreach my $condition (@conditions) {
		##############################################################
		## Search engine- and condition-dependent parameter setting	##
		##############################################################
		my $databaseName;
		my $pitFile;
		my $isolationWindow = 1.6;
		my $isolationOffset = 0.25;
		my $isolationVariation = 0.25;
		my $percentagePPI = 10;
		my $massCorrection = 2;
		my $ms2Deisotope = 1;
		my $ms2Consolidation = 10;
		my $TMTdata = 0;
		my $tagGeneration = 0;
		my $tagTolerance = 10;
		my $tagToleranceUnit = 2;
		my $ionSeries = "1 1 0 0 0 0 0 1 0";
		my $fragMassTolerance;
		my $fragMassToleranceUnit;
		my $ionLossesMS2 = "0 0 0 0";
		my $ionLossesMS1 = 0;
		my $secondSearch = 0;
		my $addNtermPeptide = "0.0000";
		my $addKLysine = "0.0000";
		my ($dynamicS, $dynamicT, $dynamicY);
		my $paramsDir;
		my $paramsName;
		## Database name and pit file
		if ($searchEngine eq "JUMP") {
			if ($condition eq "HH" || $condition eq "HL") {
				$databaseName = "/data1/database/20150201/human_ft_mc2_c0.fasta.mdx";
				$pitFile = "/data1/database/20150201/human_ft_mc2_c0.pit";
			} elsif ($condition eq "TMThh") {
				$databaseName = "/data1/database/20150201/human_ft_mc2_c0_TMT_K229.fasta.mdx";
				$pitFile = "/data1/database/20150201/human_ft_mc2_c0_TMT_K229.pit";
			} elsif ($condition eq "TMThhpho") {
				$databaseName = "/data1/database/20150201/human_ft_mc2_c0_pho_TMT_K229.fasta.mdx";
				$pitFile = "/data1/database/20150201/human_ft_mc2_c0_pho_TMT_K229.pit";
			}
		} else {
			if ($condition eq "HH" || $condition eq "HL") {
				$databaseName = "/data1/database/20150201/human_ft_mc2_c0.fasta.hdr";
				$pitFile = "/data1/database/20150201/human_ft_mc2_c0.pit";
			} elsif ($condition eq "TMThh") {
				$databaseName = "/data1/database/20150201/human_ft_mc2_c0_TMT_K229.fasta.hdr";
				$pitFile = "/data1/database/20150201/human_ft_mc2_c0_TMT_K229.pit";
			} elsif ($condition eq "TMThhpho") {
				$databaseName = "/data1/database/20150201/human_ft_mc2_c0_TMT_K229.fasta.hdr";
				$pitFile = "/data1/database/20150201/human_ft_mc2_c0_TMT_K229.pit";
			}
		}
		## Isolation window
		if ($condition eq "TMThh") {
			$isolationWindow = 1.0;
			$isolationOffset = 0.2;
			$isolationVariation = 0.2;
		} elsif ($condition eq "TMThhpho") {
			$isolationWindow = 1.4;
			$isolationOffset = 0;
			$isolationVariation = 0.2;
		}
		## PPI percentage
		if ($condition eq "TMThh" || $condition eq "TMThhpho") {
			$percentagePPI = 50;
		}
		## Mass correction and MS2-deisotope
		if ($condition eq "HL") {
			$massCorrection = 1;
			$ms2Deisotope = 0;
		}
		## MS2-consolidation
		if ($condition eq "TMThhpho") {
			$ms2Consolidation = 30;
		} elsif ($searchEngine eq "SEQUEST" && $condition eq "HL") {
			$ms2Consolidation = 20;
		}
		## SEQUEST-phospho-search
		if ($searchEngine eq "SEQUEST" && $condition eq "TMThhpho") {
			$ms2Deisotope = 0;
			$ms2Consolidation = 10000;
		}
		## TMT-data
		if ($condition eq "TMThh" || $condition eq "TMThhpho") {
			$TMTdata = 1;
		}
		## Tag generation
		if ($searchEngine eq "JUMP") {
			$tagGeneration = 1;
		}
		## Tag tolerance and ion series
		if ($condition eq "HL") {
			$tagTolerance = 0.4;
			$tagToleranceUnit = 1;
			$ionSeries = "0 1 0 0 0 0 0 1 0";
		}
		## Fragment ion mass tolerance
		if ($condition eq "HL") {
			$fragMassTolerance = 0.5;
			$fragMassToleranceUnit = 1;
		} else {
			if ($searchEngine eq "JUMP") {
				$fragMassTolerance = 15;
				$fragMassToleranceUnit = 2;
			} else {
				$fragMassTolerance = 0.015;
				$fragMassToleranceUnit = 1;
			}
		}	 
		## Ion losses and dynamic modifications for phosphoproteome
		if ($condition eq "TMThhpho") {
			$ionLossesMS2 = "1 1 1 0 ";
			$ionLossesMS1 = 1;
			$dynamicS = 79.96633;
			$dynamicT = 79.96633;
			$dynamicY = 79.96633;
		}
		## Second search
		if ($searchEngine eq "JUMP") {
			$secondSearch = 1;
		}
		## Static modification
		if ($TMTdata == 1) {
			$addNtermPeptide = 229.162932;
			$addKLysine = 229.162932;
		}
		## Params file name
		if ($searchEngine eq "JUMP") {
			$paramsName = "jump_sj_";
		} else {
			$paramsName = "jump_ss_";
		}
		if ($condition eq "HH") {
			#$paramsName = $paramsName."HH_human.params";
			$paramsName = $paramsName."HH.params";
		} elsif ($condition eq "HL") {
			$paramsName = $paramsName."HL.params";
		} elsif ($condition eq "TMThh") {
			$paramsName = $paramsName."TMThh.params";
		} elsif ($condition eq "TMThhpho") {
			$paramsName = $paramsName."TMThhpho.params";
		}
		
		##########################
		## Print .params files	##
		##########################
		if ($condition eq "HH") {
			$paramsDir = $currDir."/HH";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;
		} elsif ($condition eq "HL") {
			$paramsDir = $currDir."/HL";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;			
		} elsif ($condition eq "TMThh") {
			$paramsDir = $currDir."/TMThh";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;			
		} else {
			$paramsDir = $currDir."/TMThhpho";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;			
		}
		open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
		print PARAMS "# JUMP search parameter file (Version: $version, Date: $releaseDate)\n";
		print PARAMS "\n";
		print PARAMS "search_engine = $searchEngine										# use JUMP or SEQUEST for database search\n";
		print PARAMS "\n";
		print PARAMS "# Database settings\n";
		if ($condition eq "TMThh" or $condition eq "TMThhpho") {
			print PARAMS "database_name = $databaseName			# use .fasta.mdx for JUMP and .fasta.hdr for SEQUEST\n";
			print PARAMS "pit_file = $pitFile				# protein inference table (pit) for grouping proteins/genes\n";
		} else {
			print PARAMS "database_name = $databaseName				# use .fasta.mdx for JUMP and .fasta.hdr for SEQUEST\n";
			print PARAMS "pit_file = $pitFile						# protein inference table (pit) for grouping proteins/genes\n";
		}
		print PARAMS "peptide_tolerance = 6										# precursor mass (MH+) tolerance, default = 6 ppm after mass correction\n";
		print PARAMS "peptide_tolerance_units = 2									# 1 = Da; 2 = ppm\n";
		print PARAMS "\n";
		print PARAMS "# Preprocessing parameters\n";
		print PARAMS "first_scan_extraction = 5000									# the first scan number for search\n";
		print PARAMS "last_scan_extraction = 10000									# the last scan number for search, use a large number (e.g. 10E6) for full scans\n";
		print PARAMS "isolation_window = $isolationWindow										# +/- (isolation_window)/2 based on MS2 isolation window (e.g. 1.6 m/z)\n";
		print PARAMS "isolation_window_offset = $isolationOffset									# +/- isolation_window_offset based on MS2 isolation window offset (e.g. 0.25 m/z)\n";
		print PARAMS "isolation_window_variation = $isolationVariation								# +/- isolation_window_variation based on MS2 isolation window offset (e.g. 0.25 m/z)\n";
		print PARAMS "\n";
		print PARAMS "interscanppm = 15										# mass tolerance for interscan precursor identification\n";
		print PARAMS "intrascanppm = 10										# mass tolerance for intrascan isotopic decharging\n";
		print PARAMS "max_num_ppi = 0											# 0 = disable; 1-0 = max precursor ions selected for mixed MS2 search\n";
		print PARAMS "percentage_ppi = $percentagePPI										# minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0\n";
		print PARAMS "ppi_charge_0 = 1										# 0 = discard uncharged MS1 (charge = 0); 1 = manual charge assignment (+2 and +3)\n";
		print PARAMS "ppi_charge_1 = 1										# 0 = discard MS1 with charge +1; 1 = enable original charge +1\n";
		print PARAMS "mass_correction = $massCorrection										# 0 = no correction, 1 = MS1-based, 2 = MS1/2-based, 3 = manual correction\n";
		print PARAMS "prec_window = 3											# 0 = disable; 1-10 (Da) = mz windows for removing precursor ions\n";
		print PARAMS "MS2_deisotope = $ms2Deisotope										# 0 = disable; 1 = enable\n";
		print PARAMS "ppm = 10											# mass tolerance for MS2 decharging and deisotoping\n";
		print PARAMS "charge12_ppm = 15										# mass tolerance for merging different charged ions with the same mass\n";
		print PARAMS "ms2_consolidation = $ms2Consolidation										# maximal number of peaks retained within each 100-Da window\n";
		print PARAMS "TMT_data = $TMTdata											# 0 = disable; 1 = enable\n";
		print PARAMS "\n";
		print PARAMS "# Tagging\n";
		print PARAMS "tag_generation = $tagGeneration										# 0 = disable; 1 = enable to generate tags\n";
		print PARAMS "tag_tolerance = $tagTolerance										# mass tolerance for measuring peak distance for generating tags\n";
		print PARAMS "tag_tolerance_unit = $tagToleranceUnit										# 1 = Da; 2 = ppm\n";
		print PARAMS "tag_select_method = comb_p									# tag ranking: comb_p, hyper_p or rank_p\n";
		print PARAMS "\n";
		print PARAMS "# Database searching\n";
		print PARAMS "ion_series = $ionSeries									# a, b, c, d, v, w, x, y and z ions, respectively\n";
		print PARAMS "frag_mass_tolerance = $fragMassTolerance									# mass tolerance for MS2 ion matching\n";
		print PARAMS "frag_mass_tolerance_unit = $fragMassToleranceUnit									# 1 = Da; 2 = ppm\n";
		print PARAMS "ion_losses_MS2 = $ionLossesMS2									# 0 = disable; 1 = enable neutral loss of H2O, HPO3, H3PO4 and NH3, respectively\n";
		print PARAMS "ion_losses_MS1 = $ionLossesMS1										# 0 = disable; 1 = use precursor ion phosphate neutral loss to estimate #S/T phosphorylation\n";
		print PARAMS "ion_scoring = 1											# 1 = scoring product ions simultaneously; 2 = scoring ion series and charge states separately\n";
		print PARAMS "\n";
		print PARAMS "matching_method = hyper_p									# PSM scoring: comb_p, hyper_p, rank_p\n";
		print PARAMS "tag_search_method = 2										# 1 = exit when found; 2 = exhaustive search using tags defined by max_number_tag_for_search\n";
		print PARAMS "max_number_tags_for_search = 50									# max tags used for search unless the total number of tags is smaller than this defined value\n";
		print PARAMS "number_of_selected_result = 5									# maximal tentative PSMs in .spout ranked by Jscore\n";
		print PARAMS "number_of_detailed_result = 5									# maximal tentative PSMs in .spout ranked by pattern matching score\n";
		print PARAMS "second_search = $secondSearch										# 0 = disable; 1 = enable; for PSMs with FDR>0, perform the another round of search\n";
		print PARAMS "												# by relaxing monoisotopic mass by including M-2, M-1, M, M+1, M+2\n";
		print PARAMS "# Dynamic Modifications: SEQUEST requires no new database; but JUMP requires new database\n";
		print PARAMS "# C:  57.02146 carbamidomethylation or 71.0371 acrylamide\n";
		print PARAMS "# STY: 79.96633;  M: 15.99492; GG: 114.04293\n";
		print PARAMS "# SILAC: K:4.02511, 6.02013, 8.01420; SILAC: R:6.02013, 10.00827\n";
		print PARAMS "# TMT6-10: 229.1629321; TMT2: 225.1558327; TMT0: 224.1524779\n";
		print PARAMS "dynamic_M = 15.99492										# add each dynamic modification by one line, starting with dynamic_aa\n";
		print PARAMS "# dynamic_C = 57.02146										# add additional dynamic modification by line, starting with dynamic_aa\n";
		if (defined $dynamicS) {
			print PARAMS "dynamic_S = $dynamicS\n";
			print PARAMS "dynamic_T = $dynamicT\n";
			print PARAMS "dynamic_Y = $dynamicY\n";
		}
		print PARAMS "\n";
		print PARAMS "# Parameters for creating database (should match with the selected database)\n";
		print PARAMS "enzyme_info = Tryptic KR P									# LysC K ; ArgC R ; GluC DE ;\n";
		print PARAMS "digestion = full										# full or partial\n";
		print PARAMS "max_mis_cleavage = 2										# maximal miscleavage sites allowed for each peptide\n";
		print PARAMS "min_peptide_mass = 400.0000									# minimal mass of peptide database\n";
		print PARAMS "max_peptide_mass = 6000.0000									# maximal mass of peptide database\n";
		print PARAMS "max_modif_num = 3										# maximal modifications allowed for each peptide\n";
		print PARAMS "\n";
		print PARAMS "# Static Modification\n";
		print PARAMS "add_Nterm_peptide = $addNtermPeptide									# TMT modification or other amine labeling\n";
		print PARAMS "add_Cterm_peptide = 0.0000\n";
		print PARAMS "add_A_Alanine = 0.0000\n";
		print PARAMS "add_B_avg_NandD = 0.0000\n";
		print PARAMS "add_C_Cysteine = 0.0000										# Cys alkylation\n";
		print PARAMS "add_D_Aspartic_Acid = 0.0000\n";
		print PARAMS "add_E_Glutamic_Acid = 0.0000\n";
		print PARAMS "add_F_Phenylalanine = 0.0000\n";
		print PARAMS "add_G_Glycine = 0.0000\n";
		print PARAMS "add_H_Histidine = 0.0000\n";
		print PARAMS "add_I_Isoleucine = 0.0000\n";
		print PARAMS "add_J_user_amino_acid = 0.0000\n";
		if ($condition eq "TMThh" or $condition eq "TMThhpho") {
			print PARAMS "add_K_Lysine = $addKLysine									# TMT modification or other amine labeling\n";
		} else {
			print PARAMS "add_K_Lysine = $addKLysine										# TMT modification or other amine labeling\n";
		}
		print PARAMS "add_L_Leucine = 0.0000\n";
		print PARAMS "add_M_Methionine = 0.0000\n";
		print PARAMS "add_N_Asparagine = 0.0000\n";
		print PARAMS "add_O_Ornithine = 0.0000\n";
		print PARAMS "add_P_Proline = 0.0000\n";
		print PARAMS "add_Q_Glutamine = 0.0000\n";
		print PARAMS "add_R_Arginine = 0.0000\n";
		print PARAMS "add_S_Serine = 0.0000\n";
		print PARAMS "add_T_Threonine = 0.0000\n";
		print PARAMS "add_U_user_amino_acid = 0.0000\n";
		print PARAMS "add_V_Valine = 0.0000\n";
		print PARAMS "add_W_Tryptophan = 0.0000\n";
		print PARAMS "add_X_LorI = 0.0000\n";
		print PARAMS "add_Y_Tyrosine = 0.0000\n";
		print PARAMS "add_Z_avg_QandE = 0.0000\n";
		print PARAMS "\n";
		print PARAMS "# Other parameters\n";
		print PARAMS "simulation = 0											# 0 = disable; 1 = enable; this function used for testing the target-decoy strategy\n";
		print PARAMS "sim_MS1 = 1000											# ppm addition for MS1 decoys\n";
		print PARAMS "sim_MS2 = 5											# Da window for randomized MS2 peaks\n";
		print PARAMS "cluster = 1											# 0 = disable; 1 = enable; using master node only or entire cluster\n";
		print PARAMS "processors_used = 40										# \n";
		print PARAMS "Job_Management_System = SGE									# SGE used by current cluster; other systems (e.g. LSF & PBS) may be used\n";
		print PARAMS "temp_file_removal = 1										# 0 = disable (keep temporary files); 1 = enable (remove temporary files)\n";
		close (PARAMS);
	}
}

##############################################
## Filtering parameter generation (jump -f)	##
##############################################
foreach my $searchEngine (@searchEngines) {
	foreach my $condition (@conditions) {
		##############################################################
		## Search engine- and condition-dependent parameter setting	##
		##############################################################
		my $input;
		my $uniqueProteinPeptide = "protein";
		my $initialOutfileFDR = 5;
		my $FDR = 2;
		my $minXCorr;
		my $mods = 0;
		my $minOutfileNumXcorrFilter = 500;
		my ($oneHitMinXCorrZ1, $oneHitMinXCorrZ2, $oneHitMinXCorrZ3);
		my $paramsDir;
		my $paramsName;
		## Input name
		if ($condition eq "HH") {
			$input = "HH_human_".lc($searchEngine).":"."/data1/pipeline/release/version$version"."/SampleData/HH/HH_human_".lc($searchEngine)."/HH_human_".lc($searchEngine).".1";				
		} elsif ($condition eq "HL") {
			$input = "HL_human_".lc($searchEngine).":"."/data1/pipeline/release/version$version"."/SampleData/HL/HL_human_".lc($searchEngine)."/HL_human_".lc($searchEngine).".1";
		} elsif ($condition eq "TMThh") {
			$input = "HH_tmt10_human_".lc($searchEngine).":"."/data1/pipeline/release/version$version"."/SampleData/TMThh/HH_tmt10_human_".lc($searchEngine)."/HH_tmt10_human_".lc($searchEngine).".1";
		} else {
			$input = "HH_pho_tmt10_human_".lc($searchEngine).":"."/data1/pipeline/release/version$version"."/SampleData/TMThhpho/HH_pho_tmt10_human_".lc($searchEngine)."/HH_pho_tmt10_human_".lc($searchEngine).".1";
		}
		## Protein or peptide
		if ($condition eq "TMThhpho") {
			$uniqueProteinPeptide = "peptide";
		}
		## Initial outfile FDR
		if ($condition eq "TMThhpho") {
			$initialOutfileFDR = 10;
		}
		## FDR
		if ($condition eq "TMThhpho") {
			$FDR = 2;
		}
		## minXcorr
		if ($searchEngine eq "JUMP") {
			$minXCorr = 10;
		} else {
			$minXCorr = 1;
		}
		## Modifications
		if ($condition eq "TMThhpho") {
			$mods = "STY";
		}
		## Minimum outfile number for Xcorr filtering
		if ($condition eq "TMThhpho") {
			$minOutfileNumXcorrFilter = 200;
		}
		## One-hit wonders
		if ($searchEngine eq "JUMP") {
			$oneHitMinXCorrZ1 = 100;
			$oneHitMinXCorrZ2 = 25;
			$oneHitMinXCorrZ3 = 35;
		} else {
			$oneHitMinXCorrZ1 = 10;
			$oneHitMinXCorrZ2 = 2.5;
			$oneHitMinXCorrZ3 = 3.5;
		}
		## Params file name
		if ($searchEngine eq "JUMP") {
			$paramsName = "jump_fj_";
		} else {
			$paramsName = "jump_fs_";
		}
		if ($condition eq "HH") {
			$paramsName = $paramsName."HH.params";
		} elsif ($condition eq "HL") {
			$paramsName = $paramsName."HL.params";
		} elsif ($condition eq "TMThh") {
			$paramsName = $paramsName."TMThh.params";
		} elsif ($condition eq "TMThhpho") {
			$paramsName = $paramsName."TMThhpho.params";
		}
		
		##########################
		## Print .params files	##
		##########################
		if ($condition eq "HH") {
			$paramsDir = $currDir."/HH";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;
		} elsif ($condition eq "HL") {
			$paramsDir = $currDir."/HL";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;
		} elsif ($condition eq "TMThh") {
			$paramsDir = $currDir."/TMThh";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;
		} else {
			$paramsDir = $currDir."/TMThhpho";
			if (!-e $paramsDir) {
				system("mkdir $paramsDir");
			}
			$paramsName = $paramsDir."/".$paramsName;
		}
		open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
		print PARAMS "# JUMP filtering parameter file (Version: $version, Date: $releaseDate)\n";
		print PARAMS "\n";
		print PARAMS "# Grouping # Input of the program, using full path of the search result folders (containing dta files and out/spout files)\n";
		print PARAMS "# Information (e.g., search engine, pit table) will be directly parsed from jump.params in the input path below\n";
		print PARAMS "\n";
		print PARAMS "$input\n";
		print PARAMS "\n";
		print PARAMS "##### commonly adjusted parameters ############################################################################################################\n";
		print PARAMS "# PSMs filtered by (i) defined minimum filtering parameters, (ii) mass accuracy, and (iii) scoring\n";
		print PARAMS "# For scoring filtering, PSMs filtered by (i) initial FDR, (ii) peptide or protein categorization and mupltistep FDR filtering, and (iii) manual one hit-wonder-removal\n";
		print PARAMS "unique_protein_or_peptide = $uniqueProteinPeptide		# use protein or peptide FDR as cutoff\n";
		print PARAMS "initial_outfile_fdr = $initialOutfileFDR				# %initial FDR for score filtering; default = 5 (%)\n";
		print PARAMS "multistep_FDR_filtering = 1			# 0 = disabled; 1 = enabled\n";
		print PARAMS "FDR = $FDR						# %FDR for filtering peptides or one-hit-wonder proteins (fixed <1% FDR for proteins matched by two or more precursors)\n";
		print PARAMS "min_protein_SC = 1				# minimum spectral counts requirement for proteins\n";
		print PARAMS "one_hit_wonders_removal = 0			# keep or remove one hit wonders (-1: removal all, 0:no filter, 1: partial+fully, 2: fully)\n";
		print PARAMS "mods = $mods					# Display modified peptides and their unmodified (0:Off, K:Lys, STY: Phosphorylation, ...)\n";
		print PARAMS "modpairs = 0					# Show modified peptides pairs or just modified peptides (0:only modified peptides, 1:modified pairs)\n";
		print PARAMS "pit_file = 0					# absolute path of pit file: use updated pit file; 0 = using pit file in jump.params in the search folder\n";
		print PARAMS "\n";
		print PARAMS "# Minimum filtering parameters for removing low quality PSMs during reading process\n";
		print PARAMS "min_peptide_length = 7				# peptide length (6 can be used for small database)\n";
		print PARAMS "max_peptide_mis = 2				# maximal number of miscleavages allowed for one peptide, default=2\n";
		print PARAMS "max_peptide_mod = 3				# maximal number of modifications allowed for one peptide, M = 2, SILAC (KR) = 4, Ub = 3, Pho (STY) = 5\n";
		print PARAMS "peptide_mod_removal = 0				# 0: Off, C: Remove all C-modified peptides, STY: Remove all STY-modifed peptides\n";
		print PARAMS "peptide_aa_removal = 0				# 0: Off, M: Remove all M-containing peptides\n";
		print PARAMS "min_XCorr = $minXCorr					# XCorr (default = 1) or Jscore (default = 10)\n";
		print PARAMS "min_dCn = 0					# dCn or dJ\n";
		print PARAMS "mix_label = 0					# Remove mixed labeled peptides:  (0: None, KR: SILAC, C: ICAT, etc...)\n";
		print PARAMS "filter_contaminants = 0				# 0: Not used, 1: remove listed contaminants named with \"CON_\"\n";
		print PARAMS "12combinations = 1 1 1 1 1 1 1 1 0 0 0 0	# Trypticity and charge => FT1 FT2 FT3 FT4 PT1 PT2 PT3 PT4 NT1 NT2 NT3 NT4 # 1=yes, 0=no\n";
		print PARAMS "\n";
		print PARAMS "# Filtering PSMs by mass accuracy (no grouping, mass correction for each LC run) and matching scores (grouped)\n";
		print PARAMS "bypass_filtering = 0				# 0: NO, 1: YES bypasses all mass accuracy and dCn/XCorr filtering\n";
		print PARAMS "mass_accuracy = 1				# Mass accuracy filtering # 1=yes, 0=no\n";
		print PARAMS "mass_consideration = 1				# Mass consideration for accuracy filtering => 1:(MH), 2:(MH,MH+1), 3:(MH,MH+1,MH+2), 4:(MH,MH+1,MH+2,MH+3),\n";
		print PARAMS "						# 5:(MH,MH+1,MH+2,MH+3,MH+4), 6:(MH-1,MH,MH+1,MH+2), 7:(MH-2,MH-1,MH,MH+1,MH+2)\n";
		print PARAMS "sd_or_static = sd				# Mass accuracy cutoff based on experimental standard deviation (sd_or_static = sd)\n";
		print PARAMS "sd = 5						# or static ppm values (sd_or_static = static)\n";
		print PARAMS "static_cutoff = 6				# ppm\n";
		print PARAMS "static_cutoff_without_mass_calib = 10		# ppm; if not enough good scans, use this threshold for ppm cut without mass calibration\n";
		print PARAMS "\n";
		print PARAMS "# Filtering PSMs by matching scores (grouping by Peptide length; Trypticity; Mod; Miscleavage; Charge; deltaCn (0.01 step)\n";
		print PARAMS "# Sorting outfile into different dynamic groups until each group has sufficient outfiles (e.g. min_outfile_num_for_XCorr_filter = 500)\n";
		print PARAMS "# Filtering by assigned scan FDR based on unique peptides or proteins.\n";
		print PARAMS "# Removing false outfiles (if no charge state found, assign both +2 and +3 to make two dta files, same SC, same m/z, different charge state)\n";
		print PARAMS "# Grouping SC (1000) into unique SC (350), peptides (500), unique peptides (300), proteins (900), unique proteins (150), and protein groups (120)\n";
		print PARAMS "# Using SC FDR to predict unique protein/peptide FDR to shorten the filtering process\n";
		print PARAMS "FDR_filtering_method = group			# LDA or group (select one of the two filtering methods)\n";
		print PARAMS "min_outfile_num_for_XCorr_filter = $minOutfileNumXcorrFilter		# number of outfiles in each group for XCorr filtering; any number between 500 and 1000 is recomemded\n";
		print PARAMS "\n";
		print PARAMS "# Applyling additional filtering for one-hit-wonders\n";
		print PARAMS "one_hit_wonders_min_XCorr_z1 = $oneHitMinXCorrZ1		# minimum XCorr for peptides with charge state 1\n";
		print PARAMS "one_hit_wonders_min_XCorr_z2 = $oneHitMinXCorrZ2		# minimum XCorr for peptides with charge state 2\n";
		print PARAMS "one_hit_wonders_min_XCorr_z3 = $oneHitMinXCorrZ3		# minimum XCorr for peptides with charge state 3 or above\n";
		print PARAMS "one_hit_wonders_min_dCn = 0.1			# minimum dCn\n";
		print PARAMS "one_hit_wonders_mis = 1				# number of miscleavages allowed for one hit wonders\n";
		print PARAMS "one_hit_wonders_mods = 1			# number of modifications allowed for hit wonders, M = 1, SILAC (KR) = 3, Ub = 2, Pho (STY) = 4\n";
		print PARAMS "######################################################################################################################################################\n";
		print PARAMS "\n";
		print PARAMS "#################### To turn on pepXML generation ######################\n";
		print PARAMS "output_pepXML = 1\n";
		close (PARAMS);
	}
}

##########################################################
## Quantification parameter file generation (jump -q)	##
##########################################################
foreach my $searchEngine (@searchEngines) {	
	foreach my $condition (@conditions) {
		next if ($condition eq "HH" || $condition eq "HL");
		##############################################################
		## Search engine- and condition-dependent parameter setting	##
		##############################################################
		my $paramsDir;
		my $paramsName;
		my $idTxt;
		my $saveDir;
		if ($searchEngine eq "JUMP") {
			if ($condition eq "TMThh") {
					$paramsDir = $currDir."/TMThh";
				if (!-e $paramsDir) {
					system("mkdir $paramsDir");
				}
				$paramsName = $paramsDir."/jump_qj_HH_tmt10_human.params";
				$idTxt = "/data1/pipeline/release/version$version/SampleData/TMThh/sum_HH_tmt10_human_jump/ID.txt";
				$saveDir = "HH_tmt10_human_jump";
			} elsif ($condition eq "TMThhpho") {
				$paramsDir = $currDir."/TMThhpho";
				if (!-e $paramsDir) {
					system("mkdir $paramsDir");
				}
				$paramsName = $paramsDir."/jump_qj_HH_pho_tmt10_human.params";
				$idTxt = "/data1/pipeline/release/version$version/SampleData/TMThhpho/sum_HH_pho_tmt10_human_jump_mod/IDmod.txt";
				$saveDir = "HH_pho_tmt10_human_jump_mod";
			}
		} else {
			if ($condition eq "TMThh") {
					$paramsDir = $currDir."/TMThh";
				if (!-e $paramsDir) {
					system("mkdir $paramsDir");
				}
				$paramsName = $paramsDir."/jump_qs_HH_tmt10_human.params";
				$idTxt = "/data1/pipeline/release/version$version/SampleData/TMThh/sum_HH_tmt10_human_sequest/ID.txt";
				$saveDir = "HH_tmt10_human_sequest";
			} elsif ($condition eq "TMThhpho") {
				$paramsDir = $currDir."/TMThhpho";
				if (!-e $paramsDir) {
					system("mkdir $paramsDir");
				}
				$paramsName = $paramsDir."/jump_qs_HH_pho_tmt10_human.params";
				$idTxt = "/data1/pipeline/release/version$version/SampleData/TMThhpho/sum_HH_pho_tmt10_human_sequest_mod/IDmod.txt";
				$saveDir = "HH_pho_tmt10_human_sequest_mod";
			}
		}
		
		##########################
		## Print .params files	##
		##########################
		open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
		print PARAMS "# JUMP quantification parameter file (Version: $version, Date: $releaseDate)\n";
		print PARAMS "\n";
		print PARAMS "# Input:  ID.txt or IDmod.txt\n";
		print PARAMS "idtxt = $idTxt\n";
		if ($condition eq "TMThh") {
			print PARAMS "save_dir = $saveDir						# name of the directory for JUMPq results (prefix \"quan-\" will be added)\n";
		} elsif ($condition eq "TMThhpho") {
			print PARAMS "save_dir = $saveDir					# name of the directory for JUMPq results (prefix \"quan-\" will be added)\n";
		}	
		print PARAMS "ppi_filter = 50								# precursor peak intensity percentage threshold\n";
		print PARAMS "quan_method = MS2							# MS2 = reporter ions from  MS2; MS3 = reporter ions from MS3\n";
		print PARAMS "Top_sel = 1								# only for the MS3 method, to find MS3 scans for peptides assigned by MS2\n";
		print PARAMS "									# 1 = MS3 (used for quan) immediately after MS2 (used for ID)\n";
		print PARAMS "									# 5 = 5 MS3 followed by 5 MS2\n";
		print PARAMS "min_intensity_method = 2						# 1 = minimum, 2 = maximum, 3 = mean, 4 = median\n";
		print PARAMS "min_intensity_value = 5000						# Minimum intensity threshold\n";
		print PARAMS "ms2_ms3_tolerance = 0.02						# dalton, only for the MS3 method, to define mass tolerance for finding MS3\n";
		print PARAMS "\n";
		print PARAMS "# TMT10 reporter ions (126.127726;127.124761;127.131081;128.128116;128.134436;129.131471;129.137790;130.134825;130.141145;131.138180)\n";
		print PARAMS "# TMT8 reporter ions (126.127726;127.124761;127.131081;128.134436;129.131471;129.137790;130.141145;131.138180)\n";
		print PARAMS "# TMT6 reporter ions (126.127726;127.124761;128.134436;129.131471;130.141145;131.138180)\n";
		print PARAMS "# ITRAQ 4-plex reporter ions (114.1112;115.1082;116.1116;117.1149)\n";
		print PARAMS "# ITRAQ 8-plex reporter ions (113.1078;114.1112;115.1082;116.1116;117.1149;118.1120;119.1153;121.1220)\n";
		print PARAMS "itraq_reporters = 126.127726;127.124761;127.131081;128.128116;128.134436;129.131471;129.137790;130.134825;130.141145;131.138180\n";
		print PARAMS "itraq_peak_extraction_first_ppm = 25					# reporter ion initial mass tolerance +/- ppm\n";
		print PARAMS "itraq_peak_extraction_second_sd = 8					# SD used for identification of reporter ions\n";
		print PARAMS "itraq_peak_extraction_method = 1					# 1 = strongest intensity; 2 = closest to expected report ion mass;\n";
		print PARAMS "									# only if multiple peaks detected within mass tolerance\n";
		print PARAMS "report_impurity = 1							# 1 = Yes; 0 = No;\n";
		print PARAMS "impurity_matrix = /data1/pipeline/release/version$version/JUMPq/TMT10.ini	# impurity table for correction\n";
		print PARAMS "impurity_batch = Batch2							# Batch # of TMT reagents\n";
		print PARAMS "\n";
		print PARAMS "# The program allows multiple comparisons (e.g. comparison_groups_comp1, do not change the prefix \"comparison_groups_\")\n";
		print PARAMS "# For TMT10, input sig126, sig127N, sig127C, sig128N, sig128C, sig129N, sig129C, sig130N, sig130C, and sig131 (sorted by mass)\n";
		print PARAMS "# e.g. comparison_groups_twoGroups = sig126, sig127N, sig127C, sig128N, sig128C : sig129N, sig129C, sig130N, sig130C, sig131\n";
		print PARAMS "#      comparison_groups_threeGroups = sig126, sig127N, sig127C : sig128N, sig128C, sig129N : sig129C, sig130N, sig130C, sig131\n";
		print PARAMS "comparison_analysis = 1							# 1 = Yes; 0 = No; for group comparison\n";
		print PARAMS "comparison_groups_twoGroups = sig126, sig127N, sig127C, sig128N, sig128C : sig129N, sig129C, sig130N, sig130C, sig131\n";
		print PARAMS "\n";
		print PARAMS "Loading_Bias_Adjustment = 1						# 1 = Yes; 0 = No;\n";
		print PARAMS "min_SN_for_adjustment = 10						# define the minimal signal (SN ratio) used for adjustment\n";
		print PARAMS "percentage_outlier_removal = 10						# percenage of scans that are not included for each end\n";
		print PARAMS "loading_bias_correction_method = 1					# 1 = average; 2 = median;\n";
		print PARAMS "\n";
		print PARAMS "# Define variations by moving average\n";
		print PARAMS "# If auto_create_window is set as 1, automatically create 100 windows based on the total number of scans\n";
		print PARAMS "# If auto_create_window is set as 1, moving_window_size and move_window_step are not used\n";
		print PARAMS "auto_create_window = 0							# 1 = Yes; 0 = No; if setting 1, the scan # needs to be >500\n";
		print PARAMS "moving_window_size = 100						# A small number setting significantly increase the processing time\n";
		print PARAMS "moving_window_step = 20							# A small number setting significantly increase the processing time\n";
		print PARAMS "\n";
		print PARAMS "# Outlier removal by Dixon's Q test for summing up peptide and proteins\n";
		print PARAMS "outlier_threshold = Q95							# Input Q90, Q95 and Q99 (the larger number, the more stringent)\n";
		print PARAMS "FDR_method = BH								# Nine FDR methods can be used (i.e.  \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\",\"0\")\n";
		print PARAMS "\n";
		print PARAMS "publication_table = 1							# 1 = Yes, 0 = NO\n";
		print PARAMS "publication_comparison_group =  sig126, sig127N, sig127C, sig128N, sig128C, sig129N, sig129C, sig130N, sig130C, sig131\n";
		print PARAMS "publication_Loading_Bias_Adjustment = 0					# 1 = Yes, 0 = NO\n";
		close (PARAMS);
	}
}

##################################################
## Database parameter file generation (jump -d)	##
##################################################
my $paramsDir = $currDir;
if (!-e $paramsDir) {
	system ("mkdir $paramsDir");
}
my $paramsName = $paramsDir."/jump_d.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP database generation parameter file (Version: $version, Date: $releaseDate)\n";
print PARAMS "# JUMPd generates a search database and/or a protein inference table (PIT) using the input databases (.fasta files)\n";
print PARAMS "# The type of a database (for JUMP or SEQUEST) solely depends on the search_engine information in \"jump.params\" file below\n";
print PARAMS "# All the modification information for the database also depends on the \"jump.params\" file\n";
print PARAMS "\n";
print PARAMS "# Parameters for generating a new database ######################################################################################################################################\n";
print PARAMS "input_database1 = /data1/database/20150201/HUMAN.fasta						# Absolute path of input_database1 (i.e. .fasta file)\n";
print PARAMS "#input_database2 = /data1/database/20150201/MOUSE.fasta						# Absolute path of input_database2\n";
print PARAMS "												# More databases can be used (e.g. input_database3 = /data1/database/20150201/ECOLI.fasta)\n";
print PARAMS "output_prefix = mouse										# Prefix for a new database (and .pit) file\n";
print PARAMS "include_contaminants = 1									# 0 = do not include contaminants; 1 = include contaminants\n";
print PARAMS "input_contaminants = /data1/database/contaminants.fasta						# Absolute path of a .fasta file containing contaminants\n";
print PARAMS "decoy_generation = 1										# 0 = do not include decoys; 1 = include decoys\n";
print PARAMS "decoy_generation_method = 1									# 1 = reverse; 2 = reverse and switch every K/R with its preceding AA\n";
print PARAMS "jump.params = /usr/search/jump_search.params							# Only search engine and modification information will be obtained from the file\n";
print PARAMS "												# \"database_name\" in a jump.params file will be ignored\n";
print PARAMS "bypass_db_generation = 0									# 0 = generate database, 1 = bypass the generation of database (only pit will be generated)\n";
print PARAMS "\n";
print PARAMS "# Parameters for generating a protein inference table (PIT) (keep the prefix \"list_\") ###########################################################################################\n";
#print PARAMS "list_protein_abundance1 = /data1/database/KnowledgeTables/ProteinAbundance/Mouse_Abundance_emPAI.txt	# Absolute path of a file containing protein abundance information\n";
print PARAMS "list_protein_abundance1 = /data1/database/KnowledgeTables/ProteinAbundance/Human_Abundance_emPAI.txt	# Absolute path of a file containing protein abundance information\n";
print PARAMS "												# More abundance information can be used\n";
print PARAMS "												# (e.g. list_protein_abundance2 = /data1/database/Rat_Abundance.txt)\n";
print PARAMS "list_TFs = /data1/database/KnowledgeTables/TFs/tfs_from_TRANSFAC.txt				# Absolute path of a file containing protein annotation information, TFs (transcription factors)\n";
print PARAMS "list_oncogenes = /data1/database/KnowledgeTables/Oncogenes/oncogenes_from_literatures.txt	# Absolute path of a file containing protein annotation information, oncogenes\n";
print PARAMS "list_kinases = /data1/database/KnowledgeTables/Kinases/kinases_from_pkinfam.txt			# Absolute path of a file containing protein annotation information, kinases\n";
print PARAMS "list_GPCRs = /data1/database/KnowledgeTables/GPCRs/gpcrs.txt					# Absolute path of a file containing protein annotation information, GPCRs\n";
print PARAMS "list_epigenetic_factors = /data1/database/KnowledgeTables/EpigeneticRegulators/epigenetic_regulators.txt\n";
print PARAMS "list_spliceosomal_proteins = /data1/database/KnowledgeTables/SpliceosomalProteins/spliceosomal_proteins.txt\n";
close (PARAMS);

##################################################################
## Absolute quantification parameter file generation (jump -aq)	##
##################################################################
$paramsName = $paramsDir."/jump_aq.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP absolute quantification parameter file (Version: $version, Date: $releaseDate)\n";
print PARAMS "# JUMPaq generates a protein abundance table based on a search/filtering result\n";
print PARAMS "\n";
print PARAMS "# Database (i.e. a fasta file)\n";
print PARAMS "database = /data1/database/20150201/HUMAN.fasta\n";
print PARAMS "\n";
print PARAMS "# The Report file generated by JUMPf\n";
print PARAMS "report_file = /data1/pipeline/release/version$version/SampleData/HH/sum_HH_human_jump/publications/id_all_prot.txt\n";
print PARAMS "\n";
print PARAMS "# Output file name\n";
print PARAMS "output = protein_absolute_abundance.txt\n";
print PARAMS "\n";
print PARAMS "# Parameters used for calculating theoretical peptides\n";
print PARAMS "min_peptide_length = 6\n";
print PARAMS "enzyme = trypsin			#types of enzyme: trypsin, chymotrypsin, gluc_nahpo, gluc, lysc, argc, aspc\n";
print PARAMS "max_peptide_length = 30\n";
print PARAMS "max_peptide_hydro = 17\n";
print PARAMS "min_peptide_hydro = -24\n";
close (PARAMS);

##################################################################
## Statistical inference parameter file generation (jump -i)	##
##################################################################
$paramsName = $paramsDir."/jump_i_default.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP statistical inference parameter file (Version: $version, Date: $releaseDate)\n";
print PARAMS "\n";
print PARAMS "#-----------------------------------------------------------------------------------------------------------------------------------\n";
print PARAMS "# Expected results:\n";
print PARAMS "#\n";
print PARAMS "# 1) log 2 ratio distribution for each pair of samples:\n";
print PARAMS "# 	a) log 2 ratio distributions figure: log2ratio_distributions.pdf\n";
print PARAMS "# 	b) log 2 ratio distribution mean: log2ratio_mean_matrix.txt (typical mean for replicates is around 0)\n";
print PARAMS "# 	c) log 2 ratio distribution SD: log2ratio_SD_matrix.txt (typical SD for replicates is 0.2)\n";
print PARAMS "#\n";
print PARAMS "# 2) Heatmap for sample cluster: heatmap_XXpct.pdf (default cutting %: 0.1, 0.5, 1, 5, 10 and 20)\n";
print PARAMS "#\n";
print PARAMS "#-----------------------------------------------------------------------------------------------------------------------------------\n";
print PARAMS "\n";
print PARAMS "# input file (either protein or peptide table from jump -q results in folder 'publications')\n";
print PARAMS "input_table = /data1/pipeline/release/version11.2.0/SampleData/TMThh/quan_HH_tmt10_human_jump/publications/id_uni_prot_quan.txt		# absolute path needed;\n"; 
print PARAMS "										# protein file: /home/user/project/quan_test1/publications/id_uni_prot_quan.txt\n";
print PARAMS "										# peptide file: /home/user/project/quan_test1/publications/quan_phosphopep_uni_prot.txt\n";
print PARAMS "\n";
print PARAMS "# sample labels\n";
print PARAMS "sig126 = condition1a\n";
print PARAMS "sig127N = condition1b\n";
print PARAMS "sig127C = condition2a\n";
print PARAMS "sig128N = condition2b\n";
print PARAMS "sig128C = condition3a\n";
print PARAMS "sig129N = condition3b\n";
print PARAMS "sig129C = condition4a\n";
print PARAMS "sig130N = condition4b\n";
print PARAMS "sig130C = condition5a\n";
print PARAMS "sig131 = condition5b\n";
print PARAMS "\n";
print PARAMS "# output folder\n";
print PARAMS "output_folder = test1				# output folder suffix name; prefix always 'inf_'\n";
print PARAMS "\n";
print PARAMS "# contaminant removel option\n";
print PARAMS "remove_contaminants = 1				# remove contaminants? 1 = yes; 0 = no\n";
print PARAMS "\n";
print PARAMS "# SD / mean matrix: analysis based on the stable proteins / peptides\n";
print PARAMS "pair_cutoff_percentage = 10			# percentage/2 of data will be ignored on each end of the distribution\n";
print PARAMS "pair_cutoff_intensity = 3000			# median intensity cutoff across samples\n";
print PARAMS "\n";
print PARAMS "# heatmap / sample clustering: analysis based on the most variable proteins / peptides\n";
print PARAMS "cluster_cutoff_percentage = 0.5			# percentage of most variable proteins (defined by RSD of intensity) used for clustering\n";
print PARAMS "cluster_cutoff_intensity = 3000			# median intensity cutoff across samples\n";
close (PARAMS);

$paramsName = $paramsDir."/jump_i_custom.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP statistical inference parameter file (Version: $version, Date: $releaseDate)\n";
print PARAMS "\n";
print PARAMS "#-----------------------------------------------------------------------------------------------------------------------------------\n";
print PARAMS "# Expected results:\n";
print PARAMS "#\n";
print PARAMS "# 1) log 2 ratio distribution for each pair of samples:\n";
print PARAMS "# 	a) log 2 ratio distributions figure: log2ratio_distributions.pdf\n";
print PARAMS "# 	b) log 2 ratio distribution mean: log2ratio_mean_matrix.txt (typical mean for replicates is around 0)\n";
print PARAMS "# 	c) log 2 ratio distribution SD: log2ratio_SD_matrix.txt (typical SD for replicates is 0.2)\n";
print PARAMS "#\n";
print PARAMS "# 2) Heatmap for sample cluster: heatmap_XXpct.pdf (default cutting %: 0.1, 0.5, 1, 5, 10 and 20)\n";
print PARAMS "#\n";
print PARAMS "#-----------------------------------------------------------------------------------------------------------------------------------\n";
print PARAMS "\n";
print PARAMS "# customized input\n";
print PARAMS "# 1) Excel: make an Excel sheet with format below:\n";
print PARAMS "# customID sig126 sig127N ... sig131\n";
print PARAMS "# 2) Excel: save the file as 'Text (Tab delimited)'\n";
print PARAMS "# 3) WinSCP: upload the saved file to cluster using WinSCP\n";
print PARAMS "# 4) JUMPi: turn off all filtering in jump -i parameter file (default)\n";
print PARAMS "input_table = /home/yli4/development/JUMPi/v11.3.042915/cus_prot1.txt				# absolute path of customized input file\n";
print PARAMS "\n";
print PARAMS "# sample labels\n";
print PARAMS "sig126 = condition1a\n";
print PARAMS "sig127N = condition1b\n";
print PARAMS "sig127C = condition2a\n";
print PARAMS "sig128N = condition2b\n";
print PARAMS "sig128C = condition3a\n";
print PARAMS "sig129N = condition3b\n";
print PARAMS "sig129C = condition4a\n";
print PARAMS "sig130N = condition4b\n";
print PARAMS "sig130C = condition5a\n";
print PARAMS "sig131 = condition5b\n";
print PARAMS "\n";
print PARAMS "# output folder\n";
print PARAMS "output_folder = cu_test1				# output folder suffix name; prefix always 'inf_'\n";
print PARAMS "\n";
print PARAMS "# contaminant removel option\n";
print PARAMS "remove_contaminants = 0				# remove contaminants? 1 = yes; 0 = no\n";
print PARAMS "\n";
print PARAMS "# SD / mean matrix: analysis based on the stable proteins / peptides\n";
print PARAMS "pair_cutoff_percentage = 0			# percentage/2 of data will be ignored on each end of the distribution\n";
print PARAMS "pair_cutoff_intensity = 0			# median intensity cutoff across samples\n";
print PARAMS "\n";
print PARAMS "# heatmap / sample clustering: analysis based on the most variable proteins / peptides\n";
print PARAMS "cluster_cutoff_percentage = 100			# percentage of most variable proteins (defined by RSD of intensity) used for clustering\n";
print PARAMS "cluster_cutoff_intensity = 0			# median intensity cutoff across samples\n";
close (PARAMS);

######################################################################
## Validation of known targets parameter file generation (jump -v)	##
######################################################################
$paramsName = $paramsDir."/jump_v.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP validation (of known targets) parameter file (Version: $version, Date: $releaseDate)\n";
print PARAMS "# JUMPv extracts the identification/quantification results of input peptides/proteins\n";
print PARAMS "\n";
print PARAMS "# Peptides or proteins to be validated\n";
print PARAMS "peptide_1 = VRHDSPDPSPPR\n";
print PARAMS "peptide_2 = SSDEDATGEPK\n";
print PARAMS "#peptide_3 = VMSSLAPYNSSTSPQK\n";
print PARAMS "protein_1 = sp|Q9BRD0|BUD13_HUMAN\n";
print PARAMS "#protein_2 = tr|B2KFM4|B2KFM4_MOUSE\n";
print PARAMS "\n";
print PARAMS "# The raw quantification file generated by JUMPq\n";
print PARAMS "quan_file = /data1/pipeline/release/version$version/SampleData/TMThhpho/quan_HH_pho_tmt10_human_jump_mod/raw_quan_HH_pho_tmt10_human_jump_mod_scan.txt\n";
print PARAMS "\n";
print PARAMS "# Output file name\n";
print PARAMS "output = validation_quan_HH_pho_tmt10_human_jump_mod.txt\n";
close (PARAMS);

###################################
## customized database (jump -c) ##
###################################
$paramsName = $paramsDir."/jump_c.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP customized database (jump -c)\n";
print PARAMS "mode = junction                         # mutation, junction, 6FT\n";
print PARAMS "input_file = /home/yli4/development/JUMPg/wrapUp_052915/junction_simple_input.txt       # absolute path of input file\n";
print PARAMS "output_folder = correctChr_jun_test12\n";
print PARAMS "species = human                 # human; required for 'mutation' and 'junction' mode\n";
print PARAMS "read_coverage_cutoff = 2                # minimum number of reads that support a candidate peptide; required for '6FT' mode\n";
close (PARAMS);

##############################
## spectrum quality control ##
##############################
$paramsName = $paramsDir."/jump_qc.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP spectrum quality control (jump -qc)\n";
print PARAMS "run1 = /home/yli4/customizedDB/AD/JUMPsearch_032015/multistage/ad_pl01/ad_pl01.5\n";
print PARAMS "confident_IDtxt = /home/yli4/customizedDB/AD/JUMPsearch_032015/uniPro/sum_10fr_peptide0FDR_ppi1_test2/ID.txt\n";
print PARAMS "accepted_IDtxt = /home/yli4/customizedDB/AD/JUMPsearch_032015/multistage/sum_proNtermAcetyl_test1/ID.txt\n";
print PARAMS "PSM_recoveray_rate = 99\n";
print PARAMS "output_folder = after_proNtermAcetyl_test1\n";
close (PARAMS);

##################################################
## Localization of modification sites (jump -l)	##
##################################################

$paramsName = $paramsDir."/jump_l.params";
open (PARAMS, ">", $paramsName) or die "Cannot generate $paramsName\n";
print PARAMS "# JUMP localization parameter file (Version: $version, Date: $releaseDate)\n";
print PARAMS "# JUMPl identifies (localizes) the modification sites (e.g. phosphorylation)\n"; 
print PARAMS "\n";
print PARAMS "# Input and output settings\n";
print PARAMS "IDmod = /home/htan/1412Gilbertson/pho/sum_pho_NSC_20150206_mod/IDmod.txt	# path for the input IDmod.txt file\n";
print PARAMS "#IDmod = /home/xwang4/JUMPl_testing/test.txt\n";
print PARAMS "Output = ID.lscore								# output file name\n";
print PARAMS "\n";
print PARAMS "# Preprocessing parameters\n";
print PARAMS "peptide_score_tolerance = 10							# tolerance for peptide localization score (percentage)\n";
print PARAMS "mass_correction = 2								# 0 = no correction, 1 = MS1-based, 2 = MS2-based, 3 = manual correction\n";                                                                          
print PARAMS "isolation_window = 1.2								# +/- (isolation_window)/2 based on MS2 isolation window (e.g. 1.2 m/z)\n";
print PARAMS "first_scan_extraction = 0							# the first scan number for search\n";
print PARAMS "last_scan_extraction = 1000000							# the last scan number for search, use a large number (e.g. 10E6) for full scans\n";
print PARAMS "MS2_deisotope = 1								# 0 = disable; 1 = enable\n";
print PARAMS "ppm = 10									# mass tolerance for MS2 decharging and deisotoping\n";
print PARAMS "ms2_consolidation = 10								# maximal number of peaks retained within each 100-Da window\n";
print PARAMS "\n";
print PARAMS "# Peptide selection and scoring options\n";
print PARAMS "ion_series = 0 1 0 0 0 0 0 1 0							# a, b, c, d, v, w, x, y and z ions, respectively\n";
print PARAMS "ion_losses_MS2 = 1 0 0 0							# 0 = disable; 1 = enable neutral loss of H2O, HPO3, H3PO4 and NH3, respectively\n"; 
print PARAMS "frag_mass_tolerance = 10							# mass tolerance for MS2 ion matching\n";
print PARAMS "frag_mass_tolerance_unit = 2							# 1 = Da; 2 = PPM;\n";
print PARAMS "\n";
print PARAMS "# Dynamic modifications\n";
print PARAMS "dynamic_M = 15.99492								# add each dynamic modification by one line, starting with dynamic_AA\n"; 
print PARAMS "dynamic_S = 79.96633\n";
print PARAMS "dynamic_T = 79.96633\n";
print PARAMS "dynamic_Y = 79.96633\n";
print PARAMS "max_modif_num = 3\n";
print PARAMS "\n";
print PARAMS "# Static Modification\n";
print PARAMS "add_Nterm_peptide = 229.1629321							# TMT modification or other amine labeling\n";
print PARAMS "add_Cterm_peptide = 0.0000\n";
print PARAMS "add_A_Alanine = 0.0000\n";
print PARAMS "add_B_avg_NandD = 0.0000\n";
print PARAMS "add_C_Cysteine = 0.0000								# Cys alkylation\n";
print PARAMS "add_D_Aspartic_Acid = 0.0000\n";
print PARAMS "add_E_Glutamic_Acid = 0.0000\n";
print PARAMS "add_F_Phenylalanine = 0.0000\n";
print PARAMS "add_G_Glycine = 0.0000\n";
print PARAMS "add_H_Histidine = 0.0000\n";
print PARAMS "add_I_Isoleucine = 0.0000\n";
print PARAMS "add_J_user_amino_acid = 0.0000\n";
print PARAMS "add_K_Lysine = 229.1629321							# TMT modification or other amine labeling\n";
print PARAMS "add_L_Leucine = 0.0000\n";
print PARAMS "add_M_Methionine = 0.0000\n";
print PARAMS "add_N_Asparagine = 0.0000\n";
print PARAMS "add_O_Ornithine = 0.0000\n";
print PARAMS "add_P_Proline = 0.0000\n";
print PARAMS "add_Q_Glutamine = 0.0000\n";
print PARAMS "add_R_Arginine = 0.0000\n";
print PARAMS "add_S_Serine = 0.0000\n";
print PARAMS "add_T_Threonine = 0.0000\n";
print PARAMS "add_U_user_amino_acid = 0.0000\n";
print PARAMS "add_V_Valine = 0.0000\n";
print PARAMS "add_W_Tryptophan = 0.0000\n";
print PARAMS "add_X_LorI = 0.0000\n";
print PARAMS "add_Y_Tyrosine = 0.0000\n";
print PARAMS "add_Z_avg_QandE = 0.0000\n";
close (PARAMS);










