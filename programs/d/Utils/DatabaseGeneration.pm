#!/usr/bin/perl

## Release date: 05/01/2015
## Release version: version 12.0.0

package Utils::DatabaseGeneration;
use strict;
use warnings;
use File::Basename;
use Cwd;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();
my $JUMPscript="/home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/s/jump.pl";

sub new {
	my ($class, %arg)=@_;
    my $self = {};
    bless $self, $class;
	return $self;
}

sub generateDb {
	my ($self, $jumpParams, $dbName, $log) = @_;
	if (lc($$jumpParams{'search_engine'}) eq "jump") {
		######################################
		## JUMP search database generation	##
		######################################
		generateJumpParams($jumpParams, $dbName, $log);
		print "JUMP search database is being generated\n";
		print $log "JUMP search database is being generated\n";
		#system("perl /data1/pipeline/release/version12.0.0/JUMPsj/jump.pl -p jumpdb.params");
		#system("perl /home/yli4/customizedDB/AD/DB/AD_labelfree_pt/jump.pl -p jumpdb.params");
		#system("perl /home/yli4/customizedDB/AD/DB/pepNterm_Gln2pyroGlu/JUMPsj/jump.pl -p jumpdb.params");
		#system("perl /home/yli4/development/JUMPg/v1.0/programs/JUMPsj/jump.pl -p jumpdb.params");
		#system("perl $JUMPscript -p jumpdb.params");
		system("$JUMPscript -p jumpdb.params");
		system("rm jumpdb.params");
		my $generatedDb = $dbName.".fasta.mdx";
		if (!-e $generatedDb) {
			print "Problem in generating a JUMP database $generatedDb\n";
			print "Please check builddb.params file carefully\n";
			print $log "Problem in generating a JUMP database $generatedDb\n";
			print $log "Please check builddb.params file carefully\n";
			
			exit;
		} else {
			print "JUMP search database $generatedDb has been successfully generated\n";
			print $log "JUMP search database $generatedDb has been successfully generated\n";
		}
	} elsif (lc($$jumpParams{'search_engine'}) eq "sequest") {
		##########################################
		## SEQUEST search database generation	##
		##########################################
		generateMakedbParams($jumpParams, $dbName, $log);
		print "SEQUEST search database is being generated\n";
		print $log "SEQUEST search database is being generated\n";
		system("makedb");
		my $generatedDb = $dbName.".fasta.hdr";
		if (!-e $generatedDb) {
			print "Problem in generating a SEQUEST database $generatedDb\n";
			print "Please check builddb.params file carefully\n";
			print $log "Problem in generating a SEQUEST database $generatedDb\n";
			print $log "Please check builddb.params file carefully\n";
			exit;
		}
		## Test the generated database using SEQUEST
		generateSequestParams($jumpParams, $generatedDb, $log);
		system ("cp /data1/pipeline/release/version12.0.0/JUMPd/test.dta .");
		my $sequestResult = `sequest28single -Psequest.params test.dta`;
		system ("rm makedb.params sequest.params test.dta test.out");
		if ($sequestResult =~ /Rank\/Sp/) {
			print "Testing $generatedDb is done\n";
			print "SEQUEST search database $generatedDb has been successfully generated\n";
			print $log "Testing $generatedDb is done\n";
			print $log "SEQUEST search database $generatedDb has been successfully generated\n";
		} else {
			print "Error in testing the SEQUEST database $generatedDb\n";
			print "Please check builddb.params and jump.params files carefully\n";
			print $log "Error in testing the SEQUEST database $generatedDb\n";
			print $log "Please check builddb.params and jump.params files carefully\n";
			exit;
		}
	} else {
		print "Please set \"search_engine\" parameter in the jump.params file correctly\n";
		print "Either JUMP or SEQUEST\n";
		print $log "Please set \"search_engine\" parameter in the jump.params file correctly\n";
		print $log "Either JUMP or SEQUEST\n";
		exit;
	}
}

sub generateJumpParams {
	## $params is from jump.params
	my ($params, $dbName, $log) = @_;
	my $dbPath = getcwd()."/".$dbName.".fasta";
	if (!-e $dbPath) {
		print "Cannot find $dbPath file for the database generation\n";
		print "Please check the location of $dbName".".fasta file\n";
		print $log "Cannot find $dbPath file for the database generation\n";
		print $log "Please check the location of $dbName".".fasta file\n";
		exit;
	}
	open (PARFILE, ">", "jumpdb.params");
	foreach my $key (keys %{$params}) {
		if ($key =~ /^search_engine/) {
			print PARFILE "search_engine = JUMP\n";
		} elsif ($key =~ /^database_name/) {
			print PARFILE "database_name = ".$dbPath."\n";
=head
		#} elsif ($key =~ /^digestion/ && $$params{$key} ne "full") {
			print "JUMP only considers a fully-tryptic condition\n";
			print "Please check jump.params file and run again\n";
			print $log "JUMP only considers a fully-tryptic condition\n";
			print $log "Please check jump.params file and run again\n";
			exit;
=cut
		} else {
			print PARFILE $key." = ".$$params{$key}."\n";
		}
	}
	close (PARFILE);
}

sub generateMakedbParams {
	## $params is from jump.params
	my ($params, $dbName, $log) = @_;
	my $dbPath = getcwd()."/".$dbName.".fasta";
	open (PARFILE, ">", "makedb.params");
	print PARFILE "\[MAKEDB\]\n";
	print PARFILE "database_name = $dbPath\n";
	print PARFILE "intermediate_directory =\n";
	print PARFILE "sort_directory = \/tmp\/\n";
	print PARFILE "sort_program = \/bin\/sort\n";
	if ($$params{'enzyme_info'} eq "Tryptic KR P") {
		if ($$params{'digestion'} eq "full") {
			print PARFILE "enzyme_info = Trypsin 1 1 KR P\n";
		} elsif ($$params{'digestion'} eq "partial") {
			print PARFILE "enzyme_info = Trypsin 2 1 KR P\n";
		} else {
			print "Please set \"digestion\" parameter correctly; either full or partial\n";
			print $log "Please set \"digestion\" parameter correctly; either full or partial\n";
			exit;
		}
	} else {
		print PARFILE "enzyme_info = ".$$params{'enzyme_info'}."\n";
	}
	print PARFILE "protein_or_nucleotide_dbase = 0\n";
	print PARFILE "nucleotide_reading_frames = 0\n";
	print PARFILE "use_mono\/avg_masses = 1\n";
	print PARFILE "min_peptide_mass = ".$$params{'min_peptide_mass'}."\n";
	print PARFILE "max_peptide_mass = ".$$params{'max_peptide_mass'}."\n";
	print PARFILE "max_num_internal_cleavage_sites = ".$$params{'max_mis_cleavage'}."\n";
	print PARFILE "\[STATIC_MODIFICATIONS\]\n";
	print PARFILE "add_Cterm_peptide = ".$$params{'add_Cterm_peptide'}."\n";
	print PARFILE "add_Nterm_peptide = ".$$params{'add_Nterm_peptide'}."\n";
	print PARFILE "add_G_Glycine = ".$$params{'add_G_Glycine'}."\n";
	print PARFILE "add_A_Alanine = ".$$params{'add_A_Alanine'}."\n";
	print PARFILE "add_S_Serine = ".$$params{'add_S_Serine'}."\n";
	print PARFILE "add_P_Proline = ".$$params{'add_P_Proline'}."\n";
	print PARFILE "add_V_Valine = ".$$params{'add_V_Valine'}."\n";
	print PARFILE "add_T_Threonine = ".$$params{'add_T_Threonine'}."\n";
	print PARFILE "add_C_Cysteine = ".$$params{'add_C_Cysteine'}."\n";
	print PARFILE "add_L_Leucine = ".$$params{'add_L_Leucine'}."\n";
	print PARFILE "add_I_Isoleucine = ".$$params{'add_I_Isoleucine'}."\n";
	print PARFILE "add_X_LorI = ".$$params{'add_X_LorI'}."\n";
	print PARFILE "add_N_Asparagine = ".$$params{'add_N_Asparagine'}."\n";
	print PARFILE "add_O_Ornithine = ".$$params{'add_O_Ornithine'}."\n";
	print PARFILE "add_B_avg_NandD = ".$$params{'add_B_avg_NandD'}."\n";
	print PARFILE "add_D_Aspartic_Acid = ".$$params{'add_D_Aspartic_Acid'}."\n";
	print PARFILE "add_Q_Glutamine = ".$$params{'add_Q_Glutamine'}."\n";
	print PARFILE "add_K_Lysine = ".$$params{'add_K_Lysine'}."\n";
	print PARFILE "add_Z_avg_QandE = ".$$params{'add_Z_avg_QandE'}."\n";
	print PARFILE "add_E_Glutamic_Acid = ".$$params{'add_E_Glutamic_Acid'}."\n";
	print PARFILE "add_M_Methionine = ".$$params{'add_M_Methionine'}."\n";
	print PARFILE "add_H_Histidine = ".$$params{'add_H_Histidine'}."\n";
	print PARFILE "add_F_Phenylalanine = ".$$params{'add_F_Phenylalanine'}."\n";
	print PARFILE "add_R_Arginine = ".$$params{'add_R_Arginine'}."\n";
	print PARFILE "add_Y_Tyrosine = ".$$params{'add_Y_Tyrosine'}."\n";
	print PARFILE "add_W_Tryptophan = ".$$params{'add_W_Tryptophan'}."\n";
	print PARFILE "add_J_user_amino_acid = ".$$params{'add_J_user_amino_acid'}."\n";
	print PARFILE "add_U_user_amino_acid = ".$$params{'add_U_user_amino_acid'}."\n";
	close PARFILE;	
}

sub generateSequestParams {
	my ($params, $dbName) = @_;
	my $dbPath = getcwd()."/".$dbName;
	open (PARFILE, ">", "sequest.params");
	print PARFILE "first_database_name = $dbPath\n";
	print PARFILE "peptide_mass_tolerance = ".$$params{'peptide_tolerance'}."   # +/- peptide mass tolerance\n";
	print PARFILE "peptide_mass_units = ".$$params{'peptide_tolerance_units'}." # 0 = amu, 1 = mmu, 2 = ppm\n";
	print PARFILE "ion_series = 0 0 0 ".$$params{'ion_series'}."\n";
	print PARFILE "fragment_ion_tolerance = ".$$params{'frag_mass_tolerance'}."  # Da\n";
	print PARFILE "num_output_lines = 10\n";
	print PARFILE "num_results = 250\n";	
	print PARFILE "num_description_lines = 5\n";
	print PARFILE "show_fragment_ions = 0\n";
	print PARFILE "print_duplicate_references = 1\n";
	if ($$params{'enzyme_info'} eq "Tryptic KR P" && $$params{'digestion'} eq "full") {
		print PARFILE "enzyme_info = Trypsin 1 1 KR P\n";
	} elsif ($$params{'enzyme_info'} eq "Tryptic KR P" && $$params{'digestion'} eq "partial") {
		print PARFILE "enzyme_info = Trypsin 2 1 KR P\n";
	} else {
		print PARFILE "enzyme_info = ".$$params{'enzyme_info'}."\n";
	}
	print PARFILE "max_num_differential_per_peptide = ".$$params{'max_modif_num'}."\n";
	print PARFILE "diff_search_options = 15.99492 M\n";
	print PARFILE "term_diff_search_options = 0.000 0.000\n";
	print PARFILE "mass_type_parent = 1\n";
	print PARFILE "mass_type_fragment = 1\n";
	print PARFILE "normalize_xcorr = 0\n";
	print PARFILE "remove_precursor_peak = 0\n";
	print PARFILE "ion_cutoff_percentage = 0.0000\n";
	print PARFILE "max_num_internal_cleavage_sites = ".$$params{'max_mis_cleavage'}."\n";
	print PARFILE "protein_mass_filter = 0 0\n";
	print PARFILE "match_peak_count = 0\n";
	print PARFILE "match_peak_allowed_error = 1\n";
	print PARFILE "match_peak_tolerance = 1\n";
	print PARFILE "digest_mass_range = ".$$params{'min_peptide_mass'}." ".$$params{'max_peptide_mass'}."\n";
	print PARFILE "add_Cterm_peptide = ".$$params{'add_Cterm_peptide'}."\n";
	print PARFILE "add_Cterm_protein = 0.0000\n";
	print PARFILE "add_Nterm_peptide = ".$$params{'add_Nterm_peptide'}."\n";
	print PARFILE "add_Nterm_protein = 0.0000\n";
	print PARFILE "add_G_Glycine = ".$$params{'add_G_Glycine'}."\n";
	print PARFILE "add_A_Alanine = ".$$params{'add_A_Alanine'}."\n";
	print PARFILE "add_S_Serine = ".$$params{'add_S_Serine'}."\n";
	print PARFILE "add_P_Proline = ".$$params{'add_P_Proline'}."\n";
	print PARFILE "add_V_Valine = ".$$params{'add_V_Valine'}."\n";
	print PARFILE "add_T_Threonine = ".$$params{'add_T_Threonine'}."\n";
	print PARFILE "add_C_Cysteine = ".$$params{'add_C_Cysteine'}."\n";
	print PARFILE "add_L_Leucine = ".$$params{'add_L_Leucine'}."\n";
	print PARFILE "add_I_Isoleucine = ".$$params{'add_I_Isoleucine'}."\n";
	print PARFILE "add_X_LorI = ".$$params{'add_X_LorI'}."\n";
	print PARFILE "add_N_Asparagine = ".$$params{'add_N_Asparagine'}."\n";
	print PARFILE "add_O_Ornithine = ".$$params{'add_O_Ornithine'}."\n";
	print PARFILE "add_B_avg_NandD = ".$$params{'add_B_avg_NandD'}."\n";
	print PARFILE "add_D_Aspartic_Acid = ".$$params{'add_D_Aspartic_Acid'}."\n";
	print PARFILE "add_Q_Glutamine = ".$$params{'add_Q_Glutamine'}."\n";
	print PARFILE "add_K_Lysine = ".$$params{'add_K_Lysine'}."\n";
	print PARFILE "add_Z_avg_QandE = ".$$params{'add_Z_avg_QandE'}."\n";
	print PARFILE "add_E_Glutamic_Acid = ".$$params{'add_E_Glutamic_Acid'}."\n";
	print PARFILE "add_M_Methionine = ".$$params{'add_M_Methionine'}."\n";
	print PARFILE "add_H_Histidine = ".$$params{'add_H_Histidine'}."\n";
	print PARFILE "add_F_Phenylalanine = ".$$params{'add_F_Phenylalanine'}."\n";
	print PARFILE "add_R_Arginine = ".$$params{'add_R_Arginine'}."\n";
	print PARFILE "add_Y_Tyrosine = ".$$params{'add_Y_Tyrosine'}."\n";
	print PARFILE "add_W_Tryptophan = ".$$params{'add_W_Tryptophan'}."\n";
	print PARFILE "add_J_user_amino_acid = ".$$params{'add_J_user_amino_acid'}."\n";
	print PARFILE "add_U_user_amino_acid = ".$$params{'add_U_user_amino_acid'}."\n";
	print PARFILE "max_num_differential_AA_per_mod = ".$$params{'max_modif_num'}."\n";
	close (PARFILE);
}
1;
