#!/usr/bin/perl -I /data1/pipeline/release/version12.0.0/JUMPd/

## Release date: 05/01/2015
## Release version: version 12.0.0

use strict;
use warnings;
use Cwd;
use File::Basename;
use LWP::UserAgent;
use HTTP::Date;
use Utils::Params;
use Utils::DatabaseUtils;
use Utils::DatabaseGeneration;

my $utils = Utils::DatabaseUtils -> new();
my $db = Utils::DatabaseGeneration -> new();

##########################################################
## Make a directory to store generated databases	##
##########################################################

my ($sec, $min, $hour, $day, $mon, $year) = localtime;
$mon = sprintf("%02d", ($mon + 1));
$day = sprintf("%02d", $day);
$year = 1900 + $year;
my $dirName = "/home/jcho/testDBs/";
$dirName = $dirName.$year.$mon.$day;
system ("mkdir $dirName");
unless (-e $dirName) {
	print "$dirName does not exist or cannot be accessed\n";
	exit;
}
## Change the working directory to where the genrated databases are stored
chdir $dirName;

## Open a log file
open (my $log, ">", "updatedb.log");

##########################################################################
## Download proteome data (i.e. .fasta files) 				##
## from Uniprot repository					  	##
## Note that E.coli proteome should be downloaded through web-query	##
##########################################################################

#my @species = ("HUMAN", "MOUSE", "RAT", "DROME", "YEAST", "ARATH", "CAEEL", "DANRE", "ECOLI");
my @species = ("ECOLI");
my $baseurl = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/";
foreach my $species (@species) {
	if ($species eq "ECOLI") {
		print "Downloading ECOLI.fasta file from Uniprot repository (may take a while)\n";
		## E.coli proteome download through web-query
		my $contact = 'cjhjhj@gmail.com'; ## Set an email address here to help us debug in case of problems
		my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
		my $ecoliFile = "ECOLI.fasta";
		my $ecoliTaxon = "83333";
		my $ecoliKeyword = "keyword:181";	# complete proteome (keyword:1185 for reference proteome)
		my $ecoliQuery = "http://www.uniprot.org/uniprot/?query=organism:$ecoliTaxon+$ecoliKeyword&format=fasta&include=yes";
		my $ecoliResponse = $agent -> mirror($ecoliQuery, $ecoliFile);
		unless (-e $ecoliFile) {
			print "Problem in downloading ECOLI.fasta\n";
		}
		print "Finished downloading ECOLI.fasta\n";
	} else {
		my $fastaGz = $species.".fasta.gz";
		my $link = $baseurl.$fastaGz;
		system ("wget $link");
		unless (-e $fastaGz) {
			print "Problem in downloading $fastaGz\n";
			exit;
		}
		system ("gunzip $fastaGz");
		print "Finished downloading and unzipping $fastaGz\n";
	}
	my $fastaFile = $species.".fasta";
	
	## Generate a "jump_d" (builddb) parameter
	my $params = generateJumpdParams($fastaFile);
				
	## Generate two files
	## 1. a new .fasta file including contaminants and decoys
	## 2. a new .pit file
	my $speciesFileName = lc($species)."_contaminants_decoys";
	$utils -> generateFastaAndPit($params, $speciesFileName, $log);
	my $speciesFasta = $speciesFileName.".fasta";
	my $speciesPit = $speciesFileName.".pit";
	
	## Generate JUMP parameters for various conditions
	## Load a sample jump.params file
	my $paramsFile = "/data1/pipeline/release/version12.0.0/JUMPd/jump.params";
	my $jumpParams = Utils::Params -> new('-path' => $paramsFile);
	$jumpParams = $jumpParams -> parse_param();

	## Initialization
	my @searchEngines = ("SEQUEST", "JUMP");
	my @trypticConditions = ("ft", "pt");
	my @cysteineConditions = ("c0", "c57");
	my @tmtConditions = ("noTMT", "TMT");
	my @phoConditions = ("noPho", "pho");

	foreach my $engine (@searchEngines) {
		$$jumpParams{'search_engine'} = $engine;
		foreach my $tryptic (@trypticConditions) {
			foreach my $cysteine (@cysteineConditions) {
				foreach my $pho (@phoConditions) {
					foreach my $tmt (@tmtConditions) {
						my $dbSuffix = "";
						## Exceptions
						## 1. JUMP only considers the fully-tryptic condition
						## 2. SEQUEST does not consider the phosphorylation
						next if ($engine eq "JUMP" && $tryptic eq "pt");
						next if ($engine eq "SEQUEST" && $pho eq "pho");
						
						## Change JUMP parameters according to the conditions
						if ($tryptic eq "ft") {
							$$jumpParams{'enzyme_info'} = "Tryptic KR P";
							$$jumpParams{'digestion'} = "full";
							$$jumpParams{'max_mis_cleavage'} = "2";
							$dbSuffix = "ft_mc2_";
						} else {
							$$jumpParams{'enzyme_info'} = "Tryptic KR P";
							$$jumpParams{'digestion'} = "partial";
							$$jumpParams{'max_mis_cleavage'} = "5";
							$dbSuffix = "pt_mc5_";				
						}
						if ($cysteine eq "c0") {
							$$jumpParams{'add_C_Cysteine'} = "0.0000";
						} else {
							$$jumpParams{'add_C_Cysteine'} = "57.02146";
						}
						$dbSuffix = $dbSuffix.$cysteine;
						if ($pho eq "pho") {
							$$jumpParams{'dynamic_S'} = "79.96633";
							$$jumpParams{'dynamic_T'} = "79.96633";
							$$jumpParams{'dynamic_Y'} = "79.96633";
							$dbSuffix = $dbSuffix."_pho";				
						} else {
							if (defined $$jumpParams{'dynamic_S'}) {
								delete $$jumpParams{'dynamic_S'};
							}
							if (defined $$jumpParams{'dynamic_T'}) {
								delete $$jumpParams{'dynamic_T'};
							}
							if (defined $$jumpParams{'dynamic_Y'}) {
								delete $$jumpParams{'dynamic_Y'};
							}
						}
						if ($tmt eq "TMT") {
							$$jumpParams{'add_Nterm_peptide'} = "229.1629321";
							$$jumpParams{'add_K_Lysine'} = "229.1629321";
							$dbSuffix = $dbSuffix."_TMT_K229";
						} else {
							$$jumpParams{'add_Nterm_peptide'} = "0.0000";
							$$jumpParams{'add_K_Lysine'} = "0.0000";
						}
						my $targetDbName = lc($species)."_".$dbSuffix;
						my $targetFasta = $targetDbName.".fasta";
						my $targetPit = $targetDbName.".pit";
						system ("cp $speciesFasta $targetFasta");
						system ("cp $speciesPit $targetPit");
						$db -> generateDb($jumpParams, $targetDbName, $log);
					}
				}
			}
		}
	}
	system ("rm $speciesFasta $speciesPit");
}
close ($log);

sub generateJumpdParams {
	# When updating common databases, always include contaminants and decoys
	my ($inputFasta) = @_;
	my %p;
	$p{'include_contaminants'} = 1;
	$p{'input_contaminants'} = "/data1/database/contaminants.fasta";
	$p{'input_database1'} = getcwd()."/".basename($inputFasta);
	$p{'decoy_generation'} = 1;
	$p{'decoy_generation_method'} = 2;
	$p{'list_protein_abundance_database1'} = "/data1/database/KnowledgeTables/ProteinAbundance/Human_Abundance_emPAI.txt";
	$p{'list_protein_abundance_database2'} = "/data1/database/KnowledgeTables/ProteinAbundance/Mouse_Abundance_emPAI.txt";
	$p{'list_TFs'} = "/data1/database/KnowledgeTables/TFs/tfs_from_TRANSFAC.txt";
	$p{'list_oncogenes'} = "/data1/database/KnowledgeTables/Oncogenes/oncogenes_from_literatures.txt";
	$p{'list_kinases'} = "/data1/database/KnowledgeTables/Kinases/kinases_from_pkinfam.txt";
	$p{'list_GPCRs'} = "/data1/database/KnowledgeTables/GPCRs/gpcrs.txt";
	return (\%p);
}
