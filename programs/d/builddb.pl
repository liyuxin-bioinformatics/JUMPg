#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/d

## Release date: 05/01/2015
## Release version: version 12.0.0

##################################################
## Builddb					##
## Generation of JUMP/SEQUEST search databases	##
## Generation of PIT (protein inference table)	##
##################################################

use strict;
use warnings;
use Cwd;
use File::Basename;
use Utils::Params;
use Utils::DatabaseUtils;
use Utils::DatabaseGeneration;

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}

open (my $log, ">", "jump_d.log");
my $paramsFile = shift (@ARGV);
unless (-e $paramsFile) {
	print "Cannot open $paramsFile\n";
	print "You need to have builddb.params file in the current working directory\n";
	print $log "Cannot open $paramsFile\n";
	print $log "You need to have builddb.params file in the current working directory\n";
	exit;
}
#print "JUMPd version11.1.1\n\n";
#print $log "JUMPd version11.1.1\n\n";

##################################################
## Read builddb.params file and Initialize	##
##################################################
my $params = Utils::Params -> new('-path' => $paramsFile);
$params = $params -> parse_param();
if (!-e $$params{'jump.params'}) {
	print "Cannot open $$params{'jump.params'}\n";
	print "Please check the path for \"jump.params\" in the builddb.params file carefully\n";
	print $log "Cannot open $$params{'jump.params'}\n";
	print $log "Please check the path for \"jump.params\" in the builddb.params file carefully\n";
	exit;
}
my $jumpParams = Utils::Params -> new('-path' => $$params{'jump.params'});
$jumpParams = $jumpParams -> parse_param();
my $utils = Utils::DatabaseUtils -> new();
my $db = Utils::DatabaseGeneration -> new();

##################################################
## Automatic generation of a new database name	##
##################################################
my $dbSuffix = $utils -> generateDbName($jumpParams, $log);
my $dbName = $$params{'output_prefix'}."_".$dbSuffix;
print "\n";
print $log "\n";

# check whether mda or hdr file already exist
if ($jumpParams->{search_engine} eq 'JUMP' and -e "$dbName.fasta.mdx"
or $jumpParams->{search_engine} eq 'SEQUEST' and -e "$dbName.fasta.hdr")
{ 
	print "mdx/hdr file already generated. skipped the database generation\n"; 
	exit;
}

##################################################################################
## Read input .fasta file(s) and 						##
## generate a new .fasta file and the corresponding .pit file for the database	##
## Contaminants and decoys will be added, if necessary				##
##################################################################################
$utils -> generateFastaAndPit($params, $dbName, $log);
print "\n";
print $log "\n";

##########################################################
## Generate a new JUMP or SEQUEST search database	##
##########################################################
if ($$params{'bypass_db_generation'} == 1) {
	print "Skipped the database generation (only .fasta and .pit files are generated)\n";
#} elsif (-e "$dbName.fasta.mdx") { print "mdx file already generated. skipped the database generation\n";
} else {
	$db -> generateDb($jumpParams, $dbName, $log);
	print "\n";
	print $log "\n";
}
close ($log);

# print $dbName to a temparoray file so that -s knows 
open(OUT,">\.jump_d_tmp");
my $currentpath=getcwd;
print OUT "$currentpath\/$dbName\.fasta\.mdx\n";
print OUT "$currentpath\/$dbName\.pit\n";
close OUT;
