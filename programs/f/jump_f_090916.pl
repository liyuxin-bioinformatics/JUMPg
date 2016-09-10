#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.4/programs/f

## Release date: 11/01/2015
## Release version: version 12.1.0

use FindBin qw($Bin);
use lib "$Bin";
#print "$Bin\n";

use strict;
use idsum2::CommonUtils;
use idsum2::XMLParser;
use idsum2::IdsumHtml4_beta;
#use idsum2::GenerateHtml;
use idsum2::GenerateHtml_Tim07102014;
use Getopt::Long;
use Data::Dumper;
use Cwd;
use Storable;
use idsum2::CL;
use File::Basename;
use Statistics::Distributions;
use Clone qw(clone);
use FileHandle;
use idsum2::pepXML_parser;

## added by tim
use Sys::Hostname;
use Socket;

#################
#use File::Sync qw(fsync sync);

# idsum2 modeuls
use idsum2::ProteinGrouping;;
use idsum2::IdsumUtils6_beta;

#use idsum2_debug::ProteinGrouping;;
#use idsum2_debug::IdsumUtils6_beta;

my $gp = idsum2::ProteinGrouping->new;
my $idutils = idsum2::IdsumUtils6_beta->new();

my $debug=0;
my $xml = idsum2::XMLParser->new();
my $utils = idsum2::CommonUtils->new();
#my $idutils = Utils::IdsumUtils6_beta->new();
#my $html = Html::IdsumHtml4_beta->new();
#my $html = idsum2::GenerateHtml->new();
my $html = idsum2::GenerateHtml_Tim07102014->new();
my $pepxmlParser = idsum2::pepXML_parser->new();
my $annotationColHeader; # for pit annotation colomns
my $LARGE_GROUP_OFFSET=1000000; 	# for option 'prioritize_contaminants = 0'; add this number to contaminants
my $CONTAMINAT_PATTERN='co|CON_';	# for option 'prioritize_contaminants = 0'; match this pattern to decide whether a protein is contaminant

my ($Hydrogen_mass)=(1.007825032);


print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump filter                             ****     #
#       ****  Version 12.1.0                          ****     #
#       ****  Copyright (C) 2012 - 2015               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF

# For debug
#use idsum2_debug::IdsumUtils6_beta;
#my $idutils = idsum2_debug::IdsumUtils6_beta->new();

my $progname = $0;
$progname =~ s@(.*)/@@i;

$| = 1;
###### added by yanji
# parameter for output directory path
my $out_dir = "";
# check the number of groups of data
my $group_count = 1;
# log file name
#my $log_file = "idsum.log";
my $log_file = "jump_f.log";
open (LOGFILE, ">$log_file");
select((select(LOGFILE), $|=1)[0]);
#fsync(\*LOGFILE) or die "fsync: $!";sync();
#print LOGFILE "\n\n";
###### end of addition

# cgi paths
###################### changed by xusheng ###########################
#my $server = "proteox.genetics.emory.edu";
### changed by yanji
#my $server = "10.4.14.5";
#my $server = "spiders02.stjude.org";
#my $server = "10.4.15.169";

## Modified by Tim
#my $server = "spiderscluster.stjude.org";
my($addr)=inet_ntoa((gethostbyname(hostname))[4]);
#print "$addr\n";
my $server = $addr;
##########################################

###### end of change
#
################################################
my $showion = "http://$server/cgi-bin/displayions_html5";
#my $gp = Html::ProteinGrouping->new;

my %residhash; $utils->get_monoresidhash(\%residhash);

# Set options
my ($save_dir, $sum_dir); GetOptions("o=s" => \$save_dir, "s=s" => \$sum_dir);
$save_dir = "" if (!defined($save_dir)); $sum_dir = "sum_all" if (!defined($sum_dir));
my $paramsfile = shift @ARGV;
help() if (!defined($paramsfile));
if (!(-e $paramsfile)){
   print "$paramsfile does not exist\n";
   print LOGFILE "$paramsfile does not exist\n";
   while (!(-e $paramsfile)){ $paramsfile = idsum2::CL->ask("\nPlease enter the correct path to the parameter file.", ""); }
}

# Load parameter file
my (%paramhash, @runarray); $idutils->parse_params(\%paramhash, $paramsfile,\@runarray);
$idutils->parse_badproteinfile(\%{$paramhash{'badproteins'}}, "/usr/local/lib/perl5/Utils/TProtein.list");
if (!defined($paramhash{'one_hit_wonders_removal'})){
	print Dumper(\%paramhash); print "Please update your parameter file\n";
	print LOGFILE Dumper(\%paramhash);
	print LOGFILE "Please update your parameter file\n"; 
	exit;
}
$paramhash{'sorcererouts'} = 0;

# Make summary directories
my ($nomodsum_dir, $modsum_dir, $nomodsub_dir, $modsub_dir);
if (scalar(keys %{$paramhash{'groups'}}) > 1){
###### added by yanji
        $group_count = 2;
###### end of addition
	($nomodsum_dir, $modsum_dir) = (getcwd()."\/$sum_dir", getcwd()."\/$sum_dir\_mod");
	if (-e $nomodsum_dir) {system("rm -rf $nomodsum_dir");}
	$nomodsub_dir = $nomodsum_dir."\/fractions";	mkdir $nomodsum_dir; mkdir $nomodsub_dir; 
	if ($paramhash{'mods'} ne 0){	if (-e $modsum_dir){system("rm -rf $modsum_dir");}
	$modsub_dir = $modsum_dir."\/fractions"; mkdir $modsum_dir; mkdir $modsub_dir;}
}

# Run individual idsum on each group 
my $subgroup_percent = $paramhash{'subgroup'};
my (%sumhash, %modsumhash, $database, %dbhash,$pepXML_mode,%usedDBs);
for my $group (sort {$paramhash{'grouporder'}{$a}<=>$paramhash{'grouporder'}{$b}}keys %{$paramhash{'grouporder'}}){
undef($database);
my $members = $paramhash{'groups'}{$group};
#while (my ($group, $members) = each %{$paramhash{'groups'}}){
	my $group_dir = getcwd(); $group_dir = "$nomodsub_dir" if (defined($nomodsub_dir));
	$group_dir .= "\/sum\_$group"; $group_dir .= "\_$save_dir" if ($save_dir ne "");
	mkdir $group_dir;
	#print "$group_dir\n";

###### added by yanji
	$out_dir = $group_dir;
###### end of addition

	system("cp $paramsfile $group_dir");

#### create a dynamic link to jump
if (!(-e "$group_dir/html"))  { system("mkdir $group_dir/html"); }
open (DTAFILES, ">$group_dir/html/dtafiles");
for my $run (@$members){
    my @array = split(/\//, $run);
    my $len = @array;
    my $dtastuff = "";
    if ($len > 0) {
        $dtastuff = $array[$len - 1];
    }
    #for my $stuff (@array) {
    #    print DTAFILES "$stuff\t";
    #}
    print DTAFILES "$dtastuff\t$run\n";
}
close (DTAFILES);

	my (%peptidehash, %proteinhash, %runhash, %rundbhash);
	my (%peptidehash_1pct,%proteinhash_1pct);
	my $runnum = 0;
	my (%totalOutHash,%totalParaHash); # for all outfile in this group
	my ($ms2totalN,$outfileTotalN)=(0,0); # to count total # MS2 and outfiles

	for my $run (@$members){ # loop for each run

		# initialization
		my %del; $del{'contaminants'} = 0; $del{'length'} = 0; $del{'mixlabel'} = 0; $del{'maxmis'}=0; $del{'maxmod'}=0; $del{'DX'}=0; $del{'sharedecoy'}=0; $del{'no_required_mod'}=0; $del{'ppi_filtered'}=0;
		$del{'target'}{'contaminants'}=$del{'decoy'}{'contaminants'}=$del{'target'}{'length'}=$del{'decoy'}{'length'}=$del{'target'}{'mixlabel'}=$del{'decoy'}{'mixlabel'}=$del{'target'}{'maxmis'}=$del{'decoy'}{'maxmis'}=$del{'target'}{'maxmod'}=$del{'decoy'}{'maxmod'}=$del{'target'}{'DX'}=$del{'decoy'}{'DX'}=$del{'target'}{'trypz'}=$del{'decoy'}{'trypz'}=$del{'target'}{'no_required_mod'}=$del{'decoy'}{'no_required_mod'}=0;
    		my ($orignum, $count) = (0, 0);
		my $blank=0;
		my %ms2Scan;

		# parse jump.params
		$idutils->parse_jumpparams(\%paramhash, $run);
		if (!defined($paramhash{search_engine})) { die "Parameter search_engine not defined!!!\n"; }
		if (!defined($paramhash{pit_file})) { die "pit file not defined!!!\nPlease update your parameter file.\n"; }
		#print "search_engine: $paramhash{search_engine}\n";
		$paramhash{search_engine}=lc($paramhash{search_engine});

		# out file mode
		if ($paramhash{search_engine} eq 'sequest') {
		$idutils->parse_seqparams(\%paramhash, $run,$group_dir);

		print "\nWorking with $group:$run\n";
		print LOGFILE "\nWorking with $group:$run\n";
		print "search_engine: $paramhash{search_engine}\n";
		print "Counting dta and out files (please be patient)\n";
		print LOGFILE "Counting dta and out files (please be patient)\n";

		my (%paraHash, %outInforHash);
                my $pepxml = $run;
                if ($pepxml =~ /\/$/) { $pepxml =~ s/\/$/\.pepcXML/; }
                else  { $pepxml.='.pepcXML'; }
		#print "$pepxml\n";
                my $xml=$pepxml;
		if (-e $xml) { $pepXML_mode=1; } 
		else 
		{ 
			$xml =~ s/pepcXML$/pepXML/;
			if (-e $xml) { $pepXML_mode=1; } else { $pepXML_mode=0; }
		}
#$pepXML_mode=0;
		if ($pepXML_mode)
		{
                print "Reading pepXML ($xml)......\n";
                #$idutils->pepXML2Hashes(\%paraHash, \%outInforHash, $xml);
                $pepxmlParser->pepXML2Hashes(\%paraHash, \%outInforHash, $xml);
                my @out=keys(%outInforHash);
                $orignum=scalar(@out);
		print "There are $orignum total outfiles\n";
		print LOGFILE "There are $orignum total outfiles\n";

                for (my $i=0; $i<scalar(@out); $i++) {$out[$i].='.out';}

		# transfer outfiles from %outInforHash to %totalOutHash
		if (defined($paramhash{output_pepXML}) and $paramhash{output_pepXML})
		{
			foreach my $out (keys %outInforHash)
			{
				#$pepxmlParser->add_outfile(\%outInforHash,\%totalOutHash,$out);
				$pepxmlParser->add_outfile(\%paraHash, \%outInforHash,\%totalOutHash,$out);
			}
			%totalParaHash=%{clone(\%paraHash)};
		}		
		
		for my $outfile (@out){
			MS2_scan_count($outfile,\%ms2Scan);
                        $count++;
			#print "Parsing $outfile\r";
                        my ($db) = $idutils->pepXMLhashesTransfer_Sequest(\%paramhash, $run, $outfile, \%paraHash, \%outInforHash, \%peptidehash, \%runhash, \%del, \$blank);

                        next if (!defined($db));
			if (!(-e $db)){ $db = $paramhash{$run}{'database_name'}; }
			if (!defined($database)){
			        $database = $db;
		      } elsif ($database !~ /$db/){

		        	print "$database $db\n";
        			print "\nRuns in this group ($group) do not use the same database!!!\n";
		        	print "Please reorganized groups with only fractions that run with the same database\n";
		        	print LOGFILE "$database $db\n";
		        	print LOGFILE "\nRuns in this group ($group) do not use the same database!!!\n";
        			print LOGFILE "Please reorganized groups with only fractions that run with the same database\n";
			        exit;
      			}
		}
 
		} # end if (for pepXML mode rather than outfile mode)
		else
		{
		opendir(DIR,$run); my @out = grep {/\.out\Z/} readdir(DIR);
		seekdir (DIR,0); my @dta = grep {/\.dta\Z/} readdir(DIR);
		closedir(DIR);
    		$orignum=scalar(@out);
		printf LOGFILE "There are $orignum total outfiles.\n";
		if (scalar(@dta) != scalar(@out)){
			printf "THIS SEARCH IS MIGHT BE INCOMPLETE(%d dta files and %d out files)!!!!!\n", scalar(@dta), scalar(@out);
			printf LOGFILE "THIS SEARCH IS MIGHT BE INCOMPLETE(%d dta files and %d out files)!!!!!\n", scalar(@dta), scalar(@out);
			print "WOULD YOU LIKE TO CONTINUE? [NO]: ";
			print LOGFILE "WOULD YOU LIKE TO CONTINUE? [NO]: ";
  		chomp(my $choice = <STDIN>);	exit if ($choice !~ /^y/i);
		}
		if (scalar(@dta) == 0 || scalar(@out) == 0){ printf "     THERE ARE NO OUTPUT FILES FOR THIS RUN !!!!\n";exit;}

		# Load all outfile data
		for my $outfile (@out){
			MS2_scan_count($outfile,\%ms2Scan);

			$count++;
			print "\rGathering information from $count of $orignum outfiles ....     ";
			my ($db) = $idutils->parse_outfile(\%paramhash, \%peptidehash, \%runhash, $run, $outfile, \%del, \$blank);
			next if (!defined($db));
			if (!(-e $db)){ $db = $paramhash{$run}{'database_name'}; }
      			if (!defined($database)){
        			$database = $db;
		      	} elsif ($database !~ /$db/){
				
		        	print "$database $db\n";
	        		print "\nRuns in this group ($group) do not use the same database!!!\n";
		        	print "Please reorganized groups with only fractions that run with the same database\n";
	        		print LOGFILE "$database $db\n";
			        print LOGFILE "\nRuns in this group ($group) do not use the same database!!!\n";
        			print LOGFILE "Please reorganized groups with only fractions that run with the same database\n";
	        		exit;
	      		}
		}
		} # end else (for outfile mode rather than pepXML mode)
		} # end if ($paramhash{search_engine} eq 'sequest')
		# pepXML mode
		elsif ($paramhash{search_engine} eq 'jump') {
		$pepXML_mode=1;
		$idutils->parse_JUMPparams(\%paramhash, $run);   

		print "\nWorking with $group:$run\n";
		print LOGFILE "\nWorking with $group:$run\n";
		print "search_engine: $paramhash{search_engine}\n";
		opendir(DIR,$run); 		
		#my @out = grep {/\.out\Z/} readdir(DIR);
		#seekdir (DIR,0); my @dta = grep {/\.dta\Z/} readdir(DIR);
		#closedir(DIR);
	
		my (%paraHash, %outInforHash);
		my $pepxml = $run;
		if ($pepxml =~ /\/$/) { $pepxml =~ s/\/$/\.pepXML/; }
		else  { $pepxml.='.pepXML'; }
		my $xml=$pepxml;
		unless (-e $xml) { die "Not existed pepXML: $xml!!!\n"; }

		#my $xml="$run\/$pepxml[0]";
		#$idutils->pepXML2Hashes(\%paraHash, \%outInforHash, $xml);
		#print "testing!!!!!!!!!!\n";

		print "Reading pepXML ......\n";
		print LOGFILE "Reading pepXML ......\n";
		$pepxmlParser->pepXML2Hashes(\%paraHash, \%outInforHash, $xml);
		my @out=keys(%outInforHash);
    		$orignum=scalar(@out);
		print "There are $orignum total outfiles\n";
		print LOGFILE "There are $orignum total outfiles\n";
		$outfileTotalN+=$orignum;
		
		for (my $i=0; $i<scalar(@out); $i++) {$out[$i].='.spout';}

		# transfer outfiles from %outInforHash to %totalOutHash
		if (defined($paramhash{output_pepXML}) and $paramhash{output_pepXML})
		{
			foreach my $out (keys %outInforHash)
			{
				$pepxmlParser->add_outfile(\%paraHash, \%outInforHash,\%totalOutHash,$out);
			}
			%totalParaHash=%{clone(\%paraHash)};
		}


=head
		if (scalar(@dta) != scalar(@out)){
			printf "THIS SEARCH IS MIGHT BE INCOMPLETE(%d dta files and %d out files)!!!!!\n", scalar(@dta), scalar(@out);
			printf LOGFILE "THIS SEARCH IS MIGHT BE INCOMPLETE(%d dta files and %d out files)!!!!!\n", scalar(@dta), scalar(@out);
			print "WOULD YOU LIKE TO CONTINUE? [NO]: ";
			print LOGFILE "WOULD YOU LIKE TO CONTINUE? [NO]: ";
  		chomp(my $choice = <STDIN>);	exit if ($choice !~ /^y/i);
		}
=cut
		if (scalar(@out) == 0){ printf "     THERE ARE NO OUTPUT FILES FOR THIS RUN !!!!\n";exit;}

		for my $outfile (@out){
			MS2_scan_count($outfile,\%ms2Scan);

			$count++;
#			print "\rGathering information from $count of $orignum outfiles ....     ";
			my ($db) = $idutils->pepXMLhashesTransfer(\%paramhash, $run, $outfile, \%paraHash, \%outInforHash, \%peptidehash, \%runhash, \%del, \$blank);
			
			next if (!defined($db));
		
			if (!(-e $db)){ $db = $paramhash{$run}{'database_name'}; }
	  
			if (!defined($database)){
				$database = $db;
			}
		}
		} # end if ($paramhash{search_engine} eq 'jump')

		print "\n";
		print LOGFILE "\n";

		printf "There are %d total MS2 scans\n",scalar(keys %ms2Scan);
		printf LOGFILE "There are %d total MS2 scans\n",scalar(keys %ms2Scan);
		$ms2totalN+=scalar(keys %ms2Scan);

		#($runhash,$bottom,$roof,$step,$xcorr_dstr)
		my %xcorr_dstr;
		xcorr_TD_dstr($runhash{$run},$paramhash{min_XCorr},1.5,0.1,\%xcorr_dstr);


		my ($good, $bad) = (0,0);
# changed by xusheng because database file was changed ################### 
#		while (my ($outfile, $hash) = each %{$runhash{$run}}){ if ($$hash{'protein'} =~ /^(Random__|Decoy__)/){ $bad++; } else { $good++; }	}
		while (my ($outfile, $hash) = each %{$runhash{$run}}){ if ($$hash{'protein'} =~ /(Random__|Decoy__)/){ $bad++; } else { $good++; }     }
		if ($good == 0){ print "No outfiles have passed the initial minimal filtering.\n"; print LOGFILE "No outfiles have passed the initial minimal filtering.\n";exit; }
		print "Removing $blank outfiles containing no search results\n";
		print LOGFILE "Removing $blank outfiles containing no search results\n";
	if (defined($paramhash{ppi_included}))
	{
		print "Remove $del{ppi_filtered} scans due to ppi filtering (ppi included: $paramhash{ppi_included})\n";
		print LOGFILE "Remove $del{ppi_filtered} scans due to ppi filtering (ppi included: $paramhash{ppi_included})\n";
	}

		my $rmFDR=0;
		if (defined($paramhash{keep_only_mod_peptide})) { printf "Removing $del{no_required_mod} outfiles without required dynamic modifications ($paramhash{keep_only_mod_peptide}): $del{'target'}{'no_required_mod'} targets and $del{'decoy'}{'no_required_mod'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",100*$del{'decoy'}{'no_required_mod'}/(0.01+$del{target}{no_required_mod}); }
		$rmFDR=0;if ($del{'length'}) {$rmFDR=100*$del{'decoy'}{'length'}/(0.01+$del{'target'}{'length'});} #else {$rmFDR='na';} 
		printf "Removing $del{'length'} outfiles due to minimum peptide length ($paramhash{'min_peptide_length'}): $del{'target'}{'length'} targets and $del{'decoy'}{'length'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;
		$rmFDR=0;if ($del{'maxmis'}) {$rmFDR=100*$del{'decoy'}{'maxmis'}/(0.01+$del{target}{'maxmis'});#} else {$rmFDR='na';} 
		printf "Removing $del{'maxmis'} outfiles due to max miscleavage sites ($paramhash{'max_peptide_mis'}): $del{'target'}{'maxmis'} targets and $del{'decoy'}{'maxmis'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;}
		$rmFDR=0;if ($del{'maxmod'}) {$rmFDR=100*$del{'decoy'}{'maxmod'}/(0.01+$del{target}{'maxmod'});} #else {$rmFDR='na';} 
		printf "Removing $del{'maxmod'} outfiles due to max modification sites ($paramhash{'max_peptide_mod'}): $del{'target'}{'maxmod'} targets and $del{'decoy'}{'maxmod'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;
		$rmFDR=0;if ($del{'DX'}) {$rmFDR=100*$del{'decoy'}{'DX'}/(0.01+$del{target}{'DX'});} #else {$rmFDR='na';}
		printf "Removing $del{'DX'} outfiles due to minimum XCorr ($paramhash{'min_XCorr'}) and dCn values ($paramhash{'min_dCn'}): $del{'target'}{'DX'} targets and $del{'decoy'}{'DX'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR; 
		for (my $i=0;$i<(1.5-$paramhash{min_XCorr})/0.1; $i++) 
		{ 
			printf "  XCorr threshold %.2f: %d targets and %d decoys deleted (FDR for deleted outfiles = %.2f%%)\n",
				$paramhash{min_XCorr}+($i+1)*0.1,
				$xcorr_dstr{T}[$i]+$del{'target'}{'DX'},
				$xcorr_dstr{D}[$i]+$del{'decoy'}{'DX'},
				#200*($del{'decoy'}{'DX'}+$xcorr_dstr{D}[$i])/($del{'DX'}+$xcorr_dstr{T}[$i]+$xcorr_dstr{D}[$i]+0.01); 
				100*($del{'decoy'}{'DX'}+$xcorr_dstr{D}[$i])/($del{target}{'DX'}+$xcorr_dstr{T}[$i] + 0.01); 
		} 
		$rmFDR=0;if ($del{'modremoval'}) {$rmFDR=100*$del{'decoy'}{'modremoval'}/($del{target}{'modremoval'}+0.01); #else {$rmFDR='na';}
		printf "Removing $del{'modremoval'} outfiles due to peptide_mod_removal ($paramhash{'peptide_mod_removal'}): $del{'target'}{'modremoval'} targets and $del{'decoy'}{'modremoval'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR; }
		$rmFDR=0;if ($del{'aaremoval'}) {$rmFDR=100*$del{'decoy'}{'aaremoval'}/($del{target}{'aaremoval'}+0.01); #else {$rmFDR='na';}
		printf "Removing $del{'aaremoval'} outfiles due to peptide_aa_removal ($paramhash{'peptide_aa_removal'}): $del{'target'}{'aaremoval'} targets and $del{'decoy'}{'aaremoval'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR; }
		$rmFDR=0;if ($del{'mix_label'}) {$rmFDR=100*$del{'decoy'}{'mix_label'}/$del{target}{'mix_label'};} #else {$rmFDR='na';} 
		if ($paramhash{'mix_label'} ne 0) {printf "Mixed Peptide Filtering for ($paramhash{'mix_label'}) removed $del{'mixlabel'} outfiles: $del{'target'}{'mixlabel'} targets and $del{'decoy'}{'mixlabel'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;}
#		print "There are $del{'sharedecoy'} outfiles matched to both target and decoy proteins (remove the decoy proteins)\n";
		print "Contaminant Peptide Filtering removed $del{'contaminants'} outfiles\n" if ($paramhash{'filter_contaminants'} ne 0);
		#printf "There are %d outfiles: $good targets and $bad decoys(FDR for filtered outfiles = %.2f%%)\n", $good+$bad,200*$bad/($good+$bad);
		#

		# print in logfile
		my $rmFDR=0;
		if (defined($paramhash{keep_only_mod_peptide})) { printf  LOGFILE "Removing $del{no_required_mod} outfiles without required dynamic modifications ($paramhash{keep_only_mod_peptide}): $del{'target'}{'no_required_mod'} targets and $del{'decoy'}{'no_required_mod'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",100*$del{'decoy'}{'no_required_mod'}/(0.01+$del{target}{no_required_mod}); }
		$rmFDR=0;if ($del{'length'}) {$rmFDR=100*$del{'decoy'}{'length'}/(0.01+$del{'target'}{'length'});} #else {$rmFDR='na';} 
		printf  LOGFILE "Removing $del{'length'} outfiles due to minimum peptide length ($paramhash{'min_peptide_length'}): $del{'target'}{'length'} targets and $del{'decoy'}{'length'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;
		$rmFDR=0;if ($del{'maxmis'}) {$rmFDR=100*$del{'decoy'}{'maxmis'}/(0.01+$del{target}{'maxmis'});#} else {$rmFDR='na';} 
		printf  LOGFILE "Removing $del{'maxmis'} outfiles due to max miscleavage sites ($paramhash{'max_peptide_mis'}): $del{'target'}{'maxmis'} targets and $del{'decoy'}{'maxmis'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;}
		$rmFDR=0;if ($del{'maxmod'}) {$rmFDR=100*$del{'decoy'}{'maxmod'}/(0.01+$del{target}{'maxmod'});} #else {$rmFDR='na';} 
		printf  LOGFILE "Removing $del{'maxmod'} outfiles due to max modification sites ($paramhash{'max_peptide_mod'}): $del{'target'}{'maxmod'} targets and $del{'decoy'}{'maxmod'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;
		$rmFDR=0;if ($del{'DX'}) {$rmFDR=100*$del{'decoy'}{'DX'}/(0.01+$del{target}{'DX'});} #else {$rmFDR='na';}
		printf  LOGFILE "Removing $del{'DX'} outfiles due to minimum XCorr ($paramhash{'min_XCorr'}) and dCn values ($paramhash{'min_dCn'}): $del{'target'}{'DX'} targets and $del{'decoy'}{'DX'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR; 
		for (my $i=0;$i<(1.5-$paramhash{min_XCorr})/0.1; $i++) 
		{ 
			printf  LOGFILE "  XCorr threshold %.2f: %d targets and %d decoys deleted (FDR for deleted outfiles = %.2f%%)\n",
				$paramhash{min_XCorr}+($i+1)*0.1,
				$xcorr_dstr{T}[$i]+$del{'target'}{'DX'},
				$xcorr_dstr{D}[$i]+$del{'decoy'}{'DX'},
				#200*($del{'decoy'}{'DX'}+$xcorr_dstr{D}[$i])/($del{'DX'}+$xcorr_dstr{T}[$i]+$xcorr_dstr{D}[$i]+0.01); 
				100*($del{'decoy'}{'DX'}+$xcorr_dstr{D}[$i])/($del{target}{'DX'}+$xcorr_dstr{T}[$i] + 0.01); 
		} 
		$rmFDR=0;if ($del{'modremoval'}) {$rmFDR=100*$del{'decoy'}{'modremoval'}/($del{target}{'modremoval'}+0.01); #else {$rmFDR='na';}
		printf  LOGFILE "Removing $del{'modremoval'} outfiles due to peptide_mod_removal ($paramhash{'peptide_mod_removal'}): $del{'target'}{'modremoval'} targets and $del{'decoy'}{'modremoval'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR; }
		$rmFDR=0;if ($del{'aaremoval'}) {$rmFDR=100*$del{'decoy'}{'aaremoval'}/($del{target}{'aaremoval'}+0.01); #else {$rmFDR='na';}
		printf  LOGFILE "Removing $del{'aaremoval'} outfiles due to peptide_aa_removal ($paramhash{'peptide_aa_removal'}): $del{'target'}{'aaremoval'} targets and $del{'decoy'}{'aaremoval'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR; }
		$rmFDR=0;if ($del{'mix_label'}) {$rmFDR=100*$del{'decoy'}{'mix_label'}/$del{target}{'mix_label'};} #else {$rmFDR='na';} 
		if ($paramhash{'mix_label'} ne 0) {printf  LOGFILE "Mixed Peptide Filtering for ($paramhash{'mix_label'}) removed $del{'mixlabel'} outfiles: $del{'target'}{'mixlabel'} targets and $del{'decoy'}{'mixlabel'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;}
		print LOGFILE "Contaminant Peptide Filtering removed $del{'contaminants'} outfiles\n" if ($paramhash{'filter_contaminants'} ne 0);
#=head
		#MS2 consolidation
		#print "$runhash{$run}{'u170k_ctr_02ug.32227.1.3.out'}{protein}\n";
		rm_multiple_charge_for_single_precursor($runhash{$run},\%peptidehash,\%del);
		printf "Removing %d outfiles due to \+2\/\+3 charge assignments of uncharged precursor ions: %d targets and %d decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$del{multiCharge}{T}+$del{multiCharge}{D},$del{multiCharge}{T},$del{multiCharge}{D},100*$del{multiCharge}{D}/($del{multiCharge}{T} + 0.01);
		printf LOGFILE "Removing %d outfiles due to \+2\/\+3 charge assignments of uncharged precursor ions: %d targets and %d decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$del{multiCharge}{T}+$del{multiCharge}{D},$del{multiCharge}{T},$del{multiCharge}{D},100*$del{multiCharge}{D}/($del{multiCharge}{T} + 0.01);

		rm_multi_monoisotopic_outfiles($runhash{$run},\%peptidehash,\%del);
		printf "Removing %d outfiles due to multiple monoisotoptic assigments of the same precursur ion: %d targets and %d decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$del{multiMonoisot}{T}+$del{multiMonoisot}{D},$del{multiMonoisot}{T},$del{multiMonoisot}{D},100*$del{multiMonoisot}{D}/($del{multiMonoisot}{T} + 0.01);
		printf LOGFILE "Removing %d outfiles due to multiple monoisotoptic assigments of the same precursur ion: %d targets and %d decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$del{multiMonoisot}{T}+$del{multiMonoisot}{D},$del{multiMonoisot}{T},$del{multiMonoisot}{D},100*$del{multiMonoisot}{D}/($del{multiMonoisot}{T} + 0.01);
#=cut
		# tryptic and charge filtering
		($good, $bad) = (0,0); my %trypcharge;
		for my $combo (keys %{$paramhash{'12combinations'}}){
			next if ($paramhash{'12combinations'}{$combo} == 0);#print "$combo\n";
			$combo =~ s/([A-Z]+)(\d)//; my ($tryptic, $charge) = ($1, $2);
			if ($tryptic =~ /FT/){ $trypcharge{'2'}{$charge} = 1;
			} elsif ($tryptic =~ /PT/){ $trypcharge{'1'}{$charge} = 1;
			} elsif ($tryptic =~ /NT/){ $trypcharge{'0'}{$charge} = 1;
			}
		}
		while (my ($outfile, $hash) = each %{$runhash{$run}}){
			my ($intpep, $tryptic, $charge,$protein) = ($$hash{'intpep'}, $$hash{'tryptic'}, $$hash{'charge'},$$hash{'protein'});
			#print "$outfile:$$hash{'peptide'}, $intpep, $tryptic, $charge,$protein\n";
			if (!defined($trypcharge{$tryptic})){
				#print "Not defined trypcharge (1): $outfile,$tryptic,$charge,$$hash{peptide}\n";
				#print LOGFILE "Not defined trypcharge (1): $outfile,$tryptic,$charge,$$hash{peptide}\n";
				delete $runhash{$run}{$outfile}; delete ($peptidehash{$intpep}{'outfiles'}{$outfile}); $del{'trypz'}++;
				if ($protein =~ /Decoy/) { $del{'decoy'}{'trypz'}++; } else { $del{'target'}{'trypz'}++; }
			} else {
				if (!defined($trypcharge{$tryptic}{$charge})){
				#print "Not defined trypcharge (2): $outfile,$tryptic,$charge,$$hash{peptide}\n";
				#print LOGFILE "Not defined trypcharge (2): $outfile,$tryptic,$charge,$$hash{peptide}\n";
					delete $runhash{$run}{$outfile}; delete ($peptidehash{$intpep}{'outfiles'}{$outfile}); $del{'trypz'}++;
				if ($protein =~ /Decoy/) { $del{'decoy'}{'trypz'}++; } else { $del{'target'}{'trypz'}++; }
				}
			}
		}
		for my $intpep (keys %peptidehash){
			if (scalar(keys %{$peptidehash{$intpep}{'outfiles'}}) == 0){
				delete ($peptidehash{$intpep});
			}
		}
		($good, $bad) = (0,0);
		while (my ($outfile, $hash) = each %{$runhash{$run}}){
###################### changed by xusheng ######################
#			if ($$hash{'protein'} =~ /^Random__/ || $$hash{'protein'} =~ /^Decoy__/){ $bad++;
			if ($$hash{'protein'} =~ /Random__/ || $$hash{'protein'} =~ /Decoy__/){ $bad++;
			} else { $good++; }
		}
		if ($good == 0){
			print "No outfiles have passed the initial minimal filtering.\n"; 
			print LOGFILE "No outfiles have passed the initial minimal filtering.\n";exit;
		}
		if (!defined($del{'trypz'})) { $del{'trypz'}=0; }
		if (defined($del{'trypz'})){
			print "Removing $del{'trypz'} outfiles due to trypticity and charge (";
			print LOGFILE "Removing $del{'trypz'} outfiles due to trypticity and charge (";
			my $str = "";
			for my $tryptic (keys %trypcharge){
				if ($tryptic == 2){	$str .= "FT";	} elsif ($tryptic == 1){ $str .= "PT"; } else {	$str .= "NT";	}
				for my $charge (sort keys %{$trypcharge{$tryptic}}){
					$str .= "$charge";
			}
				$str .= " ";
			}
			$str =~ s/\s\Z//;
			if ($del{'trypz'}) {$rmFDR=100*$del{'decoy'}{'trypz'}/($del{target}{'trypz'}+0.1);} else {$rmFDR='na';} 
			if (defined($del{'trypz'})) {printf "$str): $del{'target'}{'trypz'} targets and $del{'decoy'}{'trypz'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;}
			if (defined($del{'trypz'})) {printf LOGFILE "$str): $del{'target'}{'trypz'} targets and $del{'decoy'}{'trypz'} decoys deleted (FDR for deleted outfiles = %.2f%%)\n",$rmFDR;}
			#printf "Accepting %d (out of $orignum, %.2f%%) outfiles: $good targets and $bad decoys (FDR for filtered outfiles = %.2f%%, expecting %d true targets)\n", $good+$bad,($good+$bad)*100/$orignum ,200*$bad/($good+$bad),$good-$bad;
			printf "Accepting %d (out of $orignum, %.2f%%) outfiles: $good targets and $bad decoys (FDR for filtered outfiles = %.2f%%, expecting %d true targets)\n", $good+$bad,($good+$bad)*100/$orignum ,100*$bad/($good+0.01),$good-$bad;
			printf LOGFILE "Accepting %d (out of $orignum, %.2f%%) outfiles: $good targets and $bad decoys (FDR for filtered outfiles = %.2f%%, expecting %d true targets)\n", $good+$bad,($good+$bad)*100/$orignum ,100*$bad/($good+0.01),$good-$bad;
			#print LOGFILE "$str) filtering removed $del{'trypz'} outfiles\n" if (defined($del{'trypz'}));;
                        #printf LOGFILE "There are %d outfiles: $good from target and $bad from decoy database\n", $good+$bad, $good, $bad;
		}

	}
	if ($paramhash{'combine_filtermethods'}){
		eval { store(\%runhash, "$group_dir/preMA_runhash"); };
 		print "Error writing to file: $@\n" if $@;
		print LOGFILE "Error writing to file: $@\n" if $@;
		eval { store(\%peptidehash, "$group_dir/preMA_peptidehash"); };
 		print "Error writing to file: $@\n" if $@;
		print LOGFILE "Error writing to file: $@\n" if $@;
	}	

	
	# Bypass Filtering 
	if (defined($paramhash{database})) { $database=$paramhash{database}; }
	if ($database =~ /\.hdr$/) { $database =~ s/\.hdr$//; }
	if ($database =~ /\.mdx$/) { $database =~ s/\.mdx$//; }
	#print "Creating databse hash ($database) ...\n";
	my %fprhash;
	my ($numofproteins);
	if (!defined($usedDBs{$database}))
	{
		($numofproteins) = $utils->create_dbHash($database,\%dbhash);
		$usedDBs{$database}=1;
	}
	#%dbhash = %$tempdbhash;
	if ($paramhash{'bypass_filtering'}){
		BypassFiltering(\%peptidehash, \%proteinhash, \%fprhash, \%runhash, $group);
	} else {
		# Mass Accuracy Filtering
		if ($paramhash{'mass_accuracy'}){
			print "Running Mass Accuracy Filtering for $group using only $paramhash{'mass_consideration'}{'type'}\n";
			print LOGFILE "\nRunning Mass Accuracy Filtering for $group using only $paramhash{'mass_consideration'}{'type'}\n";
			for my $run (sort keys %runhash){
=head
				my $mzxml = $run; 
				if ($mzxml =~ /\.\d+\Z/) {$mzxml =~ s/\.[\d]+\Z/.mzXML/;}
				elsif ($mzxml =~ /\.\d+\/\Z/) {$mzxml =~ s/\.[\d]+\/\Z/.mzXML/;}
				elsif ($mzxml =~ /\/$/) {$mzxml =~ s/\/$/.mzXML/;}
				else { $mzxml .= ".mzXML";} #DMD June 28, 2007 next 9 lines
				if (!(-e $mzxml)){ print "$mzxml does not exist!!!!\n"; print LOGFILE "$mzxml does not exist!!!!\n"; exit;}
=cut
				print "";
				my $source = $paramhash{'MS_Source'};
				if ($source == 0){
					#print "Survey scans in this run come from LTQ, skipping mass accuracy\n";
					#print LOGFILE "Survey scans in this run come from LTQ, skipping mass accuracy\n";
					for my $outfile (keys %{$runhash{$run}}){	$runhash{$run}{$outfile}{'skipma'} = 1;	}
					next;
				} else {
					#print "Survey scans in this run come from ORBI\n";
					#print LOGFILE "Survey scans in this run come from ORBI\n";
				}

				print "\nRunning Mass Accuracy Filtering for $run\n";
				print LOGFILE "\nRunning Mass Accuracy Filtering for $run\n";
				close LOGFILE;

				$idutils->acumass($run, $runhash{$run}, \%paramhash);
				my $total = scalar(keys %{$runhash{$run}});
				my ($good, $bad, $grandom, $brandom) = (0,0,0,0);
	#mkdir "randomoutdir";
				open (LOGFILE, ">>$log_file");
				while (my ($outfile, $hash) = each %{$runhash{$run}}){
					my $intpep = $$hash{'intpep'};
					#if ($$hash{'protein'} =~ /^Decoy__/){
					#	system ("cp $$hash{'path'} randomoutdir/.");
					#}	
					if ($$hash{'status'} == 1){
########## changed by xusheng 
#						$good++; $grandom++ if ($$hash{'protein'}=~/^Random__/ || $$hash{'protein'} =~ /^Decoy__/);
						$good++; $grandom++ if ($$hash{'protein'}=~/Random__/ || $$hash{'protein'} =~ /Decoy__/);
						$peptidehash{$intpep}{'outfiles'}{$outfile}{'ppm'} = $$hash{'correctedppm'};
						$peptidehash{$intpep}{'outfiles'}{$outfile}{'rawppm'} = $$hash{'ppm'};
					} else {
########### changed by xusheng ###################
#						$bad++; $brandom++ if ($$hash{'protein'}=~/^Random__/ || $$hash{'protein'} =~ /^Decoy__/);
						$bad++; $brandom++ if ($$hash{'protein'}=~/Random__/ || $$hash{'protein'} =~ /Decoy__/);						
						delete $peptidehash{$intpep}{'outfiles'}{$outfile};
					}
				}
				if ($bad)
				{
					#printf "Removing $bad outfiles due to precursor ion mass accuracy: %d targets and $brandom decoys (FDR for deleted outfiles = %.2f%%)\n",$bad-$brandom,200*$brandom/($bad+0.01);
					printf "Removing $bad outfiles due to precursor ion mass accuracy: %d targets and $brandom decoys (FDR for deleted outfiles = %.2f%%)\n",$bad-$brandom,100*$brandom/($bad-$brandom+0.01);
					printf LOGFILE "Removing $bad outfiles due to precursor ion mass accuracy: %d targets and $brandom decoys (FDR for deleted outfiles = %.2f%%)\n",$bad-$brandom,100*$brandom/($bad-$brandom+0.01);
				}
				else 
				{
					print "Removing no outfiles due to precursor ion mass accuracy\n";
				}
				if ($good)
				{
					#printf "Accepting $good (out of $total, %.2f%%) outfiles: %d targets and $grandom decoys (FDR for filtered outfiles = %.2f%%, expecting %d true targets)\n",$good*100/$total,$good-$grandom,200*$grandom/($good+0.01),$good-2*$grandom;
					printf "Accepting $good (out of $total, %.2f%%) outfiles: %d targets and $grandom decoys (FDR for filtered outfiles = %.2f%%, expecting %d true targets)\n",$good*100/$total,$good-$grandom,100*$grandom/($good-$grandom+0.01),$good-2*$grandom;
					printf LOGFILE "Accepting $good (out of $total, %.2f%%) outfiles: %d targets and $grandom decoys (FDR for filtered outfiles = %.2f%%, expecting %d true targets)\n",$good*100/$total,$good-$grandom,100*$grandom/($good-$grandom+0.01),$good-2*$grandom;
				}
				if ( $total-2*($grandom+$brandom) )
				{
					print "  Recovery rate of expected true PSMs = (",$good-$grandom,"-$grandom)/(",$total-$grandom-$brandom,"-",$grandom+$brandom,") = ";
					printf "%.2f%%\n",100*($good-2*$grandom)/($total-2*($grandom+$brandom));
					print LOGFILE "  Recovery rate of expected true PSMs = (",$good-$grandom,"-$grandom)/(",$total-$grandom-$brandom,"-",$grandom+$brandom,") = ";
					printf LOGFILE "%.2f%%\n",100*($good-2*$grandom)/($total-2*($grandom+$brandom));
					#printf "  Recovery rate of expected true PSMs = ($good-2*$grandom)/($total-2*%d) = %.2f%%\n",$grandom+$brandom,100*($good-2*$grandom)/($total-2*($grandom+$brandom));
					#printf LOGFILE "  Recovery rate of expected true PSMs = ($good-2*$grandom)/($total-2*%d) = %.2f%%\n",$grandom+$brandom,100*($good-2*$grandom)/($total-2*($grandom+$brandom));
				}
			}

			while (my ($peptide, $pephash) = each %peptidehash){	
				if (scalar(keys %{$$pephash{'outfiles'}}) == 0){
					delete $peptidehash{$peptide};
				}
			}
		}
		## Add in redundancy for SORCERER outfiles
		if ($paramhash{'add_redundancy'} == 1 && $paramhash{'sorcererouts'} == 1){
			my %redpephash;
			while (my ($peptide, $pephash) = each %peptidehash){
				#next if (scalar(keys %{$$pephash{'proteins'}}) > 1); DMD Feb 5 2009
				my $nomod = $peptide; $nomod =~ s/[\#\*\@\%\&\~\$\^]//g;
				$redpephash{$nomod}{'peps'}{$peptide} = 1;
			}
			my %reddbhash; my @redseqarray;
			while (my ($protein, $prohash) = each %dbhash){
######## changed by xusheng 
#				next if ($protein =~ /^Decoy/ || $protein =~ /^Random/);
				next if ($protein =~ /Decoy/ || $protein =~ /Random/);
				$reddbhash{$protein}{'annotation'} = substr $$prohash{'annotation'}, 0, 73-length($protein);
				$reddbhash{$protein}{'sequence'} = $$prohash{'sequence'};
				push (@redseqarray, $$prohash{'sequence'});
			}
			my ($start, $pepnum) = (1, scalar(keys %redpephash));
			for my $peptide (keys %redpephash){
				print "\rAdding redundancy for $start of $pepnum peptides ...";
				#my @keys = grep {$reddbhash{$_}{'sequence'} =~ m/$peptide/} keys %reddbhash;
				my @keys = grep(/$peptide/, @redseqarray);
				$start++;
				next if (scalar (@keys) < 2);
				@keys = grep {$reddbhash{$_}{'sequence'} =~ m/$peptide/} keys %reddbhash;
				for my $pep (keys %{$redpephash{$peptide}{'peps'}}){
					my %tempouthash = %{$peptidehash{$pep}{'outfiles'}};
					#print Dumper(\%tempouthash);exit;
					for my $protein (@keys){
						$peptidehash{$pep}{'proteins'}{$protein}++;
					}
					for my $file (keys %tempouthash){
						my $outfile = $tempouthash{$file}{'path'}; 
						open (IN, "+<$outfile");
						my @lines = <IN>;
						seek IN,0,0;
						my @first15 = splice(@lines, 0,15);
						for (@first15){ print IN "$_" };
						my $temp = shift @lines;
						next if ($temp =~ /\+\d+/);
						my @array = split(/([A-Z\-]\.)/, $temp, 2);
						chop($array[0]);chop($array[0]);
						printf IN " $array[0] +%d  $array[1]$array[2]", scalar(@keys)-1;
						for my $protein (@keys){
							next if ($tempouthash{$file}{'protein'} eq $protein);
							print IN "      $protein $reddbhash{$protein}{'annotation'}\n";
						}
						for (@lines){ print IN "$_" };
					}
				}
			}
		}

		#protein grouping before filtering
		#$idutils->create_proteinhash(\%proteinhash, \%peptidehash, \%dbhash);	count_nomod(\%proteinhash);
		#my ($grouphash, $group_num, $subgroup_num) = $gp->group_proteins(\%proteinhash, \%peptidehash, \%paramhash);

		# print feature hash
		my (%featurehash);
		buildFeatureHash(\%peptidehash, \%featurehash);

		# check if featurehash is empty
		if ( scalar(keys %featurehash)==0 ) { die "No scans accepted before FDR filtering!!!\n Please check your parameter and data!!!\n"; }

		#open(FTR,">feature_infor.txt");
		if (!(-e "$group_dir\/misc")) { system("mkdir $group_dir\/misc"); }
		open(FTR,">$group_dir\/misc\/feature_infor.txt");
		print FTR "outfile\txcorr\tdCn\tpeptideLength\ttryptic\tppm\tmod\tmc\tcharge1\tcharge2\tcharge3\tcharge4\ttargetDecoy\n";
		foreach my $outfile (keys %featurehash)
		{
			print FTR "$outfile\t$featurehash{$outfile}{XCorr}\t$featurehash{$outfile}{dCn}\t$featurehash{$outfile}{peptideLength}\t$featurehash{$outfile}{tryptic}\t$featurehash{$outfile}{ppm}\t$featurehash{$outfile}{mod}\t$featurehash{$outfile}{mis}\t";
			for (my $i=1; $i<=4; $i++)
			{
				if ( $featurehash{$outfile}{charge} == $i ) { print FTR "1\t"; }
				else { print FTR "0\t"; }
			}
			if ($featurehash{$outfile}{random}==1) { print FTR "decoy"; }
			else { print FTR "target"; }
			print FTR "\n";
		}
		close FTR;

		# for each peptide, assign 'best' protein
		my %best_proteins;
		#if ($use_pit)
		#{
		#assignBestProtein2Peptide(\%peptidehash, $grouphash, $decoy_grouphash, \%best_proteins);
		#assignOrder2group($grouphash, \%best_proteins);
		#assignOrder2group($decoy_grouphash, \%best_proteins);
		#}

		# LDA analysis: %featurehash{$outfile}{qvalue}
		my (%featurehash); my $LDAfile="$group_dir\/misc\/feature_infor.txt";
		#if (defined($paramhash{'FDR_filtering_method'}) and $paramhash{'FDR_filtering_method'} eq 'LDA')
		if (1)
		{
			print "\nPerforming LDA analysis (~300k outfiles per minute)\n";
			print LOGFILE "\nPerforming LDA analysis (~300k outfiles per minute)\n";
			$idutils->buildFeatureHash(\%peptidehash, \%featurehash);
			$idutils->printFeatureHash($LDAfile, \%featurehash);
			my $ldaresult=$idutils->rLDA_analysis($LDAfile,$paramhash{'search_engine'},"$Bin/idsum2/LDA.R");
			$idutils->attchLDAresult2FeatureHash($LDAfile, \%featurehash);

			# peptide FDR
			$idutils->peptide_Qvalue(\%peptidehash, \%featurehash);
		}

		# protein FDR filtering:
		print "\nRunning score filtering:\n";
		print LOGFILE "\nRunning score filtering:\n";
		my ($org_good,$org_bad)=scan_TD_num(\%peptidehash);
		printf "\nStarting with %d outfiles: $org_good targets and $org_bad decoys (FDR = %.2f%%)\n",$org_good+$org_bad,100*$org_bad/($org_good+0.01);
		printf LOGFILE"\nStarting with %d outfiles: $org_good targets and $org_bad decoys (FDR = %.2f%%)\n",$org_good+$org_bad,100*$org_bad/($org_good+0.01);

		# protein FDR filtering based on 3 groups
		# cut scans at a specific initial FDR (e.g., 1% or 5%)
		print "\nInitial filtering: cut outfiles at a specific initial FDR ($paramhash{initial_outfile_fdr}%):\n";
		print LOGFILE"\nInitial filtering: cut outfiles at a specific initial FDR ($paramhash{initial_outfile_fdr}%):\n";
		my ($pepfpr,$peptide_decoy,$peptide_target,$adjust_pepfpr) = $idutils->find_XCorrdCn_cutoffs(\%peptidehash, \%paramhash, $paramhash{initial_outfile_fdr}, \%featurehash, 0, 1);
		$idutils->create_proteinhash(\%proteinhash, \%peptidehash, \%dbhash);
		my ($current_good,$current_bad)=scan_TD_num(\%peptidehash);
		printf "  Recovery rate of expected true PSMs after initial filtering: ($current_good-$current_bad)/($org_good-$org_bad) = %.2f%%\n",($current_good-$current_bad)*100/($org_good-$org_bad+0.01);
		printf LOGFILE "  Recovery rate of expected true PSMs after initial filtering: ($current_good-$current_bad)/($org_good-$org_bad) = %.2f%%\n",($current_good-$current_bad)*100/($org_good-$org_bad+0.01);

		if (!defined($paramhash{multistep_FDR_filtering}) or $paramhash{multistep_FDR_filtering}==1 ) # whether perform mupltistep_FDR_filtering
		{
		# clone hashes for protein q value estimation
		my %peptidehash1 = %{clone(\%peptidehash)};
		my %runhash1 = %{clone(\%runhash)};
		my %featurehash1 = %{clone(\%featurehash)};
		my %fprhash1 = %{clone(\%fprhash)};
		my %proteinhash1 = %{clone(\%proteinhash)};

		if ( $paramhash{unique_protein_or_peptide} eq 'protein' )
		{
		# split %peptidehash
		#print "peptidehash: ",scalar(keys %peptidehash),"\n";
		print "\nSplitting scans into 3 categories according to protein assignment (please be patient)\n";
		print LOGFILE"\nSplitting scans into 3 categories according to protein assignment (please be patient)\n";
		my (%pep1hash,%pep2hash,%pep3hash);
		my ($scFDR1,$scFDR2,$scFDR3)=splitpeptidehash(\%peptidehash,\%proteinhash,\%pep1hash,\%pep2hash,\%pep3hash);
		#print "pep1hash: $scFDR1\n";
		#print "pep2hash: $scFDR2\n";
		#print "pep3hash: $scFDR3\n";
		

		# FDR_filtering
		my ($orig_initial_outfile_fdr,$orig_FDR)=($paramhash{initial_outfile_fdr},$paramhash{FDR});
		my (%pro1hash,%pro2hash,%pro3hash,%fpr1hash,%fpr2hash,%fpr3hash);
		print "\n1) most confident category: FDR filtering for scans from proteins with 3 peptides (filter the outfiles to protein FDR < 1%):\n";
		print LOGFILE "\n1) most confident category: FDR filtering for scans from proteins with 3 peptides (filter the outfiles to protein FDR < 1%):\n";
		#$paramhash{initial_outfile_fdr}=$paramhash{FDR}=1;
		XcorrDcnFiltering(\%pep3hash, \%pro3hash, \%fpr3hash, \%runhash, $group, 0, \%featurehash, 1);
		print "\n2) median confident category: FDR filtering for scans from proteins with 2 peptides (filter the outfiles to protein FDR < 1%):\n";
		print LOGFILE "\n2) median confident category: FDR filtering for scans from proteins with 2 peptides (filter the outfiles to protein FDR < 1%):\n";
		#$paramhash{initial_outfile_fdr}=$paramhash{FDR}=1;
		XcorrDcnFiltering(\%pep2hash, \%pro2hash, \%fpr2hash, \%runhash, $group, 0, \%featurehash, 1);
		print "\n3) least confident category: FDR filtering for scans from proteins with 1 peptide:\n";
		print LOGFILE "\n3) least confident category: FDR filtering for scans from proteins with 1 peptide:\n";
		($paramhash{initial_outfile_fdr},$paramhash{FDR})=($orig_initial_outfile_fdr,$orig_FDR);
		my ($scanFDRcut,$proteinFDR)=XcorrDcnFiltering(\%pep1hash, \%pro1hash, \%fpr1hash, \%runhash, $group, 0, \%featurehash, 1);
		#$paramhash{initial_outfile_fdr}=$orig_initial_outfile_fdr; # recover the original initial_outfile_fdr

		# merge 3 hashes
		print "\nMerging results from 3 outfile categories:\n";
		print LOGFILE "\nMerging results from 3 outfile categories:\n";

		print "(Originally starting with $ms2totalN MS2 scans and $outfileTotalN outfiles)\n";
		print LOGFILE "(Originally starting with $ms2totalN MS2 scans and $outfileTotalN outfiles)\n";


		mergePeptidehashes(\%peptidehash,\%pep1hash,\%pep2hash,\%pep3hash);
		undef(%proteinhash);
		$idutils->create_proteinhash(\%proteinhash, \%peptidehash, \%dbhash);
		count_nomod(\%proteinhash);

		my $proteinFDR=printFilteringResults(\%peptidehash,\%proteinhash,\%fprhash,$org_good,$org_bad);

		UpdateProteinQvalue_pephashVersion(\%proteinhash,\%pep3hash,1);
		UpdateProteinQvalue_pephashVersion(\%proteinhash,\%pep2hash,1);
		UpdateProteinQvalue_pephashVersion(\%proteinhash,\%pep1hash, $proteinFDR);
=head # deactivate on 7/19/16
		# 1% protein FDR estimation
		print "\nCalculating protein q-value  (please be patient)\n";
		print LOGFILE "\nCalculating protein q-value  (please be patient)\n";
		if ( $proteinFDR > 1+0.9  ) # need more filtering to achieve ~1% pro. FDR
		{
			my ($T3,$D3)=protein_TD_num(\%pep3hash);
			my ($T2,$D2)=protein_TD_num(\%pep2hash);
			$orig_initial_outfile_fdr=$paramhash{initial_outfile_fdr};
			$paramhash{initial_outfile_fdr}=$scanFDRcut;
			#($scanFDRcut,$proteinFDR)=XcorrDcnFiltering(\%pep1hash, \%pro1hash, \%fpr1hash, \%runhash, $group, 0, \%featurehash, 0,1,$T2+$T3,$D2+$D3); # no printing
			($scanFDRcut,$proteinFDR)=XcorrDcnFiltering(\%pep1hash, \%pro1hash, \%fpr1hash, \%runhash, $group, 0, \%featurehash, 1,1,$T2+$T3,$D2+$D3); # print for debug
			UpdateProteinQvalue_pephashVersion(\%proteinhash,\%pep1hash,1);
			$paramhash{initial_outfile_fdr}=$orig_initial_outfile_fdr;
		}
		else # ~1% pro. FDR already achieved
		{
			#print "~1% pro. FDR already achieved\n";
			UpdateProteinQvalue_pephashVersion(\%proteinhash,\%pep1hash,1);
		}

		# 3 more levels
		pseudo_protein_Qvalue(\%proteinhash,\%proteinhash1,\%peptidehash1,\%featurehash1,\%paramhash,$scanFDRcut/10,0.1);
		pseudo_protein_Qvalue(\%proteinhash,\%proteinhash1,\%peptidehash1,\%featurehash1,\%paramhash,$scanFDRcut/100,0.01);
		pseudo_protein_Qvalue(\%proteinhash,\%proteinhash1,\%peptidehash1,\%featurehash1,\%paramhash,0,0);
=cut
#=head
		# build 1FDR hashes to print 1% FDR results
		my (%peptidehash_1pct,%proteinhash_1pct);
		build_1pct_hash(\%peptidehash,\%proteinhash,\%dbhash,\%peptidehash_1pct,\%proteinhash_1pct);
=head
		foreach my $pep (keys %peptidehash)
		{
			my $proteinQ=100;
			foreach my $pro (keys %{$peptidehash{$pep}{proteins}})
			{
	#if (!defined($proteinhash{$pro})) { die "not defined $pro in proteinhash\n"; }
	if (!defined($proteinhash{$pro})) { next; }
				if (defined($proteinhash{$pro}{Qvalue}) and  
				$proteinQ>$proteinhash{$pro}{Qvalue} )
				{ $proteinQ=$proteinhash{$pro}{Qvalue}; }
			}
			$peptidehash{$pep}{proteinQ}=$proteinQ;
		}
		foreach my $pep (keys %peptidehash)
		{
			if ( $peptidehash{$pep}{proteinQ}<=1 )
			{ $peptidehash_1pct{$pep}=clone($peptidehash{$pep}); }
		}

		$idutils->create_proteinhash(\%proteinhash_1pct, \%peptidehash_1pct, \%dbhash);
=cut
		# deactivate on 4/27/15
		#print "For 1% unique protein FDR:\n";
		#print LOGFILE "For 1% unique protein FDR:\n";
		#printFilteringResults(\%peptidehash_1pct,\%proteinhash_1pct,\%fprhash,$org_good,$org_bad);
#=cut
		} # end of 'if ( $paramhash{unique_protein_or_peptide} eq 'protein' )'
		elsif ( $paramhash{unique_protein_or_peptide} eq 'peptide' )
		{
			#print "\nSplitting scans into 3 categories according to multiple precursor ions and modification forms for each peptide (please be patient)\n";
			#print LOGFILE "\nSplitting scans into 3 categories according to multiple precursor ions and modification forms for each peptide (please be patient)\n";
			print "\nThen three layers of scan categorization and into three groups:\n";
			print LOGFILE "\nThen three layers of scan categorization and into three groups:\n";
			my (%pep1hash,%pep2hash,%pep0hash,%pep0_1hash,%pep0_2hash);
			splitpeptidehash_multiPPI(\%peptidehash,\%pep1hash,\%pep2hash,\%pep0hash,\%pep0_1hash,\%pep0_2hash);
			# FDR_filtering
			#my (%pro1hash,%pro2hash,%fpr1hash,%fpr2hash,%pro0hash,%fpr0hash);
			my (%pro1hash,%pro2hash,%fpr1hash,%fpr2hash,%pro0_1hash,%fpr0_1hash,%pro0_2hash,%fpr0_2hash);
			print "\nLayer 1 (group 1): peptides matched by precursor ions of at least 2 different charge states (e.g. doubly charged 601 m/z and triply charged 401 m/z, only charge state used in the program)\n";
			print LOGFILE "\nLayer 1 (group 1): peptides matched by precursor ions of at least 2 different charge states (e.g. doubly charged 601 m/z and triply charged 401 m/z, only charge state used in the program)\n";
			XcorrDcnFiltering(\%pep2hash, \%pro2hash, \%fpr2hash, \%runhash, $group, 0, \%featurehash, 1);
			print "\nIn addition, include other PSMs that carry the exact same sequences with the above peptides, but have differentially modified residues (e.g. multiple modification forms)\n";
			print LOGFILE "\nIn addition, include other PSMs that carry the exact same sequences with the above peptides, but have differentially modified residues (e.g. multiple modification forms)\n";
			print "\nLayer 2 (group 2): peptides matched by precursor ions of only one charge state, and carrying multiple modification forms:\n";
			print LOGFILE "\nLayer 2 (group 2): peptides matched by precursor ions of only one charge state, and carrying multiple modification forms:\n";
			XcorrDcnFiltering(\%pep1hash, \%pro1hash, \%fpr1hash, \%runhash, $group, 0, \%featurehash, 1);
			print "\nLayer 3 (group 3): peptides matched by precursor ions of only one charge state, and carrying single modification form, but multiple scans.\n";
			print LOGFILE "\nLayer 3 (group 3): peptides matched by precursor ions of only one charge state, and carrying single modification form, but multiple scans.\n";
			#XcorrDcnFiltering(\%pep0hash, \%pro0hash, \%fpr0hash, \%runhash, $group, 0, \%featurehash, 1);
			XcorrDcnFiltering(\%pep0_2hash, \%pro0_2hash, \%fpr0_2hash, \%runhash, $group, 0, \%featurehash, 1);
			print "\nLayer 3 (group 4): peptides matched by precursor ions of only one charge state, and carrying single modification form, and single scan.\n";
			print LOGFILE "\nLayer 3 (group 4): peptides matched by precursor ions of only one charge state, and carrying single modification form, and single scan.\n";
			XcorrDcnFiltering(\%pep0_1hash, \%pro0_1hash, \%fpr0_1hash, \%runhash, $group, 0, \%featurehash, 1);


			# Merge results from 2 categories 
			print "\nMerging results from 4 outfile categories:\n";
			print LOGFILE "\nMerging results from 4 outfile categories:\n";

			print "(Originally starting with $ms2totalN MS2 scans and $outfileTotalN outfiles)\n";
			print LOGFILE "(Originally starting with $ms2totalN MS2 scans and $outfileTotalN outfiles)\n";

			mergePeptidehashes(\%peptidehash,\%pep1hash,\%pep2hash,\%pep0_1hash,\%pep0_2hash);
			undef(%proteinhash);
			$idutils->create_proteinhash(\%proteinhash, \%peptidehash, \%dbhash);
			count_nomod(\%proteinhash);
			my $proteinFDR=printFilteringResults(\%peptidehash,\%proteinhash,\%fprhash,$org_good,$org_bad);
			UpdateProteinQvalue(\%proteinhash, \%proteinhash,$proteinFDR);

		# build 1FDR hashes to print 1% FDR results
		my (%peptidehash_1pct,%proteinhash_1pct);
		build_1pct_hash(\%peptidehash,\%proteinhash,\%dbhash,\%peptidehash_1pct,\%proteinhash_1pct);
=head
		foreach my $pep (keys %peptidehash)
		{
			if ( $peptidehash{$pep}{Qvalue}<=1 ) 
			{ $peptidehash_1pct{$pep}=clone($peptidehash{$pep}); }
		}
		$idutils->create_proteinhash(\%proteinhash_1pct, \%peptidehash_1pct, \%dbhash);
=cut
		#print "\nFor 1% unique peptide FDR:\n";
		#print LOGFILE "\nFor 1% unique peptide FDR:\n";
		#printFilteringResults(\%peptidehash_1pct,\%proteinhash_1pct,\%fprhash,$org_good,$org_bad);
		}
		} # whether perform mupltistep_FDR_filterin
		else { print "\nBypass multistep FDR filtering\n"; }

#=cut # for debug		
	foreach my $pro (keys %proteinhash) { my $seq=$proteinhash{$pro}{sequence}; if (length($seq)==0) { die "protein $pro sequence: $seq\n\n"; } }

		if ($paramhash{'combine_filtermethods'}){
			print "\n";
			print "Rerunning data without mass accuracy filtering\n";
			print LOGFILE "\nRerunning data without mass accuracy filtering\n";
			my (%peptidehash2, %proteinhash2, %fprhash2, %runhash2);
	  	 	my $hash = retrieve("$group_dir/preMA_runhash");
			%runhash2 = %$hash;
			$hash = retrieve("$group_dir/preMA_peptidehash");
			%peptidehash2 = %$hash;
			my $printResults=1;
			XcorrDcnFiltering(\%peptidehash2, \%proteinhash2, \%fprhash2, \%runhash2, $group, 1, \%featurehash, $printResults);
			print "Combining data from the different filtering methods ... ";
			print LOGFILE "Combining data from the different filtering methods ... ";
			while (my ($peptide, $pephash) = each %peptidehash2){
				if (!defined($peptidehash{$peptide})){
					%{$peptidehash{$peptide}} = %$pephash;
				} else {
					while (my ($outfile, $outhash) = each %{$$pephash{'outfiles'}}){
						if (!defined($peptidehash{$peptide}{'outfiles'}{$outfile})){
							%{$peptidehash{$peptide}{'outfiles'}{$outfile}} = %$outhash;
							$peptidehash{$peptide}{'proteins'}{$$outhash{'protein'}}++;
						}
					}
				}
			}
			my %tempproteinhash;
			$idutils->create_proteinhash(\%tempproteinhash, \%peptidehash, \%dbhash);
			%proteinhash = %tempproteinhash;
			count_nomod(\%proteinhash);
			printf "Finally Accepted %d proteins and %d peptides\n", scalar(keys %proteinhash), scalar(keys %peptidehash);
			printf LOGFILE "Finally Accepted %d proteins and %d peptides\n", scalar(keys %proteinhash), scalar(keys %peptidehash);
		}
#=cut # for debug		
	}#end of bypassfiltering's else

#=head #for debug
	#emPAI calculation
	if ($paramhash{'abundance_index'} == 1){
		$idutils->calc_emPAI(\%paramhash, \%proteinhash);
	}
	#Group proteins
	my %onehitdels;
	my %onehitpeptides;

	#  one hit wonders removal
	if ($paramhash{'one_hit_wonders_removal'} > 0){
		my %tempproteinhash = %proteinhash;
		for my $protein (keys %tempproteinhash){
			next if (!defined($tempproteinhash{$protein}{'total'}));
			#next if ($tempproteinhash{$protein}{'total'} > 1); # old criteria: peptide level
			next if ($tempproteinhash{$protein}{'occurrence'} > 2); # new criteria: PSM level 
			for my $peptide (keys %{$tempproteinhash{$protein}{'peptides'}}){
		
				my $peptide_num = scalar(keys %{$peptidehash{$peptide}{'proteins'}});
				#Check for multiple charge state
				next if (multiple_charge($peptidehash{$peptide}{'outfiles'}));
				#Get best outfile to do XCorr and dCn cutoff
				my ($XCorr, $dCn, $charge) = best_outfile($peptidehash{$peptide}{'outfiles'});
				
				next if (!defined($dCn));
				if ($dCn <= $paramhash{'one_hit_wonders_min_dCn'}){
					delete  $tempproteinhash{$protein};
				          delete $proteinhash{$protein};
        				  $onehitdels{$protein} = $peptidehash{$peptide}{'orig_peptide'};
					$onehitpeptides{$peptide} = $peptidehash{$peptide}{'orig_peptide'};
					if($protein =~ /Decoy/)
					{
						$protein = "##" . $protein;
					}
					
         				 delete $peptidehash{$peptide}{'proteins'}{$protein};
					next;
				}
				if ($charge == 1){
					if ($XCorr <= $paramhash{'one_hit_wonders_min_XCorr_z1'}){
						delete  $tempproteinhash{$protein};
	         				 delete $proteinhash{$protein};
  	     					   $onehitdels{$protein} = $peptidehash{$peptide}{'orig_peptide'};
						   $onehitpeptides{$peptide} = $peptidehash{$peptide}{'orig_peptide'};
                                       		   if($protein =~ /Decoy/)
                                         	   {
                                               		 $protein = "##" . $protein;
                                        	   }	
    	   					   delete $peptidehash{$peptide}{'proteins'}{$protein};
						next;
					} 
				} elsif ($charge == 2){
					if ($XCorr <= $paramhash{'one_hit_wonders_min_XCorr_z2'}){
						delete  $tempproteinhash{$protein};
	     					     delete $proteinhash{$protein};
  	     					   $onehitdels{$protein} = $peptidehash{$peptide}{'orig_peptide'};
						$onehitpeptides{$peptide} = $peptidehash{$peptide}{'orig_peptide'};
                                      		 if($protein =~ /Decoy/)
                                        	{
                                                	$protein = "##" . $protein;
                                        	}
						delete $peptidehash{$peptide}{'proteins'}{$protein};
						next;
					}

				} 
				else 
				{
					if ($XCorr <= $paramhash{'one_hit_wonders_min_XCorr_z3'}){
						delete  $tempproteinhash{$protein};
	      					    delete $proteinhash{$protein};
  	   					     $onehitdels{$protein} = $peptidehash{$peptide}{'orig_peptide'};
							$onehitpeptides{$peptide} = $peptidehash{$peptide}{'orig_peptide'};
                                       			if($protein =~ /Decoy/)
                                        		{
                                                		$protein = "##" . $protein;
                                        		}
							delete $peptidehash{$peptide}{'proteins'}{$protein};
							next;
					}
 
				}

				if (($tempproteinhash{$protein}{'peptides'}{$peptide}{'tryptic'} == 0) || ($paramhash{'one_hit_wonders_removal'} == 2 && $tempproteinhash{$protein}{'peptides'}{$peptide}{'tryptic'} == 1) || ($peptidehash{$peptide}{'mis'} > $paramhash{'one_hit_wonders_mis'}) || ($peptidehash{$peptide}{'mod'} > $paramhash{'one_hit_wonders_mods'})){
					delete  $tempproteinhash{$protein};
					delete $proteinhash{$protein};
					$onehitdels{$protein} = $peptidehash{$peptide}{'orig_peptide'};
					$onehitpeptides{$peptide} = $peptidehash{$peptide}{'orig_peptide'};
                                        if($protein =~ /Decoy/)
                                        {
                                             $protein = "##" . $protein;
                                        }
					delete $peptidehash{$peptide}{'proteins'}{$protein};
					next;
				}
			}
		}
		printf "\nRemoving %d one-hit wonders\n", scalar(keys %onehitdels);
		printf LOGFILE "\nRemoving %d one-hit wonders\n", scalar(keys %onehitdels);
	}
	 elsif ($paramhash{'one_hit_wonders_removal'} == -1){
		my %tempproteinhash = %proteinhash;
		for my $protein (keys %tempproteinhash){
			if (!defined($tempproteinhash{$protein}{'total'})){
				print Dumper($tempproteinhash{$protein});print "undefined protein in tempproteinhash $protein\n"; exit;
			}
			#next if ($tempproteinhash{$protein}{'total'} > 1); # old criteria: peptide level
			next if ($tempproteinhash{$protein}{'occurrence'} > 2); # new criteria: PSM level
			for my $peptide (keys %{$tempproteinhash{$protein}{'peptides'}}){
				
				delete $tempproteinhash{$protein};
				delete $proteinhash{$protein};
				$onehitdels{$protein} = $peptidehash{$peptide}{'orig_peptide'};
				delete $peptidehash{$peptide}{'proteins'}{$protein};
			}
		}
		printf "Deleted %d peptides (one-hit wonder)\n", scalar(keys %onehitdels);
		printf LOGFILE "Deleted %d peptides (one-hit wonder)\n", scalar(keys %onehitdels);
	}

	# min_protein_SC # added on 5/26/15 as required for Rapid LC-MS/MS analysis
	if ( defined($paramhash{'min_protein_SC'}) and $paramhash{'min_protein_SC'}>1 )
        {
                my %smallSCproteins; my $k=0;
                my %tempproteinhash = %proteinhash;
                for my $protein (keys %tempproteinhash){
                        if (!defined($tempproteinhash{$protein}{'occurrence'})){
                                print Dumper($tempproteinhash{$protein});
                                print "undefined protein in tempproteinhash $protein\n";
                                exit;
                        }
                        next if ($tempproteinhash{$protein}{'occurrence'} >= $paramhash{'min_protein_SC'});
                        for my $peptide (keys %{$tempproteinhash{$protein}{'peptides'}}){
                                delete $tempproteinhash{$protein};
                                delete $proteinhash{$protein};
                                $smallSCproteins{$protein}='';
                                delete $peptidehash{$peptide}{'proteins'}{$protein};
                        }
                }
                printf "Deleted %d proteins due to minimum spctral counts requirement (%d)\n", scalar(keys %smallSCproteins),$paramhash{'min_protein_SC'};
                printf LOGFILE "Deleted %d proteins due to minimum spctral counts requirement (%d)\n", scalar(keys %smallSCproteins),$paramhash{'min_protein_SC'};
	}

############    dded by xusheng on 8/8/2012 ############################
               if($paramhash{'keep_decoy'} eq '1')
               {
                   my $txtfile="IDwDecoy.txt";
                   save_Filtering($group_dir,\%proteinhash,\%peptidehash,$txtfile);
               }

    #protein and peptide fdr calculation (after removing one hit wonder)
    my ($bad,$good)=(0,0);

      my ($proteinfpr);
     ($proteinfpr,$good,$bad)=pro_fpr(\%proteinhash);
=head
    if ( $paramhash{one_hit_wonders_removal}==0 )
    {
      print "Total protein FDR: $proteinfpr% ($bad decoy(s); $good target(s))\n";
      print LOGFILE "Total protein FDR: $proteinfpr%  ($bad decoy(s); $good target(s))\n";
    }
    else
    {
      printf "Accepting %d total proteins: %d targets and %d decoys (FDR = %.2f%%)\n",$bad+$good,$good,$bad,$proteinfpr;
      printf LOGFILE "Accepting %d total proteins: %d targets and %d decoys (FDR = %.2f%%)\n",$bad+$good,$good,$bad,$proteinfpr;
      #print "Accepting %d total proteins after deleting one hit wonder: $proteinfpr% ($bad decoy(s); $good target(s))\n";
      #print LOGFILE "Total protein FDR after deleting one hit wonder: $proteinfpr%  ($bad decoy(s); $good target(s))\n";
    }
=cut
	$fprhash{'protein_bad'} = $bad;

    ($bad,$good)=(0,0);
  # if one-hit-wonder are removed, the final protein and peptide FDR FDRs will be re-calculated
  #if ( $paramhash{one_hit_wonders_removal} !=0 ){
  if ( $paramhash{one_hit_wonders_removal} !=0 or defined($paramhash{'min_protein_SC'}) and $paramhash{'min_protein_SC'}>1){
    my %temp_peptidehash;
    foreach my $peptide (keys %peptidehash)
    {
################ remove the one hit wonder peptide except for two proteins shared the same peptide
        if ((scalar keys %{$peptidehash{$peptide}{'proteins'}}) eq '0' )
        {
               delete $peptidehash{$peptide};
               next;
        }
	foreach my $protein (keys %{$peptidehash{$peptide}{'proteins'}})
	{
		$temp_peptidehash{$peptide} = $protein;
	}
    }
    foreach my $peptide (keys %temp_peptidehash)
    {	
	if($temp_peptidehash{$peptide} =~ /Decoy/)
	{
		$bad++;
	}
	else
	{
		$good++;
	}
    }

        my $peptidefpr=sprintf("%.2f",$bad/($good)*100);
	my $total_protein = scalar(keys %proteinhash);
	my $total_peptide = scalar keys (%peptidehash);
    printf "  Unique peptides: $good targets and $bad decoys (FDR = %.2f%%)\n", $peptidefpr;
    printf LOGFILE "  Unique peptides: $good targets and $bad decoys (FDR = %.2f%%)\n", $peptidefpr;
    #print "Total peptide FDR after deleting one hit wonder: $peptidefpr% ($bad decoy(s); $good target(s))\n";
    #print LOGFILE "Total peptide PDR after deleting one hit wonder: $peptidefpr% ($bad decoy(s); $good target(s))\n";
 # }

#	($peptidefpr,$good,$bad)=pep_fpr(\%peptidehash);
	
	$fprhash{'final_peptide_fpr'}=$peptidefpr;
    $fprhash{'final_protein_fpr'}=$proteinfpr;
    	
    #    my %tempproteinhash;
     #    $idutils->create_proteinhash(\%tempproteinhash, \%peptidehash, \%dbhash);
      #   %proteinhash = %tempproteinhash;
	# uniq scan FPR
        my $uniqScanD=0; my $uniqScanT=0;
        while (my ($intpep, $intpephash) = each %peptidehash)
        {
                my %chargeHash;
                while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}})
                {
                        $chargeHash{$$outhash{charge}}='';
                }
                $$intpephash{'uniq_scan_num'}=scalar(keys %chargeHash);
                if ( $$intpephash{random}==1 ) { $uniqScanD+=$$intpephash{'uniq_scan_num'}; }
                else { $uniqScanT+=$$intpephash{'uniq_scan_num'}; }
        }

        #printf "Unique precursor ion FDR: %.2f%% (%d decoy(s); %d target(s))\n", $uniqScanD*100/$uniqScanT,$uniqScanD,$uniqScanT;
        #printf LOGFILE "Unique scan FPR:%.2f%% (%d decoy(s); %d target(s))\n", $uniqScanD*100/$uniqScanT,$uniqScanD,$uniqScanT;
 #       printf LOGFILE "Unique scan FPR = %.2f: %d (%d decoys) unique scans passed\n", $uniqScanD*100/$uniqScanT,$uniqScanT+$uniqScanD,$uniqScanD;

	#foreach my $pro (keys %proteinhash) { my $seq=$proteinhash{$pro}{sequence}; if (length($seq)==0) { die "protein $pro sequence: $seq\n\n"; } }
  } # end of 'if ( $paramhash{one_hit_wonders_removal} !=0 )'


	# protein grouping (final grouping)
	#my ($grouphash, $group_num, $subgroup_num) = $gp->group_proteins(\%proteinhash, \%peptidehash, \%paramhash);
	#my ($grouphash, $decoy_grouphash);
	my ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num)=update_grouphash(\%proteinhash, \%peptidehash);
	#foreach my $pro (keys %$grouphash) { print $pro,',', $grouphash->{$pro}->{group},';'; }
        #my %best_proteins;
        #assignBestProtein2Peptide(\%peptidehash, $grouphash, $decoy_grouphash, \%best_proteins); #print scalar(keys %best_proteins),"\n";
        #assignOrder2group($grouphash, \%best_proteins);
        #assignOrder2group($decoy_grouphash, \%best_proteins);
        #my ($group_num, $subgroup_num)=get_group_num($grouphash);
        #my ($decoy_group_num, $decoy_subgroup_num)=get_group_num($decoy_grouphash);

	#if ( $paramhash{one_hit_wonders_removal} !=0 ){
	if ( $paramhash{one_hit_wonders_removal} !=0 or defined($paramhash{'min_protein_SC'}) and $paramhash{'min_protein_SC'}>1) {
	printf "  Unique proteins: %d targets and %d decoys (FDR = %.2f%%)\n",$subgroup_num,$decoy_group_num,$decoy_group_num*100/$subgroup_num;
	printf LOGFILE "  Unique proteins: %d targets and %d decoys (FDR = %.2f%%)\n",$subgroup_num,$decoy_group_num,$decoy_group_num*100/$subgroup_num;
	}

	if (!(-e "$group_dir\/misc")) { system("mkdir $group_dir\/misc"); }
	open (OUT, ">$group_dir\/misc\/idsum.db");
	for my $protein (keys %proteinhash){
		print OUT ">$protein ";
		print OUT "$dbhash{$protein}{'annotation'}\n";
		print OUT "$dbhash{$protein}{'sequence'}\n";
	}

	#Add Pseudo_reverse
	if ($paramhash{'bypass_filtering'} != 1){
		for my $protein (keys %proteinhash){
			my $hash = $dbhash{$protein};
  		my $sequence = $$hash{'sequence'};
	  	my $peparray = $utils->get_peptides($sequence, 1, "trypsin");
  		my @revpeparray = @$peparray;
  		for (@revpeparray){
	    	my $mid = $_; $mid =~ s/([a-zA-Z\-]+)\.([A-Za-z]+)\.([a-zA-Z\-]+)/$2/;
  	  	chop($mid);
    		my $revmid = reverse $mid;
	    	$sequence =~ s/$mid/$revmid/;
  		}
	  	print OUT ">Decoy__$protein $$hash{'annotation'}\n";
  		print OUT "$sequence\n";
		}
	}
	close OUT;

	#Generate html file
	%{$fprhash{'params'}} = %paramhash;
	#print "Creating html files for $group\n";
	#print LOGFILE "Creating html files for $group\n";
########## added by xusheg ########################
#	in order to put it into html directory
#	$group_dir="";
##############################################

	# final Unique protein FPR calculation
	#foreach my $pro (keys %$grouphash) { print $pro,',', $grouphash->{$pro}->{group},';'; }
	#foreach my $pro (keys %proteinhash) { my $seq=$proteinhash{$pro}{sequence}; if (length($seq)==0) { die "protein $pro sequence: $seq\n\n"; } }
	print "\nGenerating HTML files ...\n";
	my ($unique_protein_fpr,$subgroupnum) = $html->gen_IDHtml(\%proteinhash, \%peptidehash, \%fprhash, $group_dir, $database, $grouphash); # html deactivate for debug
	$html->gen_IDwGHtml(\%proteinhash, \%peptidehash, \%fprhash, $group_dir, $database, $grouphash); # html deactivate for debug

	#printf "\nUnique protein FPR = %.2f%% ($fprhash{'protein_bad'} decoy(s); $subgroupnum unique protein(s))\n",$unique_protein_fpr;

	#print "\n";

	# print pepXML for accepted outfiles
	if ($pepXML_mode and defined($paramhash{output_pepXML}) and $paramhash{output_pepXML})
        {
                print "Generating pepXML for accepted outfiles ...\n";

		my (%acceptedOutHash);
                # filter %totalOutHash: only accepted outfiles retained
                buildAcceptedPSMpepXML(\%totalOutHash,\%peptidehash,\%acceptedOutHash,\%totalParaHash);

		# delete %totalOutHash
		undef %totalOutHash;

		# print filtered %acceptedOutHash
		$pepxmlParser->printPepXML(\%totalParaHash,\%acceptedOutHash,"$group_dir/html/accepted_PSM",10);

                # save filtered %acceptedOutHash to file
                store \%acceptedOutHash, "$group_dir/html/accepted_PSM.hash";

		# print fractionated hash
		printFractionHash(\%acceptedOutHash,\%totalParaHash,$group_dir);

                print "\n";
        }



	#printf LOGFILE "\nUnique protein FPR = %.2f%% ($fprhash{'protein_bad'} decoy(s); $subgroupnum unique protein(s))\n",$unique_protein_fpr;
#=cut  #for debug
#=head #for debug
	####ADDED DMD July 9, 2010#########Peptide position in Protein
	rm_zombie(\%proteinhash,$grouphash);
	pep_pos(\%proteinhash);
##### comment out ###################
	coverage(\%proteinhash);

	#foreach my $pro (keys %$grouphash) { print $pro,',', $grouphash->{$pro}->{group},';'; }
	save_Hashes("$group_dir\/misc\/idsum.db", $group, $group_dir, \%proteinhash, \%peptidehash,$grouphash);
	#print "\n";
	
	%{$sumhash{$group}{'protein'}} = %proteinhash;
	%{$sumhash{$group}{'peptide'}} = %peptidehash;
	$sumhash{$group}{'database'} = $database;

	#Create modified peptide and protein hash if necessary
	my $mods = $paramhash{'mods'};
	if ($mods ne 0){
		print "Creating display with only modified peptides\n";
		print LOGFILE "Creating display with only modified peptides\n";
		my $modgroup_dir = getcwd(); $modgroup_dir = "$modsub_dir" if (defined($modsub_dir));
		$modgroup_dir .= "\/sum\_$group"; $modgroup_dir .= "\_$save_dir" if ($save_dir ne "");
		$modgroup_dir .= "_mod";

		my %modpeptidehash;
		for my $peptide (keys %peptidehash){
			if ($peptide =~ /[$mods][\#\@\%\&\~\$\^\*]/){
				%{$modpeptidehash{$peptide}} = %{$peptidehash{$peptide}};
			}
		}
		if (scalar(keys %modpeptidehash) > 0){
			mkdir $modgroup_dir; 
			system("cp $paramsfile $modgroup_dir");
			if ($paramhash{'modpairs'} == 1){
				my $pairs = 0;
				for my $peptide (keys %modpeptidehash){
					my $nomod = $peptide; $nomod =~ s/([$mods])[\#\@\%\&\~\$\^\*]/$1/g;
					if (defined($peptidehash{$nomod})){
						%{$modpeptidehash{$nomod}} = %{$peptidehash{$nomod}};
						$pairs++;
					}
				}
				print "Modified peptide pairs: $pairs";
				print LOGFILE "Modified peptide pairs: $pairs";
			}
########################	print "\n";
				
				
			my %modproteinhash;
			$idutils->create_proteinhash(\%modproteinhash, \%modpeptidehash, \%dbhash);
			count_nomod(\%modproteinhash);

			my %modfprhash;
			
			my ($peptidefpr,$proteinfpr,$good,$bad);
			($peptidefpr,$good,$bad)=pep_fpr(\%modpeptidehash);
			($proteinfpr,$good,$bad)=pro_fpr(\%modproteinhash);
			
			$modfprhash{'peptide_fpr'}= $peptidefpr;
			$modfprhash{'protein_fpr'}= $proteinfpr;
###### added by yanji	
			$modfprhash{'protein_bad'}= $bad;
			#print "TTTTTTTTTT$bad\sTTTTTTTTTTTT\n";
			$modfprhash{'final_peptide_fpr'}=$peptidefpr;
			$modfprhash{'final_protein_fpr'}=$proteinfpr;
###### end of addition
	
			#Create database with only proteins in protein hash
			open (OUT, ">$modgroup_dir/idsum_mod.db");
			for my $protein (keys %proteinhash){
				print OUT ">$protein ";
					print OUT "$dbhash{$protein}{'annotation'}\n";
				print OUT "$dbhash{$protein}{'sequence'}\n";
			}
			#Add Keratin and trypsin
#			open (IN, "</home/xwang4/database/static/keratryp6.db");
#			while (<IN>){ print OUT "$_"; }
#			close IN;	close OUT;
			#my ($modgrouphash, $modgroup_num, $modsubgroup_num) = $gp->group_proteins(\%modproteinhash, \%modpeptidehash, \%paramhash);
			#my $modgrouphash = clone($grouphash);
			#my $decoy_modgrouphash = clone($decoy_grouphash);
			my ($modgrouphash, $decoy_modgrouphash, $modgroup_num, $modsubgroup_num, $decoy_modgroup_num, $decoy_modsubgroup_num)=update_grouphash(\%modproteinhash,\%modpeptidehash);

    #    		my ($grouphash, $group_num, $subgroup_num) = $gp->group_proteins(\%proteinhash, \%peptidehash, \%paramhash);


			
				#Generate html file
			print "\nCreating html files for $group modification\n";
			print LOGFILE "\nCreating html files for $group modification\n";	
#			$html->gen_IDHtml(\%modproteinhash, \%modpeptidehash, \%fprhash, $modgroup_dir, $database, $modgroup_num, $modsubgroup_num, 0, 1);
            $html->gen_IDmodHtml(\%modproteinhash, \%modpeptidehash, \%modfprhash, $modgroup_dir, $database, $modgrouphash);
            $html->gen_IDwGmodHtml(\%modproteinhash, \%modpeptidehash, \%modfprhash, $modgroup_dir, $database, $modgrouphash);


			pep_pos(\%modproteinhash);
			coverage(\%modproteinhash);
			save_Hashes("$modgroup_dir/idsum_mod.db", $group, $modgroup_dir, \%modproteinhash, \%modpeptidehash,$modgrouphash, 1);
			print "\n";
		
			%{$modsumhash{$group}{'protein'}} = %modproteinhash;
			%{$modsumhash{$group}{'peptide'}} = %modpeptidehash;
			$modsumhash{$group}{'database'} = $database;
		} else {
			print "No modified peptides found within $group\n";
			print LOGFILE "No modified peptides found within $group\n";
		}
	}
#=cut #for debug
} # end of group processing

# Summarize all groups into one 
if (scalar(keys %{$paramhash{'groups'}})>1){
	my @folders;
	for my $group (keys %{$paramhash{'groups'}}){
		push (@folders, $group);
	}
	@folders = @runarray;
	print "Creating Summary files for all groups ...\n";
	print LOGFILE "\nCreating Summary files for all groups ...\n";
	system("cp $paramsfile $nomodsum_dir");
	 
	#Recreate a new summary peptide hash
	my (%proteinhash, %peptidehash);
	$idutils->create_sumpeptidehash(\%sumhash, \%peptidehash);
	$idutils->create_proteinhash(\%proteinhash, \%peptidehash, \%dbhash);

	# Added total_nomod, unique_nomod, and shared_nomod DMD 5/24/05
	count_nomod(\%proteinhash);

    my ($peptidefpr,$good,$bad)=pep_fpr(\%peptidehash);

	my %tempfprhash;
	my %fprhash;

	$idutils->count_fpr(\%proteinhash, \%tempfprhash, 1);
	
        my ($rate1, $rate2, $rate3) = (0,0,0);
        $rate1 = ($tempfprhash{'r1'}/($tempfprhash{'t1'}-$tempfprhash{'r1'}))*100 if ($tempfprhash{'t1'} != 0);
        $rate2 = ($tempfprhash{'r2'}/($tempfprhash{'t2'}-$tempfprhash{'r2'}))*100 if ($tempfprhash{'t2'} != 0);
        $rate3 = ($tempfprhash{'r3'}/($tempfprhash{'t3'}-$tempfprhash{'r3'}))*100 if ($tempfprhash{'t3'} != 0);

                # Mark filtered outfiles
        printf "  Unique peptides: %d targets and %d decoys (FDR = %.2f%%)\n",$good,$bad,$peptidefpr;
        printf "  Protein (1pep hit) fpr = %.2f, Protein (2pep hit) fpr = %.2f, Protein (3+pep hit) fpr = %.2f\n", $rate1, $rate2, $rate3;
        printf LOGFILE "  Unique peptides: %d targets and %d decoys (FDR = %.2f%%)\n",$good,$bad,$peptidefpr;
        printf LOGFILE "  Protein (1pep hit) fpr = %.2f, Protein (2pep hit) fpr = %.2f, Protein (3+pep hit) fpr = %.2f\n", $rate1, $rate2, $rate3;
#        printf "  Total Peptide FDR = %.2f: %d (%d decoy(s)) peptide passed \n",$peptidefpr,$good+$bad,$bad;
	
	$fprhash{'protein_bad'} = $tempfprhash{'r'};
	$fprhash{'protein_fpr'}= sprintf("%.2f",($tempfprhash{'r'}/($tempfprhash{'t'}))*100);
	$fprhash{'peptide_fpr'}= $peptidefpr;
######################################################################



	#saving PeptideHash and ProteinHash for CGI use
	open (OUT, ">$nomodsum_dir\/sum\.db");
	for my $protein (keys %proteinhash){
		print OUT ">$protein ";
		print OUT "$dbhash{$protein}{'annotation'}\n";
		print OUT "$dbhash{$protein}{'sequence'}\n";
	}
	#Add Keratin and trypsin
#	open (IN, "</home/xwang4/database/static/keratryp6.db");
#	while (<IN>){ print OUT "$_"; }
#	close IN;	close OUT;

#	my ($group_num, $subgroup_num) = $html->group_proteins(\%proteinhash, \%peptidehash, $subgroup_percent, \%paramhash);
	
	#my ($grouphash, $group_num, $subgroup_num) = $gp->group_proteins(\%proteinhash, \%peptidehash, \%paramhash);
	my ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num)=update_grouphash(\%proteinhash, \%peptidehash);

	print "\nCreating html files for all groups\n";
	print LOGFILE "\nCreating html files for all groups\n";
			
	#Generate html file
#	$html->gen_SumIDHtml(\%sumhash, \%proteinhash, \%peptidehash,\%fprhash,\@runarray, $nomodsum_dir, $database, $group_num, $subgroup_num);

        my ($unique_proten_fpr,$subgroupnum) = $html->gen_sumIDHtml(\%sumhash,\%proteinhash, \%peptidehash, \%fprhash, \@runarray, $nomodsum_dir, $database, $grouphash);
        $html->gen_sumIDwGHtml(\%sumhash,\%proteinhash, \%peptidehash, \%fprhash, \@runarray, $nomodsum_dir, $database, $grouphash);
	printf "  Unique proteins: %d targets and %d decoys (FDR = %.2f%%)\n",$subgroupnum,$fprhash{'protein_bad'}, $unique_proten_fpr;
	printf LOGFILE "  Unique proteins: %d targets and %d decoys (FDR = %.2f%%)\n",$subgroupnum,$fprhash{'protein_bad'}, $unique_proten_fpr;
	#printf "  Unique protein FPR = %.2f%% ($fprhash{'protein_bad'} decoy(s) $subgroupnum unique protein(s))\n", $unique_proten_fpr;
        #printf LOGFILE "Unique Protein FPR = %.2f%% ($fprhash{'protein_bad'} decoy(s) $subgroupnum unique protein(s))\n", $unique_proten_fpr;
	pep_pos(\%proteinhash);
	print "\nPrinting ID\.txt ...\n";
	save_Hashes("$nomodsum_dir\/sum\.db", $nomodsum_dir, $nomodsum_dir, \%proteinhash, \%peptidehash,$grouphash);
	
	# Modification 
	my $mods = $paramhash{'mods'};
	if ($mods ne 0){
		print "\nCreating Summary files for all groups (modified peptide only) ...\n";
		print LOGFILE "\nCreating Summary files for all groups (modified peptide only) ...\n";
		#Recreate a new summary peptide hash
		my (%modproteinhash, %modpeptidehash);
		$idutils->create_sumpeptidehash(\%modsumhash, \%modpeptidehash);
		if (scalar(keys %modpeptidehash)>0){
			system ("mkdir $modsum_dir") if (!(-e $modsum_dir));
			system("cp $paramsfile $modsum_dir");
			$idutils->create_proteinhash(\%modproteinhash, \%modpeptidehash, \%dbhash);

			my (%modfprhash,%modtempfprhash);
#			$idutils->count_fpr(\%modproteinhash, \%modtempfprhash, 1);
			my ($peptidefpr,$proteinfpr,$good,$bad);
			($peptidefpr,$good,$bad)=pep_fpr(\%modpeptidehash);
			($proteinfpr,$good,$bad)=pro_fpr(\%modproteinhash);

#			my $proteinfpr= ($modtempfprhash{'r'}/($modtempfprhash{'t'}-$modtempfprhash{'r'}))*100;				
			$modfprhash{'peptide_fpr'}= $peptidefpr;
			$modfprhash{'protein_fpr'}= $proteinfpr;
###### add and change by yanji			
		        #$modfprhash{'protein_bad'}= $modtempfprhash{'r'};
			$modfprhash{'protein_bad'} = $bad;			
			$modfprhash{'final_peptide_fpr'}= $peptidefpr;
 	                $modfprhash{'final_protein_fpr'}= $proteinfpr;
###### end of addition			
			
			# Added total_nomod, unique_nomod, and shared_nomod DMD 5/24/05
			count_nomod(\%modproteinhash);
		
			#saving PeptideHash and ProteinHash for CGI use
			open (OUT, ">$modsum_dir\/sum\_mod\.db");
			for my $protein (keys %proteinhash){
				print OUT ">$protein ";
				print OUT "$dbhash{$protein}{'annotation'}\n";
				print OUT "$dbhash{$protein}{'sequence'}\n";
			}
			#Add Keratin and trypsin
#			open (IN, "</home/xwang4/database/static/keratryp6.db");
#			while (<IN>){ print OUT "$_"; }
#			close IN;	close OUT;
#			my ($modgroup_num, $modsubgroup_num) = $html->group_proteins(\%modproteinhash, \%modpeptidehash, $subgroup_percent, \%paramhash);

        		#my ($modgrouphash, $modgroup_num, $modsubgroup_num) = $gp->group_proteins(\%modproteinhash, \%modpeptidehash, \%paramhash);
			my ($modgrouphash, $decoy_modgrouphash, $modgroup_num, $modsubgroup_num, $decoy_modgroup_num, $decoy_modsubgroup_num)=update_grouphash(\%modproteinhash, \%modpeptidehash);

	
			#Generate html file
#			$html->gen_SumIDHtml(\%modsumhash, \%modproteinhash, \%modpeptidehash,\%fprhash, \@folders, $modsum_dir, $database, $modgroup_num, $modsubgroup_num, 0, 1);
      			  $html->gen_sumIDmodHtml(\%modsumhash,\%modproteinhash, \%modpeptidehash, \%modfprhash, \@folders, $modsum_dir, $database, $modgrouphash);
     			   $html->gen_sumIDwGmodHtml(\%modsumhash,\%modproteinhash, \%modpeptidehash, \%modfprhash, \@folders, $modsum_dir, $database, $modgrouphash);


			pep_pos(\%modproteinhash);
			print "\n";
			save_Hashes("$modsum_dir\/sum\_mod\.db", $modsum_dir, $modsum_dir, \%modproteinhash, \%modpeptidehash,$modgrouphash, 1);
		} else {
			print "No modified peptides in all groups.\n";
			print LOGFILE "No modified peptides in all groups.\n";
		}
	}	
}

###### added by yanji, copy idsum output files to /var/www/html
my $current_user = qx[whoami];
chomp($current_user);
my @out_path = split/\//, $out_dir;

## idsum output folder or directory, not full path
my $output_folder = pop @out_path;
#=head #for debug
if ($group_count > 1) {
	$output_folder = pop @out_path;
	$output_folder = pop @out_path;
	#generate path for sum_all 
	$out_dir = "";
	shift(@out_path);
	foreach (@out_path) {
		$out_dir = $out_dir."\/".$_;
	}
	$out_dir = $out_dir."\/".$output_folder;
}

#system("sudo cp -r ./$output_folder /var/www/html/$current_user");
# parse ID.html to csv
if (0) # deactivate since -q don't need it
{
	print "\nGenerating CSV files\n";
	system ("perl /usr/local/bin/html2csv.pl ./$output_folder/ID.html > ./$output_folder/ReportID.csv");
	system ("sed -i 's/?/ /g' ./$output_folder/ReportID.csv");
	system ("sed -i 's/Legend//g' ./$output_folder/ReportID.csv");
}
if (!(-e "./$output_folder/html")) { system("mkdir ./$output_folder/html"); }
system ("mv ./$output_folder/ID.html ./$output_folder/html");
system ("mv ./$output_folder/misc/Params_id.data ./$output_folder/html");
system ("mv ./$output_folder/misc/Peptides_id.data ./$output_folder/html");
system ("mv ./$output_folder/misc/Proteins_id.data ./$output_folder/html");

# parse IDwG.html to csv
if (0) # deactivate since -q don't need it
{
	system ("perl /usr/local/bin/html2csv.pl ./$output_folder/IDwG.html > ./$output_folder/ReportID_wG.csv");
	system ("sed -i 's/?/ /g' ./$output_folder/ReportID_wG.csv");
	system ("sed -i 's/Legend//g' ./$output_folder/ReportID_wG.csv");
}
system ("mv ./$output_folder/IDwG.html ./$output_folder/html");

system ("mv $log_file ./$output_folder/");

# Core report
my $corereport_ID = 0;
my $corereport_SC_comp1 = "0";

my $core_report = "CoreReport.csv";

open(PARAM, "$paramsfile");
while(<PARAM>)
{
	chomp;
	$_ =~ s/\s//g;
	if ($_ =~ /^corereport_ID\=1/)
	{
		$corereport_ID = 1;
	} elsif ($_ =~ /^corereport_SC_comp1/)
	{
		my @t1 = split/\#/, $_;
		my @t2 = split/\=/, $t1[0];
		$corereport_SC_comp1 = $t2[1];
		last;
	}
}
close PARAM;
#print "TTTTTT$corereport_SC_comp1\n";

if ($corereport_ID == 1)
{
print "\nGenerating core reports\n";
	gen_core_report($core_report, $output_folder, $corereport_SC_comp1);
}

#system("cp -r ./$output_folder /var/www/html/$current_user");

# path under /var/www/html
my $var_path = $out_dir;

if ($var_path =~ /^\/home/) {
	$var_path =~ s/^\/home/\/var\/www\/html/;
} else {
	$var_path = "\/var\/www\/html\/".$current_user.$var_path;
}

my $basename = dirname($var_path);

#print "ttttt$out_dir\t$basename\n";

unless (-e $basename) {
	system("mkdir -p $basename");	
}

#system(qq(cp -r "./$output_folder/* $var_path" >/dev/null 2>&1));
#system("ln -s ./$output_folder/* $var_path");
system(qq(ln -s $out_dir/ $basename >/dev/null 2>&1));

my $mod_output_folder = $output_folder."_mod";
#if(defined($modsum_dir))
if (-e $mod_output_folder)
{
	# parse ID.html to csv
	if (0)
	{
		system ("perl /usr/local/bin/html2csv.pl ./$mod_output_folder/IDmod.html > ./$mod_output_folder/ReportID.csv");
		system ("sed -i 's/?/ /g' ./$mod_output_folder/ReportID.csv");
		system ("sed -i 's/Legend//g' ./$mod_output_folder/ReportID.csv");
	}
	if (-e "./$output_folder/IDmod.html") { system ("mv ./$output_folder/IDmod.html ./$output_folder/html"); }

	# parse IDwG.html to csv
	if (0)
	{
		system ("perl /usr/local/bin/html2csv.pl ./$mod_output_folder/IDwGmod.html > ./$mod_output_folder/ReportID_wG.csv");
		system ("sed -i 's/?/ /g' ./$mod_output_folder/ReportID_wG.csv");
		system ("sed -i 's/Legend//g' ./$mod_output_folder/ReportID_wG.csv");
	}
	if (-e "./$output_folder/IDwGmod.html") { system ("mv ./$output_folder/IDwGmod.html ./$output_folder/html"); }

	$corereport_ID=0;
	if ($corereport_ID == 1)
	{
		$core_report = "CoreReport_mod.csv";
		gen_core_report($core_report, $mod_output_folder, $corereport_SC_comp1);
	}

	#system("cp -r ./$mod_output_folder /var/www/html/$current_user");
	my $var_path_mod = $out_dir."_mod";
	
	system(qq(ln -s $var_path_mod $basename >/dev/null 2>&1));
	
	system ("cp ./$output_folder/$log_file ./$mod_output_folder/");
}
#=cut #for debug
#system ("mv $log_file ./$output_folder/"); # activate when for debug
#system ("mv \.idsumTmp\/Grouping_for_XCorr_filtering.txt ./$output_folder/; rm -r \.idsumTmp");
#if (-e "\.idsumTmp") {system ("rm -r \.idsumTmp");}
#if (-e "\.idsumTmp") {system ("mv  \.idsumTmp ./$output_folder/");}
if (-e "\.idsumTmp") {system ("mv  \.idsumTmp ./$output_folder/ >/dev/null 2>&1");}
#system ("mv feature_infor.txt ./$output_folder/");
#if (-e 'feature_infor_LDA.txt') { system ("mv feature_infor_LDA.txt ./$output_folder/"); }
print "\nDone!\n\n\n";
print LOGFILE "\nDone!\n\n\n";
close LOGFILE;

###### end of addition


#########################################################################################
#	SUBROUTINES									#
#########################################################################################
sub gen_core_report
{
	my ($core_report, $output_folder, $corereport_SC_comp1)=@_;
	
	my $start_record = 0;
	my %column_hash;
	my @column_array;
	
	open (COREREPORT, ">$core_report");
	
	open (REPORT, "./$output_folder/ReportID.csv");
	while(<REPORT>)
	{
		chomp;
		#print "$_\n";
		#my %column_hash;
		#my @column_array;
		my $data_row = $_;
		
		if ($data_row =~ /^Group/) { $start_record = 1; } # print "555555$data_row\n";}
		my @tmp_row = split/\,/, $data_row;
		#my $tmp_row_count = 0; foreach (@tmp_row) { print "555555\t$tmp_row_count\t$_\n"; $tmp_row_count++;}
		
		my $abundance_value = pop(@tmp_row);
		
		my $mw_value = pop(@tmp_row);		
				
		#print "RRRRRR$abundance_value\t$mw_value\n";
		
		# annotation
		my @t2_row = split/\,\"/, $data_row;
		if (defined($t2_row[1])) { $t2_row[1] =~ s/\"//g; } #$t2_row[1] =~ s/\,/ /g; }
		#print "$t2_row[1]\n";
		
		if ($start_record == 1)
		{
			#$mw_value =~ s/\"//g;
			#$mw_value =~ s/\s+//g;
			#$abundance_value = sprintf("%.1f", $abundance_value);
			#$mw_value = sprintf("%.0f", $mw_value);
			#my $tmp_row_count = 0; foreach (@tmp_row) { print "555555\t$tmp_row_count\t$_\n"; $tmp_row_count++;}
			
			my $header = "Group,Accession #,Protein Annotation,SC,TP,MW (KD),Abundance";			
				
			if ($group_count > 1)
			{
				$header = "Group,Accession #,Protein Annotation";
				#my $array_end = $#tmp_row - 2;
				#print "777777\t$#tmp_row - 3\n";
				for (my $i = 3; $i <= $#tmp_row - 1; $i++)
				{
					$tmp_row[$i] =~ s/\"//g;
					#print "333333$i\t$tmp_row[$i]\n";
					if ($tmp_row[$i] !~ /Total Peptide/)
					{
						#print "111111$i\t$tmp_row[$i]\n";
						$header = $header.",".$tmp_row[$i];
						my $tmp_name = $tmp_row[$i];
						if ($tmp_name =~ /SC /)
						{
							$tmp_name =~ s/SC //;
							$column_hash{$tmp_name} = $i;
							#print "222222$tmp_name, $i, $column_hash{$tmp_name}\n";
						} 
						
						push (@column_array, $i);
					}					
				}
				
				$header = $header.",MW (KD),Abundance,p value";
				#print "TTTTTT$header\n";
			}
			
			print COREREPORT "$header\n";
			$start_record++;
		} elsif ($start_record > 1)
		{
			$mw_value =~ s/\"//g;
			$mw_value =~ s/\s+//g;
			$abundance_value = sprintf("%.1f", $abundance_value);
			$mw_value = sprintf("%.0f", $mw_value);
			
			if ($group_count > 1)
			{
				print COREREPORT "$tmp_row[0],$tmp_row[1],\"$t2_row[1]\",";
				foreach (@column_array)
				{
					#print "EEEEEE$_\n";
					print COREREPORT "$tmp_row[$_],";
				}
				#my $abundance_value = pop(@tmp_row);
				#my $mw_value = pop(@tmp_row);
				print COREREPORT "$mw_value,$abundance_value,";
				
				# p value calculation
				if ($corereport_SC_comp1 eq "0")
				{
					print COREREPORT "\n";
				} else
				{
					$corereport_SC_comp1 =~ s/\:/\,/;
					#print "TTT$corereport_SC_comp1\n";
					my @tmp_comp = split/\,/, $corereport_SC_comp1;
					my @comp_array;
					
					my $null_value = 0;
					foreach (@tmp_comp)
					{
						#print "RRR$_\n";
						my $tmp_element = $column_hash{$_};
						#print "RRRRRRR$tmp_element\n";
						push (@comp_array, $tmp_row[$tmp_element]);
						if ($tmp_row[$tmp_element] == 0) { $null_value = 1; }
					}
					
					if ($null_value == 0)
					{
						my $pvalue = chi($comp_array[0], $comp_array[1], $comp_array[2], $comp_array[3]);
						$pvalue = sprintf "%.2E", $pvalue;
						#$pvalue =~ s/e/E/g;
						#print "RRRRRR$pvalue\n";
						print COREREPORT "$pvalue\n";
					} else { print COREREPORT "\n"; }
				}				
			} else 
			{
				print COREREPORT "$tmp_row[0],$tmp_row[1],\"$t2_row[1]\",$tmp_row[2],$tmp_row[3],$mw_value,$abundance_value\n"; #$tmp_row[7],$tmp_row[8]\n";
			}
		}	
	}
	close REPORT;
	
	close COREREPORT;
	
	#print "MMMMMM$output_folder\n";
	system ("mv $core_report ./$output_folder/");
}

sub pep_fpr
{
	my ($peptidehashref)=@_;
	my %peptidehash = %$peptidehashref;
	
	my ($bad,$good)=(0,0);
	my %temp_peptidehash;
	foreach my $peptide (keys %peptidehash)
    	{
################ remove the one hit wonder peptide except for two proteins shared the same peptide
        	if ((scalar keys %{$peptidehash{$peptide}{'proteins'}}) eq '0' )
        	{
              	 	delete $peptidehash{$peptide};
               		next;
        	}
        	foreach my $protein (keys %{$peptidehash{$peptide}{'proteins'}})
        	{
                	$temp_peptidehash{$peptide} = $protein;
        	}
    	}
    	foreach my $peptide (keys %temp_peptidehash)
    	{
       		 if($temp_peptidehash{$peptide} =~ /Decoy/)
        	{
               		 $bad++;
        	}
       		 else
        	{
               		 $good++;
        	}
    	}

=head
    foreach my $peptide (keys %peptidehash)
    {
################ remove the one hit wonder peptide except for two proteins shared the same peptide		
		if ((scalar keys %{$peptidehash{$peptide}{'proteins'}}) eq '0' )
		{
			delete $peptidehash{$peptide};
			next;
		}
############################
		foreach my $protein_name (keys %{$peptidehash{$peptide}{'proteins'}})
		{
			if($protein_name =~ /Decoy/)
			{
				$bad++;
			}
			else
			{
				$good++;
			}
		}
	}
	
    my $total_peptide = scalar keys (%peptidehash);
=cut
	my $peptidefpr=sprintf("%.2f",$bad/($good)*100);	
	return 	($peptidefpr,$good,$bad);
}

sub pro_fpr
{
	my ($proteinhashref)=@_;
	my %proteinhash = %$proteinhashref;
	
	my ($bad,$good)=(0,0);
    foreach my $protein (keys %proteinhash)
    {
			if($protein =~ /Decoy/)
			{
				$bad++;
			}
			else
			{
				$good++;
			}
	}
	
    my $total_protein = scalar keys (%proteinhash);
	my $proteinfpr=sprintf("%.2f",$bad/($good)*100);	
	return 	($proteinfpr,$good,$bad);
}

sub rm_zombie{
	my ($hash,$grouphash) = @_;
	while (my ($pro, $prohash) = each %$hash){
		my $proseq = $$prohash{'sequence'};
		if ($proseq eq '') {
			delete $$hash{$pro};
			delete $grouphash->{$pro};
		}
	}
}

sub pep_pos{

	open (LOGFILE, ">>$log_file");	

	my ($hash) = @_;
	print "\nAdding peptide position to protein hash\n";
	print LOGFILE "\nAdding peptide position to protein hash\n";
	my $total = scalar(keys %$hash);
	my $num = 0;
	while (my ($pro, $prohash) = each %$hash){
		$num++;
		#print "\rWorking on $num of $total:  $pro                     ";
		my $proseq = $$prohash{'sequence'};
		if ($proseq eq '') {die "Zero length protein:$pro\n\n";}
		#if ($proseq eq '') {delete $$hash{$pro}; next;} # deal with zombie protein term [could be dangerous]
		while (my ($pepseq, $pephash) = each %{$$prohash{'peptides'}}){
#			print $pepseq,"\t",$pephash,"\n";
			$pepseq =~ s/[\*\#\@\%\&\~\$\^]//g;
			my $fpos = index($proseq, $pepseq) + 1;
			my $lpos = $fpos + length($pepseq) - 1;
			my $position = "AA$fpos"."to"."AA$lpos";
			$$pephash{'peppos'} = $position;
		}
	}
}

sub coverage{
        open (LOGFILE, ">>$log_file");

	my ($hash) = @_;
	print "\nCalculating protein coverage \n";
	print LOGFILE "\nCalculating protein coverage \n";
	while (my ($pro, $prohash) = each %$hash){
		my $Seq = $$prohash{'sequence'};
		for my $peptide (sort keys %{$$prohash{'peptides'}}){
			my $seq = $peptide;
			$seq =~ s/[\*\#\@\%\&\~\$\^]//g;
			$Seq =~ s/($seq)/\L$1/ig;
		}
		##print "$protein\n$Seq\n";
		my $num = $Seq =~ s/([a-z])/$1/g;
  	##printf "$num vs %d\n", length($Seq);exit;
  	if ($Seq eq '') {die "Zero length protein:$pro\n$Seq\n";}
  	$$prohash{'coverage'} = $num/length($Seq)*100;
	}
}

sub multiple_charge{
	my ($hash) = @_;
	my %charges;
	for my $outfile (keys %$hash){
		my $charge = $$hash{$outfile}{'charge'};
		$charges{$charge}++;
	}
	if (scalar(keys %charges) > 1){ return 1; } else { return 0; }
}

sub best_outfile{
	my ($hash) = @_;
	my $file = "";
	for my $outfile (sort {$$hash{$b}{'XCorr'}<=>$$hash{$a}{'XCorr'}} keys %$hash){
		$file = $outfile;
		last;
	}
	return ($$hash{$file}{'XCorr'}, $$hash{$file}{'dCn'}, $$hash{$file}{'charge'});
}


sub modSiteCount
{
my ($IDmod,$modSites)=@_;

my (%pephash,%prohash,$line);

open(IN,"$IDmod");
$line=<IN>;
$line=<IN>;

# build %pephash
while(<IN>)
{
        chomp;
        my @t=split(/\;/,$_);
        #my ($pep,$pro,$pos)=($t[0],$t[1],$t[$#t]);
        my ($pep,$pro,$pos)=($t[0],$t[1],$t[$#t-1]);
	next if ($pro =~ m/Decoy/);

        $pep =~ s/\.//g; chop($pep); $pep=reverse($pep); chop($pep); $pep=reverse($pep);

        $pephash{$pep}{'protein'}{$pro}{$pos}='';
}
close IN;

# calculate mod site position for each peptide
foreach my $pep (keys %pephash)
{
        my @mod=split(//,$modSites);
        my %mods; foreach my $s (@mod) { $mods{$s}=''; }

        #print "$pep\t";
        my $p=$pep;
        while ($p =~ /[\@\%\&\~\$\^\#\*\^\~]/)
        {
                my $s=substr($p,$-[0]-1,1);
                if (defined($mods{$s})) { my $pos=$-[0]-1; $pephash{$pep}{'modsite'}{$pos}=''; }# print "$pos,"; }
                $p =~ s/[\@\%\&\~\$\^\#\*\^\~]//;
        }
        #print "\n";
}

# build %prohash
foreach my $pep (keys %pephash)
{
        my $firstPro=1;
        foreach my $pro (keys %{$pephash{$pep}{protein}})
        {
                my $firstLocation=1;
                foreach my $pos (keys %{$pephash{$pep}{protein}{$pro}})
                {
                        $pos =~ /^AA(\d+)toAA(\d+)$/;
                        my $start=$1;
			if ( $start<1 ) { next; }
                        foreach my $ms (keys %{$pephash{$pep}{modsite}})
                        {
                                my $t=$start+$ms;
                                if (!defined($prohash{$pro}{$t}))
                                {
                                        $prohash{$pro}{$t}=$firstPro;
                                        unless ($firstLocation)
                                        {
                                                $prohash{$pro}{$t}=$prohash{$pro}{$t}*2;
                                                #die "multi-location: $pep,$pro\n";
                                        }
                                }
                        }

                        $firstLocation=0;
                }
                $firstPro=0;
        }
}

# calculate phospho sites
my $count=0;
my $count1=0;
foreach my $pro (keys %prohash)
{
        #print "$pro\t";
        foreach my $ms (keys %{$prohash{$pro}})
        {
                #print "$ms($prohash{$pro}{$ms}),";
                if ( $prohash{$pro}{$ms} ) {$count++;}
                if ( $prohash{$pro}{$ms}==1 ) {$count1++;}
        }
        #print "\n";
}

return (scalar(keys %pephash),$count);
=head
open (LOGFILE, ">>$log_file");
print "\n\nModified peptides and sites information:\n";
print "  There are ",scalar(keys %pephash)," unique modified peptides (",$paramhash{'mods'},")\n";
print "  There are $count unique modified sites (",$paramhash{'mods'},")\n";
print LOGFILE "\n\nModified peptides and sites information:\n";
print LOGFILE "  There are ",scalar(keys %pephash)," unique modified peptides (",$paramhash{'mods'},")\n";
print LOGFILE "  There are $count unique modified sites (",$paramhash{'mods'},")\n";
close LOGFILE;
=cut
}


sub save_Hashes{
	my ($db, $group, $dir, $proteinhash, $peptidehash, $grouphash,$mod) = @_;

  my $file = "id.data"; my $txtfile = "ID.txt";
	if (defined($mod)){
  	$txtfile = "IDmod.txt";
	}

  if (!(-e "$dir\/misc")) { system("mkdir $dir\/misc"); }	
  eval {
    store(\%paramhash, "$dir\/misc\/Params_$file");
  };
  print "Error writing to file: $@\n" if $@;
  eval {
    store($proteinhash, "$dir\/misc\/Proteins_$file");
  };
  print "Error writing to file: $@\n" if $@;
  eval {
    store($peptidehash, "$dir\/misc\/Peptides_$file");
  };
  print "Error writing to file: $@\n" if $@;

  # print anything accepted at the user specified level (usually loose, e.g., 5%FDR)
  print_IDtxt($peptidehash,$proteinhash,$grouphash,$paramhash{unique_protein_or_peptide},100,"$dir/$txtfile",$db);
  print_IDtxt($peptidehash,$proteinhash,$grouphash,$paramhash{unique_protein_or_peptide},1,"$dir/ID\_$paramhash{unique_protein_or_peptide}_1FDR.txt",$db);
  print_IDtxt($peptidehash,$proteinhash,$grouphash,$paramhash{unique_protein_or_peptide},0,"$dir/ID\_$paramhash{unique_protein_or_peptide}_0FDR.txt",$db);
=head	
  open (OUT, ">$dir/$txtfile");
	print OUT "Database=$db\n";
  
###### changed by yanji, add rention time, RT, and intensity in the header line of ID.txt
  print OUT "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;pos;precursor_peak_intensity_percentage\n";
 # print OUT "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;RT;intensity(totIonCurrent);group;subgroup;unique;tryptic;pos\n";
###### end of change

  for my $peptide (keys %$peptidehash){
    for my $outfile (keys %{$$peptidehash{$peptide}{'outfiles'}}){
      for my $protein (keys %{$$peptidehash{$peptide}{'proteins'}}){
				next if (!defined($$proteinhash{$protein}));
				#print Dumper($$proteinhash{$protein}{'peptides'}{$peptide});exit;
				print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'peptide'};";
	  		print OUT "$protein;$$peptidehash{$peptide}{'outfiles'}{$outfile}{'path'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'MH'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'expMH'};";
				if (!defined($$peptidehash{$peptide}{'outfiles'}{$outfile}{'ppm'})){
          print OUT "N/A;";
        } else {
          print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'ppm'};";
          #print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'rawppm'};"; # for ppm debug
        }
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'XCorr'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'dCn'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'ions'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'red'};";
###### added by ynaji, print rention time and intensity
               #         print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'rt'};";
                #        print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'prec_int'};";
####### end of change
#
	
			#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'mis'};";
	  		#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'mod'};";
	  		#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'tryptic'};";
				#my $nomod = $peptide; $nomod =~ s/[\@\%\&\~\$\^\*\#]//g; 
	  		#printf OUT "%d;", length($nomod);
#	  		print OUT "$$proteinhash{$protein}{'group'};";
#	  		print OUT "$$proteinhash{$protein}{'subgroup'};";
              print OUT "$grouphash->{$protein}->{'group'};";
                       print OUT "$grouphash->{$protein}->{'subgroup'};";


				print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'unique'};";
				print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'tryptic'};";
				if (!defined($$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'})){
					print "$peptide\n";
					print Dumper($$proteinhash{$protein}{'peptides'}{$peptide});exit;
				}
				#print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'}\n";
				print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'};";
				if (defined($$peptidehash{$peptide}{'outfiles'}{$outfile}{'precursor_peak_intensity_percentage'})) { print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'precursor_peak_intensity_percentage'}"; }
				else { print OUT "NA"; }
				print OUT "\n";
			}
    }
  }
=cut
	# modified Sites information
	my $modifiedSites=0;
	if (defined($mod)) 
	{ 
		my ($aPep,$aSite)=modSiteCount("$dir/IDmod.txt",$paramhash{'mods'});
		$modifiedSites=$aSite;

		open (LOGFILE, ">>$log_file");

		print "\n\nModified peptides and sites information:\n";
		print "  There are ",$aPep," unique modified peptides (",$paramhash{'mods'},")\n";
		print "  There are $aSite unique modified sites (",$paramhash{'mods'},")\n";
		print LOGFILE "\n\nModified peptides and sites information:\n";
		print LOGFILE "  There are ",$aPep," unique modified peptides (",$paramhash{'mods'},")\n";
		print LOGFILE "  There are $aSite unique modified sites (",$paramhash{'mods'},")\n";

		my @tmp=split(//,$paramhash{'mods'});
		foreach my $t (@tmp)
		{
			($aPep,$aSite)=modSiteCount("$dir/IDmod.txt",$t);
			print "    Modified sites for $t: $aSite\n";
			print LOGFILE "    Modified sites for $t: $aSite\n";
			#print "    There are $aSite unique modified sites (",$t,")\n";
			#print LOGFILE "    There are $aSite unique modified sites (",$t,")\n";
		}


		close LOGFILE;
	}

	############################################################################################
	# CORE Report Printout
	############################################################################################
	if (defined($paramhash{'core_report'})){
		if ($paramhash{'core_report'} == 1){
			if (defined($nomodsum_dir)){ # batch idsum, create report that incorporates all samples
				if (-e("$nomodsum_dir\/CoreReport.html")){
					open (OUT, ">>$nomodsum_dir\/CoreReport.html");
				} else {
					open (OUT, ">$nomodsum_dir\/CoreReport.html");
					print OUT "<HTML>\n";
				}
				if ($dir eq $nomodsum_dir){ # summary page
					open (OUT2, ">$nomodsum_dir\/CoreSummaryReport.html");
					print OUT2 "<HTML>\n";
					print OUT2 "<HEAD><Center><Font Size=6>Summary of All Samples<\/Font><\/Center><\/HEAD>\n";
				if (-e "$dir\/IDwG.html"){
					open (HTM, "<$dir\/IDwG.html");
				} elsif (-e "$dir\/IDwGmod.html"){
					open (HTM, "<$dir\/IDwGmod.html");
				}
				#	open (HTM, "<$dir\/IDwG.html");
					my ($body, $table) = (0, 0);
					while (<HTM>){
						last if (/\/HTML/);
						$body++ if (/BODY/);
						if (/\<\/TABLE/ && $table == 0){ $table++; next; }
						next if ($body < 1 || $table < 1);
						print OUT2 $_;
					}
					print OUT2 "<\/HTML>\n";
					print OUT "<\/HTML>\n";
				} else { # fraction page
					open (IND, ">$dir\/CoreReport.html");
					print IND "<HTML>\n";
					print IND "<HEAD><Center><Font Size=6>$group<\/Font><\/Center><\/HEAD>\n";
					print OUT "<HEAD><Center><Font Size=6>$group<\/Font><\/Center><\/HEAD>\n";
				if (-e "$dir\/IDwG.html"){
					open (HTM, "<$dir\/IDwG.html");
				} elsif (-e "$dir\/IDwGmod.html"){
					open (HTM, "<$dir\/IDwGmod.html");
				}
				#	open (HTM, "<$dir\/IDwG.html");
					my ($body, $table) = (0, 0);
					print "$dir\/IDwG.html\n";
					while (<HTM>){
						last if (/\/HTML/);
						$body++ if (/BODY/);
						if (/\<\/TABLE/ && $table == 0){ $table++; next; }
						next if ($body < 1 || $table < 1);
						print OUT $_;
						print IND $_;
					}
					print IND "<\/HTML>\n";
					print OUT "<BR><BR>\n";
				}
			} else { # single idsum, create report for single sample
				if (-e "$dir\/IDwG.html"){
					open (HTM, "<$dir\/IDwG.html");
				} elsif (-e "$dir\/IDwGmod.html"){
					open (HTM, "<$dir\/IDwGmod.html");
				}
				
				open (OUT, ">$dir\/CoreReport.html");
				print OUT "<HTML>\n";
				print OUT "<HEAD><Center><Font Size=6>$group<\/Font><\/Center><\/HEAD>\n";
				my ($body, $table) = (0, 0);
				while (<HTM>){
					last if (/\/HTML/);
					$body++ if (/BODY/);
					if (/\<\/TABLE/ && $table == 0){ $table++; next; }
					next if ($body < 1 || $table < 1);
					print OUT $_;
				}
				print OUT "<\/HTML>\n";
			}	
		}
	}
	############################################################################################
	# Publication Table and Link Printout
	############################################################################################
	#return if (!defined($nomodsum_dir));
	#if (($dir eq $nomodsum_dir) || ($dir eq $modsum_dir)){
		if (defined($paramhash{'publication_table'})){
			if ($paramhash{'publication_table'}>0){
				my $pubdir = "$dir\/publicationtable"; mkdir "$pubdir";
				my $specdir = "$pubdir\/spectra";
				if ($paramhash{'publication_table'} > 1) {
					system("rm -r $specdir") if (-e "$specdir"); mkdir "$specdir";
				}
				my $tabletxt = "$pubdir\/TableS.txt";
				open (OUT, ">$tabletxt");
				print OUT "Table S. Proteins\/Peptides Identified in Experiment\n";
				print OUT "Listed are total number of spectral counts \(\#SC\), the number of unique peptides assigned on each protein \(\#Pep\)";
				print OUT ", sequencing coverage \(\%\), the mass error measured in Orbitrap \(deltaMass\), the SEQUEST matching scores \(XCorr ";
				print OUT "and deltaCn\), and the link to assigned spectra. The labeled heavier K or R residues are marked.\n";
				print OUT "\#Accession;\#SC;\#Pep;Coverage\(\%\);Peptide Sequences;Tryptic;Unique;Position;Outfile;ExpMS1mz;deltaMass \(ppm\);Charge;XCorr;deltaCn;Link\n";
				my $totallinks = 0;
				for my $protein (sort {$$proteinhash{$b}{'coverage'}<=>$$proteinhash{$a}{'coverage'}} keys %$proteinhash){
					$totallinks += $$proteinhash{$protein}{'total'};
				}
				my $linknum = 0;
				for my $protein (sort {$$proteinhash{$b}{'coverage'}<=>$$proteinhash{$a}{'coverage'}} keys %$proteinhash){
					#if ($$proteinhash{$protein}{'total_nomod'} != $$proteinhash{$protein}{'total'}){ print Dumper($$proteinhash{$protein});exit; }
					for my $peptide (keys %{$$proteinhash{$protein}{'peptides'}}){
						#print Dumper($$peptidehash{$peptide});exit;
						my $best = find_best_outfile($$peptidehash{$peptide});
						%{$$proteinhash{$protein}{'peptides'}{$peptide}{'best_hash'}} = %{$$peptidehash{$peptide}{'outfiles'}{$best}};
						$$proteinhash{$protein}{'peptides'}{$peptide}{'best_out'} = $best;
					}
					for my $peptide (sort {$$proteinhash{$protein}{'peptides'}{$b}{'best_hash'}{'XCorr'}<=>
																 $$proteinhash{$protein}{'peptides'}{$a}{'best_hash'}{'XCorr'}
													 } keys %{$$proteinhash{$protein}{'peptides'}}){
						$linknum++;
						print "\rWorking on publication table: $linknum of $totallinks  ...            ";
						my %tmphash = %{$$peptidehash{$peptide}{'outfiles'}};
						for my $outfile (sort {$tmphash{$b}{'XCorr'}<=> $tmphash{$a}{'XCorr'}} keys %tmphash){
							my @parts = split('\.', $outfile);
							#print "$outfile\n";	print "@parts\n";	print Dumper($tmphash{$outfile});exit;
							my $best = $$proteinhash{$protein}{'peptides'}{$peptide}{'best_out'};
							my $peptidenum = $$proteinhash{$protein}{'total_nomod'};
							print OUT "$protein;$$proteinhash{$protein}{'occurrence'};$$proteinhash{$protein}{'total_nomod'};";
							printf OUT "%.2f;", $$proteinhash{$protein}{'coverage'};
							print OUT "$tmphash{$outfile}{'peptide'};";
							print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'tryptic'};";
							print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'unique'};";
							print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'};";
							print OUT "$outfile;";
							printf OUT "%.4f;", (($tmphash{$outfile}{'MH'} - 1.00728)/$tmphash{$outfile}{'charge'})+1.00728;
							if (defined($tmphash{$outfile}{'ppm'})){
								printf OUT "%.2f;", $tmphash{$outfile}{'ppm'};
							} else {
								print OUT "0;";
							}
							print OUT "$tmphash{$outfile}{'charge'};";
							printf OUT "%.2f;", $tmphash{$outfile}{'XCorr'};
							printf OUT "%.2f;", $tmphash{$outfile}{'dCn'};
							my ($showpep, $pep) = get_linkpeptide($tmphash{$outfile}{'peptide'});
							my $link = $peptidenum.$showpep.$parts[0].$parts[1];
							print OUT "spectra/$link.html\n";
							next if (-e "$specdir\/$link.html");
							my $string = get_modstring($pep, $tmphash{$outfile}{'run'});
							my $dta = $tmphash{$outfile}{'path'}; $dta =~ s/\.out\Z/\.dta/;
							my $poststring = "$showion?Dta=$dta&MassType=1&NumAxis=1&$string";
							if ($paramhash{'publication_table'} > 1){
								my @lines = qx(echo mod_perl rules | POST \'$poststring\');
								open (OUT2, ">$specdir\/$link.html");
								my $origgif = "";
								my $newgif = "$specdir\/$link.gif";
								for (@lines){
									if (/IMG\sSRC/){
									$_ =~ s/(\/tmp\/[0-9]+\.gif)/$link.gif/; # Use for local links
									$origgif = $1;
									print OUT2 $_;
								} else {
									print OUT2 $_;
								}
								}
								system("cp $origgif $newgif");
								system("chmod 775 $newgif");
							}
						}
					}
				}
			}
		}
	#}
	#
	#print publish table
	#foreach my $pro (keys %$grouphash) { print $pro,',', $grouphash->{$pro}->{group},';'; }
	rm_zombie($proteinhash,$grouphash);
	print "\nGenerating publication tables\n";
	#print "Generating publication tables\n";
	if (-e "$dir/publications") {}	else {system("mkdir $dir/publications");}
	pickBestScanForPeptide($peptidehash);
	printPublishedTable($proteinhash,$peptidehash,$grouphash,"$dir\/publications\/id_uni_prot.txt",1);
	printPublishedTable($proteinhash,$peptidehash,$grouphash,"$dir\/publications\/id_all_prot.txt",0);

	# print publish peptide table
	if (defined($mod))
	{
		printPublishedTable_peptide($proteinhash,$peptidehash,$grouphash,"$dir\/publications\/id_uni_pep.txt",1,$modifiedSites);
		printPublishedTable_peptide($proteinhash,$peptidehash,$grouphash,"$dir\/publications\/id_all_pep.txt",0,$modifiedSites);
	}
	else
	{
		printPublishedTable_peptide($proteinhash,$peptidehash,$grouphash,"$dir\/publications\/id_uni_pep.txt",1);
		printPublishedTable_peptide($proteinhash,$peptidehash,$grouphash,"$dir\/publications\/id_all_pep.txt",0);
	}

	if (-e "$dir/simplified_report") {} else {system("mkdir $dir/simplified_report");}
	if ( $paramhash{unique_protein_or_peptide} eq 'protein' )
	{
		system("cp $dir\/publications\/id_uni_prot.txt $dir/simplified_report");
	}
	else
	{
		system("cp $dir\/publications\/id_uni_pep.txt $dir/simplified_report");
	}

	# print 1% FDR table
	my (%peptidehash_1pct,%proteinhash_1pct);
	my ($grouphash_1pct)=build_1pct_hash($peptidehash,$proteinhash,\%dbhash,\%peptidehash_1pct,\%proteinhash_1pct);
	rm_zombie(\%proteinhash_1pct,$grouphash_1pct);
	pickBestScanForPeptide(\%peptidehash_1pct);
	printPublishedTable(\%proteinhash_1pct,\%peptidehash_1pct,$grouphash_1pct,"$dir\/publications\/id_uni_prot_1FDR.txt",1);
	printPublishedTable(\%proteinhash_1pct,\%peptidehash_1pct,$grouphash_1pct,"$dir\/publications\/id_all_prot_1FDR.txt",0);
	# print publish peptide table
	if (defined($mod))
	{
		printPublishedTable_peptide(\%proteinhash_1pct,\%peptidehash_1pct,$grouphash_1pct,"$dir\/publications\/id_uni_pep_1FDR.txt",1,$modifiedSites);
		printPublishedTable_peptide(\%proteinhash_1pct,\%peptidehash_1pct,$grouphash_1pct,"$dir\/publications\/id_all_pep_1FDR.txt",0,$modifiedSites);
	}
	else
	{
		printPublishedTable_peptide(\%proteinhash_1pct,\%peptidehash_1pct,$grouphash_1pct,"$dir\/publications\/id_uni_pep_1FDR.txt",1);
		printPublishedTable_peptide(\%proteinhash_1pct,\%peptidehash_1pct,$grouphash_1pct,"$dir\/publications\/id_all_pep_1FDR.txt",0);
	}

	# print 1% FDR databasae
	if (!(-e "$dir\/misc")) { system("mkdir $dir\/misc"); }
	print_DB(\%proteinhash_1pct,$grouphash_1pct,"$dir\/misc\/uni_prot_1FDR.fasta",1);
	print_DB(\%proteinhash_1pct,$grouphash_1pct,"$dir\/misc\/all_prot_1FDR.fasta",0);
}

sub get_modstring{
  my ($peptide, $run) = @_;
  # Get modification parameters for displayions cgi #
  my $modstring = "";
  while (my ($key, $value) = each %{$paramhash{$run}{'staticmods'}}){ $value += $residhash{$key}; $modstring .= "Mass$key=$value&"; }
  my $num = 0; my %mods;
  while (my ($key, $value) = each %{$paramhash{$run}{'dynamicmods'}}){  $num++; $modstring .= "DMass$num=$value&"; $mods{$key} = $num; }
  if ($num > 0){
    $modstring .= "DSite=";
    my @peparray = split("", $peptide);
    for (my $index = 0; $index < scalar(@peparray); $index++){
      next if ($peparray[$index] =~ /[\#\*\@\%\&\~\$\^]/); my $next = $index+1;
      if ($next >= scalar(@peparray)){ $modstring .= 0; last; }
      if ($peparray[$next] !~ /[\#\*\@\%\&\~\$\^]/){ $modstring .= 0;
      } else {
        if ($peparray[$index] =~ /[STY]/){
          $modstring .= $mods{'STY'};
        } else {
          $modstring .= $mods{$peparray[$index]};
        }
      }
    }
    my $nomodpep = $peptide; $nomodpep =~ s/[\#\@\%\&\~\$\^\*]//g;
    $modstring .= "&Pep=$nomodpep";
  }
  return $modstring;
}


sub save_Filtering
{
        my ($group_dir,$proteinhash,$FilteredHashRef,$txtfile)=@_;

        open (OUT, ">$group_dir/$txtfile");
        print OUT "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;unique;tryptic\n";
        for my $peptide (keys %$FilteredHashRef)
        {
                for my $outfile (keys %{$$FilteredHashRef{$peptide}{'outfiles'}})
                {
                        for my $protein (keys %{$$FilteredHashRef{$peptide}{'proteins'}})
                        {
                        #       print $protein,"\n";
                                $protein =~ s/\#\#//g;
                                next if (!defined($$proteinhash{$protein}));

				
			next if (!defined($peptide) || ! defined($outfile));

###################### xusheng on 8/30/2012 #####################
#				print $peptide,"\t",$outfile,"\n";
########################################

                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'peptide'};";
                                print OUT "$protein;$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'path'};";
                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'MH'};";
                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'expMH'};";
                                if (!defined($$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'ppm'}))
                                {
                                        print OUT "N/A;";
                                }
                                else
                                {
                                        print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'ppm'};";
                                }
                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'XCorr'};";
                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'dCn'};";
                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'ions'};";
                                print OUT "$$FilteredHashRef{$peptide}{'outfiles'}{$outfile}{'red'};";
                              #  print OUT "$$proteinhash{$protein}{'group'};";
                              #  print OUT "$$proteinhash{$protein}{'subgroup'};";
                                print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'unique'};";
                                print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'tryptic'};";
                                print OUT "\n";
                              #  if (!defined($$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'}))
                               # {
                                #        print "$peptide\n";
                                #        print Dumper($$proteinhash{$protein}{'peptides'}{$peptide});exit;
                               # }
                               # print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'}\n";
                        }
                }
        }
        close(OUT);


}


sub get_linkpeptide{
	my ($pep) = @_;
	$pep =~ s/[A-Z\-]\.//; $pep =~ s/\.[A-Z\-]//;
	my @array = split('', $pep);
	my $showpep = "";
	for (my $i=0; $i<scalar(@array); $i++){
		next if ($array[$i] =~ /[\#\@\%\&\~\$\^\*]/);
		my $j = $i+1;
		if (defined($array[$j])){
			if ($array[$j] =~ /[\#\@\%\&\~\$\^\*]/){ $showpep .= "$array[$i]_" }
			else { $showpep .= $array[$i]; }
		} else {
			$showpep .= $array[$i];
		}
	}
	return ($showpep, $pep);
}

sub find_best_outfile{
	my ($hash) = @_;
	my %temphash;
	for my $outfile (keys %{$$hash{'outfiles'}}){
		if ($$hash{'outfiles'}{$outfile}{'charge'} == 3){ 
			$temphash{'outfiles'}{$outfile}{'XCorr'} = $$hash{'outfiles'}{$outfile}{'XCorr'} - 1.5; 
		} elsif ($$hash{'outfiles'}{$outfile}{'charge'} == 4){ 
			$temphash{'outfiles'}{$outfile}{'XCorr'} = $$hash{'outfiles'}{$outfile}{'XCorr'} - 2.5; 
		} else {
			$temphash{'outfiles'}{$outfile}{'XCorr'} = $$hash{'outfiles'}{$outfile}{'XCorr'};
		}
	} 
	
	for my $outfile (sort {$temphash{'outfiles'}{$b}{'XCorr'}<=>$temphash{'outfiles'}{$a}{'XCorr'}}keys %{$temphash{'outfiles'}}){
		return $outfile;
	}
}

sub count_nomod{
	my ($proteinhash) = @_;

	# Added total_nomod, unique_nomod, and shared_nomod DMD 5/24/05
 	while (my ($protein, $prohash) = each %$proteinhash){
   	my (%pcount_total, %pcount_unique, %pcount_shared);
   	for my $peps (keys %{$$prohash{'peptides'}}){
     	my $origpep = $peps;
     	$peps =~ s/C[\*\@\%\&\~\$\^\#]/C/g;
     	$peps =~ s/M[\*\@\%\&\~\$\^\#]/M/g;
     	$pcount_total{$peps} = $peps if (!defined($pcount_total{$peps}));
     	if ($$prohash{'peptides'}{$origpep}{'unique'}){
       	$pcount_unique{$peps} = $peps if (!defined($pcount_unique{$peps}));
     	} else {
       	$pcount_shared{$peps} = $peps if (!defined($pcount_shared{$peps}));
     	}
   	}
   	$$prohash{'total_nomod'} = scalar(keys %pcount_total);
   	$$prohash{'unique_nomod'} = scalar(keys %pcount_unique);
   	$$prohash{'shared_nomod'} = scalar(keys %pcount_shared);
 	}
}

sub BypassFiltering{
	open (LOGFILE, ">>$log_file");

	my ($peptidehash, $proteinhash, $fprhash, $runhash, $group) = @_;
	print "\n";
	print "Creating proteinhash ...";
	print LOGFILE "\nCreating proteinhash ...";
	$idutils->create_proteinhash($proteinhash, $peptidehash, \%dbhash);
	print "\n";
	count_nomod(\%$proteinhash);
	my %tempfprhash;
	$idutils->count_fpr(\%$proteinhash, \%tempfprhash, $paramhash{'bypass_filtering'});
	if ($tempfprhash{'t'} == 0){
		print "No proteins were identified for this group ($group)!!!\n";exit;
		print LOGFILE "No proteins were identified for this group ($group)!!!\n";exit;
	}
	my ($rate1, $rate2, $rate3) = (0,0,0);
	$rate1 = ($tempfprhash{'r1'}/$tempfprhash{'t1'})*100 if ($tempfprhash{'t1'} != 0);
	$rate2 = ($tempfprhash{'r2'}/$tempfprhash{'t2'})*100 if ($tempfprhash{'t2'} != 0);
	$rate3 = ($tempfprhash{'r3'}/$tempfprhash{'t3'})*100 if ($tempfprhash{'t3'} != 0);
	printf "     1pep fpr = %.2f, 2pep fpr = %.2f, 3+pep fpr = %.2f\n", $rate1, $rate2, $rate3;
	printf LOGFILE "     1pep fpr = %.2f, 2pep fpr = %.2f, 3+pep fpr = %.2f\n", $rate1, $rate2, $rate3;

	%$fprhash = %tempfprhash;
	
	print "\n";
	if ($$fprhash{'t1'} != 0){
		printf "FPR of proteins matched with only 1 peptide:       %.2f%% out of $$fprhash{'t1'} proteins\n", ($$fprhash{'r1'}/$$fprhash{'t1'})*100;
		printf LOGFILE "FPR of proteins matched with only 1 peptide:       %.2f%% out of $$fprhash{'t1'} proteins\n", ($$fprhash{'r1'}/$$fprhash{'t1'})*100;
	}
	if ($$fprhash{'t2'} != 0){
		printf "FPR of proteins matched with 2 peptides:           %.2f%% out of $$fprhash{'t2'} proteins\n", ($$fprhash{'r2'}/$$fprhash{'t2'})*100;	
		printf LOGFILE "FPR of proteins matched with 2 peptides:           %.2f%% out of $$fprhash{'t2'} proteins\n", ($$fprhash{'r2'}/$$fprhash{'t2'})*100;
	} else {
		print "  No proteins with 2 peptides were identified!!!\n";
		print LOGFILE "  No proteins with 2 peptides were identified!!!\n";
	}
	if ($$fprhash{'t3'} != 0){
		printf "FPR of proteins matched with at least 3 peptides:  %.2f%% out of $$fprhash{'t3'} proteins\n", ($$fprhash{'r3'}/$$fprhash{'t3'})*100;
		printf LOGFILE "FPR of proteins matched with at least 3 peptides:  %.2f%% out of $$fprhash{'t3'} proteins\n", ($$fprhash{'r3'}/$$fprhash{'t3'})*100;
	} else {
		print "No proteins with 3 peptides were identified!!!\n";
		print LOGFILE "No proteins with 3 peptides were identified!!!\n";
	}
	my $total_decoy_number = ($$fprhash{'r1'}+$$fprhash{'r2'}+$$fprhash{'r3'});
	my $total_target_number = ($$fprhash{'t1'}+$$fprhash{'t2'}+$$fprhash{'t3'}) - ($$fprhash{'r1'}+$$fprhash{'r2'}+$$fprhash{'r3'});
#        printf "     Total protein FPR:  %.2f%%\n", (($$fprhash{'r1'}+$$fprhash{'r2'}+$$fprhash{'r3'})/(($$fprhash{'t1'}+$$fprhash{'t2'}+$$fprhash{'t3'})-($$fprhash{'r1'}+$$fprhash{'r2'}+$$fprhash{'r3'})))*100;

	printf "     Total protein FPR:  %.2f%%\n", $total_decoy_number/$total_target_number*100;
	printf LOGFILE "     Total protein FPR:  %.2f%%\n", $total_decoy_number/$total_target_number*100;
#        printf "     Total protein FPR:  %.2f%%\n", (($$fprhash{'r1'}+$$fprhash{'r2'}+$$fprhash{'r3'})/(($$fprhash{'t1'}+$$fprhash{'t2'}+$$fprhash{'t3'})-($$fprhash{'r1'}+$$fprhash{'r2'}+$$fprhash{'r3'})))*100;
	print "\n";
	
	my $outnum = 0;
	while (my ($peptide, $pephash) = each %$peptidehash){
		$outnum += scalar(keys %{$$pephash{'outfiles'}});
	}
	printf "Accepted %d proteins, %d peptides, and $outnum outfiles\n", scalar(keys %$proteinhash), scalar(keys %$peptidehash);
	printf LOGFILE "Accepted %d proteins, %d peptides, and $outnum outfiles\n", scalar(keys %$proteinhash), scalar(keys %$peptidehash);
#	printf "Total peptide FPR: %.2f",;
	print "\n";
}

sub XcorrDcnFiltering
{
	open (LOGFILE, ">>$log_file");

	my ($peptidehash, $proteinhash, $fprhash, $runhash, $group, $rerun, $featurehash, $printResults,$proteinQ,$acpT,$acpD) = @_;
	$proteinQ |= 0;
	my $epsilon=0.9;
	#$paramhash{'fpr_layer3'}
	
	my ($pepfpr,$peptide_decoy,$peptide_target,$adjust_pepfpr);
	# Loop for False Positive Rate 
	#my $min_fpr = $paramhash{'min_peptide_fpr'};
	my $min_fpr = $paramhash{'initial_outfile_fdr'};
	#my $min_fpr = 5;
	$paramhash{inflated_index}=1.5;
	my ($fpr3, $fpr2, $fpr1);
	my $uniqProFPR=100;
	my $tmpFolder='.idsumTmp';
	unless (-e $tmpFolder) { system("mkdir $tmpFolder"); }
	open(OUT,">\$tmpFolder\/Grouping_for_XCorr_filtering.txt"); close OUT;

	if ( !(defined($paramhash{FDR})) ) { die "FDR not defined!!!\n Cannot do FPR filtering.\n"; }
	if ( !(defined($paramhash{inflated_index})) ) { die "inflated_index not defined!!!\n Cannot do FPR filtering.\n"; }

	my %tmppeptidehash = %{clone($peptidehash)};#my %orig_peptidehash=%tmppeptidehash;
	#if ($printResults)
	if (0)
	{
		print "\nRunning score filtering for $group (grouping $paramhash{min_outfile_num_for_XCorr_filter} outfiles by peptide length, trypticity, mod, miscleavage, charge and dCn)\n";
		print LOGFILE "\nRunning score filtering for $group \n";
	}

	# count target and decoy numbers across all runs in a group
	my ($org_good,$org_bad)=scan_TD_num($peptidehash);
=head
	my ($org_good,$org_bad)=(0,0);
	foreach my $run (keys %{$runhash})
	{
		foreach my $outfile (keys %{$$runhash{$run}})
		{
			next if ($$runhash{$run}{$outfile}{status}<1);
			if ( $$runhash{$run}{$outfile}{protein} =~ /Decoy/ ) { $org_bad++; }
			else { $org_good++; }
		}
	}
=cut
	if ($printResults)
	{
	printf "\nStarting with %d outfiles: $org_good targets and $org_bad decoys (target outfiles FDR = %.2f%%)\n",$org_good+$org_bad,100*$org_bad/($org_good+0.01);
	printf LOGFILE "\nStarting with %d outfiles: $org_good targets and $org_bad decoys (target outfiles FDR = %.2f%%)\n",$org_good+$org_bad,100*$org_bad/($org_good+0.01);
	}
	$min_fpr = ($paramhash{'initial_outfile_fdr'}< 100*$org_bad/($org_good+0.01))? $paramhash{'initial_outfile_fdr'}:(100*$org_bad/($org_good+0.01));
	my $pre_min_fpr=$min_fpr;my $pre_over_fpr=0;my $over_filtering=0;
	my $loop_times=0;

	my ($scanFDRcut,$proteinFDR);
	my $orig_groupSize=$paramhash{min_outfile_num_for_XCorr_filter};
	while (1)
	{
		$loop_times++;
		#if ($loop_times>10) 
		if (0) 
		{ 
			if ($printResults) {print "Too many tries ... Cannot find optimal FPR cutoff\n"; last;} 
			else { $loop_times=0; $paramhash{FDR}=$min_fpr=0; next; }
		}
		
		# set group size = max(1/$min_fpr * 2,$orig_groupSize}
		my $allowedMaxDecoy_eachGroup=5;
		my $minGroupsize=int(100/($min_fpr+0.00001))*$allowedMaxDecoy_eachGroup;

		if ($paramhash{FDR_filtering_method} eq 'group_dyn') {
			$paramhash{min_outfile_num_for_XCorr_filter}=($minGroupsize>$orig_groupSize)?$minGroupsize:$orig_groupSize;
		}

		if ($printResults)
		{
		printf "\nTesting FDR for target outfiles (group size = %d): %.4f%% (= decoys/targets)\n",$paramhash{min_outfile_num_for_XCorr_filter},$min_fpr;
		printf LOGFILE "\nTesting FDR for target outfiles (group size = %d): %.4f%% (= decoys/targets)\n",$paramhash{min_outfile_num_for_XCorr_filter},$min_fpr;
		}
		#printf OUT "\nMinimum Peptide FPR (might take minutes): %.2f\n",$min_fpr;
		# XCorr and dCn filtering
		#my %temppeptidehash = %$peptidehash;my %orig_peptidehash=%temppeptidehash;
		#print "temppeptidehash:",scalar(keys %tmppeptidehash),"\n";

		$scanFDRcut=$min_fpr;
		($pepfpr,$peptide_decoy,$peptide_target,$adjust_pepfpr) = $idutils->find_XCorrdCn_cutoffs(\%tmppeptidehash, \%paramhash, $min_fpr, $featurehash, 0, $printResults);

		# record peptide FDR if 'group'
		#if ($paramhash{FDR_filtering_method} ne 'LDA')
		if (0)
		{
			foreach my $p (keys %tmppeptidehash)
			{
				if ( !defined($tmppeptidehash{$p}{Qvalue}) || $tmppeptidehash{$p}{Qvalue}>$adjust_pepfpr ) 
				{  
					$tmppeptidehash{$p}{Qvalue}=$adjust_pepfpr;
				}
			}
		}

		if ($printResults){
                printf "  Unique peptides: %d targets and %d decoys (FDR = %.2f%%)\n",$peptide_target,$peptide_decoy,$adjust_pepfpr;
                printf LOGFILE "  Unique peptides: %d targets and %d decoys (FDR = %.2f%%)\n",$peptide_target,$peptide_decoy,$adjust_pepfpr;
		}

		# Creation of protein hash
		my %tempproteinhash;
		$idutils->create_proteinhash(\%tempproteinhash, \%tmppeptidehash, \%dbhash);
		%$proteinhash = %tempproteinhash;
		count_nomod(\%$proteinhash);

		# protein grouping fo runiq protein FPR estimation
		#my ($grouphash, $group_num, $subgroup_num) = $gp->group_proteins(\%$proteinhash, \%tmppeptidehash, \%paramhash);
		#my ($grouphash, $decoy_grouphash);
		my ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num)=update_grouphash(\%$proteinhash, \%tmppeptidehash);
		#my %best_proteins;
                #assignBestProtein2Peptide(\%tmppeptidehash, $grouphash, $decoy_grouphash, \%best_proteins); #print scalar(keys %best_proteins),"\n";
                #assignOrder2group($grouphash, \%best_proteins);
                #assignOrder2group($decoy_grouphash, \%best_proteins);
		#my ($group_num, $subgroup_num)=get_group_num($grouphash);
		#my ($decoy_group_num, $decoy_subgroup_num)=get_group_num($decoy_grouphash);
=head
		my ($group_num, $subgroup_num, %acceptedGroup);
		$subgroup_num=0;
		foreach my $pro (keys %$grouphash) { if($grouphash->{$pro}->{'order'}==1) {$subgroup_num++; $acceptedGroup{$grouphash->{$pro}->{group}}=''; }} $group_num=scalar(keys %acceptedGroup);
		my ($decoy_group_num,$decoy_subgroup_num);
		$decoy_subgroup_num=0; undef(%acceptedGroup);
		foreach my $pro (keys %$decoy_grouphash) { if($decoy_grouphash->{$pro}->{'order'}==1) {$decoy_subgroup_num++; $acceptedGroup{$decoy_grouphash->{$pro}->{group}}='';  } } $decoy_group_num=scalar(keys %acceptedGroup);
		#print "$subgroup_num, $decoy_subgroup_num\n";
		#map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;
		#$group_num = $grouphash->{(sort {$grouphash->{$b}->{'group'} <=> $grouphash->{$a}->{'group'}} keys %$grouphash)[0]}{'group'};
=cut
		# Count the false positive rate for protein hash at 3 levels
		my %tempfprhash;
		$idutils->count_fpr(\%$proteinhash, \%tempfprhash, $paramhash{'bypass_filtering'},1-$printResults);
		if ($tempfprhash{'t'} == 0 || $subgroup_num==0){
			#print "No proteins were identified for this group ($group)!!!\n";exit;
			#print LOGFILE "No proteins were identified for this group ($group)!!!\n";exit;
		}
		#$uniqProFPR=$tempfprhash{'r'}*100/$subgroup_num;
		$uniqProFPR=$decoy_subgroup_num*100/($subgroup_num+0.01); $proteinFDR=$uniqProFPR;
		if ($printResults){
                printf "  Unique proteins: %d targets and %d decoys (FDR = %.2f%%)\n", $subgroup_num,$decoy_subgroup_num,$uniqProFPR;
                printf LOGFILE "  Unique proteins: %d targets and %d decoys (FDR = %.2f%%)\n", $subgroup_num,$decoy_subgroup_num,$uniqProFPR;
                printf "  Protein groups (genes): %d targets and %d decoys (FDR = %.2f%%)\n", $group_num,$decoy_group_num,$decoy_group_num*100/($group_num+0.01);
                printf LOGFILE "  Protein groups (genes): %d targets and %d decoys (FDR = %.2f%%)\n", $group_num,$decoy_group_num,$decoy_group_num*100/($group_num+0.01);
		}

		my $adjust_index;
		if ($min_fpr) 
		{
			if ( $paramhash{unique_protein_or_peptide} eq 'protein' ) { $adjust_index=$uniqProFPR/$min_fpr; }
			elsif ( $paramhash{unique_protein_or_peptide} eq 'peptide' ) { $adjust_index=$adjust_pepfpr/$min_fpr;  }
		}#$min_fpr=2*$paramhash{FDR}/$adjust_index;

		my $pass = 0;
		%$fprhash = %tempfprhash;
		my $targetFPR;
		if ( $paramhash{unique_protein_or_peptide} eq 'protein' ) { $targetFPR=$uniqProFPR;}
		elsif ( $paramhash{unique_protein_or_peptide} eq 'peptide' ) { $targetFPR=$adjust_pepfpr;}

		#print "$proteinQ, $targetFPR,  $paramhash{FDR}, $paramhash{unique_protein_or_peptide}, $paramhash{FDR_tolerance}\n";
		if ( defined($acpT) and defined($acpD) ) # for 1% protein FDR estimation
		{
			if ( 100*$acpD-(1+$epsilon)*$acpT < (1+$epsilon)*$subgroup_num-100*$decoy_subgroup_num ) {%$peptidehash=%{clone(\%tmppeptidehash)}; last;}
			else 
			{
				if ( $min_fpr>1 ) { $min_fpr-=1; }
				elsif ( $min_fpr>0.1 ) { $min_fpr-=0.1; }
				else { $min_fpr-=0.01; }
				next;
			}
		}
		else # normal run
		{
		# final check: 3 cases
		if ( ($proteinQ==0 and $targetFPR > $paramhash{FDR} + $paramhash{FDR}/10) 
           or ($proteinQ and $targetFPR > $paramhash{FDR} + $paramhash{FDR_tolerance}) )
		{
			#print "FDR too large: $targetFPR > $paramhash{FDR}\n";
			#print "peptidehash:",scalar(keys %$peptidehash),"\n";
			#print "temppeptidehash:",scalar(keys %tmppeptidehash),"\n";
			%$peptidehash=%{clone(\%tmppeptidehash)};

			$pre_min_fpr=$min_fpr;
			if ($min_fpr == 0 || $targetFPR==0){
				$fpr3 = "Not Achieved" if (!defined($fpr3));
				$fpr2 = "Not Achieved" if (!defined($fpr2));
				$fpr1 = "Not Achieved" if (!defined($fpr1));
				last;
			}
=head
			if ($over_filtering)
			{
				$min_fpr=($min_fpr+$pre_over_fpr)/2;
			}
			else
			{
				my $tmpFpr=$paramhash{inflated_index}*$paramhash{FDR}/$adjust_index;
				if ($tmpFpr>$min_fpr) { $min_fpr=$min_fpr * 0.9; }
				else { $min_fpr=$tmpFpr; }
			}
=cut
			if ( $min_fpr-1>1 ) { $min_fpr-=1; }
			#elsif ( $min_fpr-1 < 1) { $min_fpr=min(1,$min_fpr-0.1); }
			#elsif ( $min_fpr-1 < 1) { $min_fpr=(1<$min_fpr-0.1)?1:($min_fpr-0.1); }
			elsif ( $min_fpr > 1) { $min_fpr=1; }
			elsif ( $min_fpr>0.2 ) { $min_fpr-=0.1; }
			elsif ( $min_fpr>0.1 ) { $min_fpr=0.1; }
			else { $min_fpr-=0.01; }

			next;
		}
=head
		elsif ( $proteinQ and $paramhash{FDR}-$paramhash{FDR_tolerance}>$targetFPR)
		{
			# if over at the intial setting, exit
			#if ($min_fpr == $paramhash{min_peptide_fpr}) { %$peptidehash=%{clone(\%tmppeptidehash)}; last; }
			if ($min_fpr == $paramhash{initial_outfile_fdr}) { %$peptidehash=%{clone(\%tmppeptidehash)}; last; }
			$over_filtering=1;
			#print "peptidehash:",scalar(keys %$peptidehash),"\n";
			#print "temppeptidehash:",scalar(keys %tmppeptidehash),"\n";

			%tmppeptidehash=%{clone($peptidehash)};
			#print "peptidehash:",scalar(keys %$peptidehash),"\n";
			#print "temppeptidehash:",scalar(keys %tmppeptidehash),"\n";

			#print "pre_min_fpr: $pre_min_fpr; min_fpr: $min_fpr\n";
			$pre_over_fpr=$min_fpr;
			$min_fpr=($min_fpr+$pre_min_fpr)/2;
			next;
		}
=cut
		else {%$peptidehash=%{clone(\%tmppeptidehash)}; last;}
		}
	}
	# assign back orig_groupSize
	$paramhash{min_outfile_num_for_XCorr_filter}=$orig_groupSize;

	print "\n";
	#$peptidehash=\%temppeptidehash;
	#print "peptidehash:",scalar(keys %$peptidehash),"\n";
	#print "temppeptidehash:",scalar(keys %tmppeptidehash),"\n";
	#print "orig_peptidehash:",scalar(keys %orig_peptidehash),"\n";

	#%$peptidehash=%tmppeptidehash;
	#print "peptidehash:",scalar(keys %$peptidehash),"\n";

	# check whether decoys should be deleted (note that $min_fpr has been changed again after the last round of filtering)
############    dded by xusheng on 8/17/2012 ############################
=head
        if($paramhash{'keep_decoy'} eq '0')
        {
		($pepfpr,$peptide_decoy,$peptide_target,$adjust_pepfpr) = $idutils->find_XCorrdCn_cutoffs(\%$peptidehash, \%paramhash, $min_fpr, 1);
	}
	else
	{
               ($pepfpr,$peptide_decoy,$peptide_target,$adjust_pepfpr) = $idutils->find_XCorrdCn_cutoffs(\%$peptidehash, \%paramhash, $min_fpr, 0);
	}
=cut
#####################added by xusheng on 8/8/2012 #################
 #       my %temp_peptidehash = %$peptidehash;
####################################################################
	if ($printResults){
	my %mouthash;my %ttlPro;my %finalTDhash;$finalTDhash{T}=$finalTDhash{D}=0;
	while (my ($peptide, $pephash) = each %$peptidehash){
		for my $pro  (keys %{$$pephash{'proteins'}}){
			$ttlPro{$pro}++;
		}
		#$ttlPro+=scalar(keys %{$$pephash{proteins}});
		for my $outfile (keys %{$$pephash{'outfiles'}}){
			$mouthash{$outfile}++;
			if ( $$pephash{random} ) { $finalTDhash{D}++;}
			else { $finalTDhash{T}++;}
		}
	}
	#printf "Total outfiles in peptidehash after filtering: %d\n",scalar(keys %mouthash);
	#printf "Total proteins in peptidehash after filtering: %d\n",scalar(keys %ttlPro);
	printf "Recovery rate of expected true PSMs after running score filtering = ($finalTDhash{T}-$finalTDhash{D})/($org_good-$org_bad) = %.2f%%\n",($finalTDhash{T}-$finalTDhash{D})*100/($org_good-$org_bad+0.1);
	printf LOGFILE "\nRecovery rate of expected true PSMs after running score filtering = ($finalTDhash{T}-$finalTDhash{D})/($org_good-$org_bad) = %.2f%%\n",($finalTDhash{T}-$finalTDhash{D})*100/($org_good-$org_bad+0.1);

	if ($paramhash{'mass_accuracy'} == 1 && $rerun == 0){
		my $total = 0;
		my ($M0, $M1, $M2, $M3, $M4, $M6, $M7, $MX) = (0,0,0,0,0,0,0,0);
		for my $hash (values %$runhash){
			while (my ($outfile, $outhash) = each %$hash){
				last if (defined($$outhash{'skipma'}));
				next if (!defined($mouthash{$outfile}));
				if ($$outhash{'Mtype'} == -1){ $MX++;
				} elsif ($$outhash{'Mtype'} == 7){ $M7++;
				} elsif ($$outhash{'Mtype'} == 6){ $M6++;
				} elsif ($$outhash{'Mtype'} == 0){ $M0++;
				} elsif ($$outhash{'Mtype'} == 1){ $M1++;
				} elsif ($$outhash{'Mtype'} == 2){ $M2++;
				} elsif ($$outhash{'Mtype'} == 3){ $M3++;
				} elsif ($$outhash{'Mtype'} == 4){ $M4++;
				}
				$total++;
			}
		}
		#print "Distribution of sequenced isotopes: M-2=$M7 M-1=$M6 M=$M0, M+1=$M1, M+2=$M2, M+3=$M3, M+4=$M4\n" if ($total>0);
		#print LOGFILE "Distribution of sequenced isotopes: M-2=$M7 M-1=$M6 M=$M0, M+1=$M1, M+2=$M2, M+3=$M3, M+4=$M4\n" if ($total>0);
	}
	# Creation of protein hash
	my %tempproteinhash;
	$idutils->create_proteinhash(\%tempproteinhash, $peptidehash, \%dbhash);
	#my %temppeptidehash = %$peptidehash;
	#$idutils->create_proteinhash(\%tempproteinhash, \%temppeptidehash, \%dbhash);
	%$proteinhash = %tempproteinhash;
	count_nomod(\%$proteinhash);

	#print "Removing one-hit wonders if needed\n\n";
	#print LOGFILE "\nRemoving one-hit wonders if needed\n\n";
=head
	print "\n****************************************************\n";
	print "Recalculation of FDR after removing decoy peptides\n";
        print "****************************************************\n";
         print LOGFILE "\n****************************************************\n";
         print LOGFILE "Recalculation of FDR after removing decoy peptides\n";
         print LOGFILE "****************************************************\n";
=cut
	if (0){
	if ($$fprhash{'t1'} != 0){
                printf "FDR of total proteins matched with only 1 peptide:       %.2f%% out of $$fprhash{'t1'} proteins\n", ($$fprhash{'r1'}/($$fprhash{'t1'}))*100;
		printf LOGFILE "FDR of total proteins matched with only 1 peptide:       %.2f%% out of $$fprhash{'t1'} proteins\n", ($$fprhash{'r1'}/($$fprhash{'t1'}))*100;
	}
	if ($$fprhash{'t2'} != 0){
                printf "FDR of total proteins matched with 2 peptides:           %.2f%% out of $$fprhash{'t2'} proteins\n", ($$fprhash{'r2'}/($$fprhash{'t2'}))*100;
		printf LOGFILE "FDR of total proteins matched with 2 peptides:           %.2f%% out of $$fprhash{'t2'} proteins\n", ($$fprhash{'r2'}/($$fprhash{'t2'}))*100;
	} else {
		print "No proteins with 2 peptides were identified!!!\n";
		print LOGFILE "No proteins with 2 peptides were identified!!!\n";
	}
	if ($$fprhash{'t3'} != 0){
                printf "FDR of total proteins matched with at least 3 peptides:  %.2f%% out of $$fprhash{'t3'} proteins\n", ($$fprhash{'r3'}/($$fprhash{'t3'}))*100;
		printf LOGFILE "FDR of total proteins matched with at least 3 peptides:  %.2f%% out of $$fprhash{'t3'} proteins\n", ($$fprhash{'r3'}/($$fprhash{'t3'}))*100;
	} else {
		print "No proteins with 3 peptides were identified!!!\n";
		print LOGFILE "No proteins with 3 peptides were identified!!!\n";
	}
	}
=head
        printf "Total protein FPR:  %.2f%% ($$fprhash{'r'} decoy(s); $$fprhash{'t'} target(s))\n", ($$fprhash{'r'}/($$fprhash{'t'}))*100;
	printf "Total peptide FPR:  %.2f%% ($peptide_decoy decoy(s); $peptide_target target(s))\n", $adjust_pepfpr;
        printf LOGFILE "Total protein FPR:  %.2f%% ($$fprhash{'r'} decoy(s); $$fprhash{'t'} target(s))\n", ($$fprhash{'r'}/($$fprhash{'t'}))*100;
        printf LOGFILE "Total peptide FPR:  %.2f%% ($peptide_decoy decoy(s); $peptide_target target(s))\n", $adjust_pepfpr;
=cut
#	print "\n";
#	print "Peptide FPRs:\n";
###### changed by yanji
	#print "  For proteins matched with only 1 peptide:      $fpr1\n";
	#print "  For proteins matched with 2 peptides:          $fpr2\n";
	#print "  For proteins matched with at least 3 peptides: $fpr3\n";
###### end of change
#	printf "Accepted %d proteins and %d peptides       ", scalar(keys %$proteinhash), scalar(keys %$peptidehash);

############### added by xusheng on 8/8/2012 ##########################
#	my $accepted_peptide_num = scalar(keys %$peptidehash);
#	my $total_peptide_num = scalar(keys %temp_peptidehash);
#	printf "Total peptide FPR: %.2f%%\n",($total_peptide_num - $accepted_peptide_num)/$total_peptide_num*100;
###########################################################################
###### end of change
#        printf "Accepted %d proteins and %d peptides\n", scalar(keys %$proteinhash), scalar(keys %$peptidehash);
	#print "\n";
	}
	return ($scanFDRcut,$proteinFDR);
}

sub chi
{
        my ($c1,$c2,$t1,$t2) = @_;
        #print "($c1,$c2,$t1,$t2)\n";

        my $c = $c1 + $c2;
        my $t = $t1 + $t2;

        my $ave = ($c + $t)/2;

        my $x = 2*($c*(log($c/$ave)) + $t*(log($t/$ave)));
        #print "$x\n";

        my $pval = Statistics::Distributions::chisqrprob (1,$x);
}

sub help{
	print "\n";
  print "  USAGE: \n";
	#print "    idsum idsum.params (Please find latest parameter file in /spiders/samples_parameters/idsum/)\n\n";
	print "    $progname idsum.params\n\n (Latest parameter file can be found in /spiders/samples_parameters/idsum/)\n\n";
=head	print " idsum.params = contains parameters that are used in the program (default from /lab/user_data/params folder)\n\n";     
  print "  OPTIONS:\n";
	print "    -o  [Save Directory] -> set Output Directory, default is idsum.infodir\n\n";
	print "    -s  [Save Sum Directory] -> set Sum Output Directory, default is Sum_ID.infodir\n\n";
=cut
	exit;
}

sub xcorr_TD_dstr
{
	my ($runhash,$bottom,$roof,$step,$xcorr_dstr)=@_;

	for (my $i=0; $i<=($roof-$bottom)/$step; $i++) {$$xcorr_dstr{T}[$i]=$$xcorr_dstr{D}[$i]=0; }

	foreach my $outfile (keys %$runhash)
	{
		#print "$$runhash{$outfile}{XCorr}<=$bottom $$runhash{$outfile}{XCorr}>=$roof\n";
		next if ( $$runhash{$outfile}{XCorr}<=$bottom || $$runhash{$outfile}{XCorr}>=$roof );
		if ( $$runhash{$outfile}{protein} =~ /Decoy/ ) { $$xcorr_dstr{D}[int(($$runhash{$outfile}{XCorr}-$bottom)/$step)]++; }
		else  { $$xcorr_dstr{T}[int(($$runhash{$outfile}{XCorr}-$bottom)/$step)]++; }
	}

	for (my $i=1; $i<=($roof-$bottom)/$step; $i++)
	{
		$$xcorr_dstr{T}[$i]+=$$xcorr_dstr{T}[$i-1];
		$$xcorr_dstr{D}[$i]+=$$xcorr_dstr{D}[$i-1];
	}
}

sub rm_multiple_charge_for_single_precursor
{
	my ($runhash,$peptidehash,$del)=@_;
	$$del{multiCharge}{T}=$$del{multiCharge}{D}=0;

        # if there are multiple outfiles for one scan (due to un-predicted charge state), pick the best one
        #         # establish scanhash
        my %scanhash;
	while ( my ($outfile,$outhash) = each %{$runhash} )
	{
                        my $scan = $$outhash{'scan'};
#                        next if (!defined($outfile)||!defined($scan));
                        $scanhash{$scan}{$outfile}{'peptide'} = $$outhash{'intpep'};
                        $scanhash{$scan}{$outfile}{'charge'} = $$outhash{'charge'};
                        $scanhash{$scan}{$outfile}{'XCorr'} = $$outhash{'XCorr'};
        }
 	#if there are multiple outfiles (e.g., with different charges) for a scan, pick the best one
 	my $num = 0;
        for my $hash (values %scanhash){# %runhash{$scan}{$outfile}
        #while (my ($run,$runhash) = each %scanhash){# %runhash{$scan}{$outfile}
                #for my $hash (values %$runhash){# %hash{$outfile}
                        if (scalar(keys %$hash) > 1)# multiple outfiles for one scan: with different charges
                        {
                                my $bestXCorr = 0;
                                my $bestoutfile;
                                for my $outfile (keys %$hash)
                                {
                                        my $orig_XCorr = $$hash{$outfile}{'XCorr'};
=head
                                        if ($$hash{$outfile}{'charge'} == 3)
                                        {
                                                if ($bestXCorr<$orig_XCorr-1)
                                                {
                                                        $bestoutfile = $outfile;
                                                        $bestXCorr = $orig_XCorr-1;
                                                }
                                        }
                                        else
                                        {
=cut
                                                if ($bestXCorr <= $orig_XCorr)
                                                {
                                                        $bestoutfile = $outfile;
                                                        $bestXCorr = $orig_XCorr;
                                                }
                                        #}
                                }
                                for my $outfile (keys %$hash){
                                        next if ($outfile eq $bestoutfile);
                                        my $peptide = $$hash{$outfile}{'peptide'};
                                        $num++;

					if (!defined($$runhash{$outfile}{protein})) {die '?????',"$outfile,$$peptidehash{$peptide}{outfiles}{$outfile}{potein}\n";}
					if ( $$runhash{$outfile}{protein} =~ /Decoy/ ) { $$del{multiCharge}{D}++; }
					else  { $$del{multiCharge}{T}++; }

                                        delete $$peptidehash{$peptide}{'outfiles'}{$outfile};
					delete $$runhash{$outfile};
                                }
                        }#end of if
                #}
        }

}

sub rm_multi_monoisotopic_outfiles
{
	my ($runhash,$peptidehash,$del)=@_;
	$$del{multiMonoisot}{T}=$$del{multiMonoisot}{D}=0;

	my %scanhash;
	while ( my ($outfile,$outhash) = each %{$runhash} )
	{
		my $scan = $$outhash{'scan'};
		$scan =~ m/^(\d+)\.(\d+)$/;
		$scan=$1;
		my $ppi=$2;
		my @tmp=split(//,$ppi);
		$ppi=$tmp[0];

		$scanhash{$scan}{$ppi}{$outfile}{'peptide'} = $$outhash{'intpep'};
		$scanhash{$scan}{$ppi}{$outfile}{'XCorr'} = $$outhash{'XCorr'};
	}
	#open(OUT,">multiMonoisot_deleted");
	for my $ppihash (values %scanhash)
	{
		for my $hash (values %$ppihash)
		{
			if ( scalar(keys %$hash) > 1 )
			{
				my $bestXCorr = 0;
				my $bestoutfile;
				for my $outfile (keys %$hash)
				{
					my $orig_XCorr = $$hash{$outfile}{'XCorr'};
					if ($bestXCorr <= $orig_XCorr)
					{
						$bestoutfile = $outfile; $bestXCorr = $orig_XCorr;
					}
				}
				for my $outfile (keys %$hash)
				{
					#print OUT "$outfile\t$$hash{$outfile}{'XCorr'}\t$bestoutfile\n";
					next if ($outfile eq $bestoutfile);
					my $peptide = $$hash{$outfile}{'peptide'};

					if ( $$runhash{$outfile}{protein} =~ /Decoy/ ) { $$del{multiMonoisot}{D}++; }
					else { $$del{multiMonoisot}{T}++; }

					delete $$peptidehash{$peptide}{'outfiles'}{$outfile};
					delete $$runhash{$outfile};
				}
			}
		}
	}
	#close OUT;
}

sub MS2_scan_count
{
	my ($outfile,$ms2Scan)=@_;

	my $scan;
	if ($outfile =~ m/\.(\d+)\.(\d+)\.(\d)\.out$/) { $scan=$1; }
	elsif ($outfile =~ m/\.(\d+)\.(\d+)\.(\d)\.spout$/) { $scan=$1; }
	else {die "Wrong outfile format: $outfile!!! \n";}

	if (!defined($$ms2Scan{$scan})) { $$ms2Scan{$scan}=''; }
} 

sub buildAcceptedPSMpepXML
{
        my ($outInforHash,$peptidehash,$acpOut,$paraHash)=@_;

        # store accepted outfiles in a hash
        while ( my ($pep, $hash) = each %$peptidehash )
        {
                foreach my $out (keys %{$$hash{'outfiles'}})
                {
                        $out =~ s/\.out$//;
                        $out =~ s/\.spout$//;
                        $pepxmlParser->add_outfile($paraHash, $outInforHash,$acpOut,$out);
                }
        }
}

sub printFractionHash
{
	my ($acceptedOutHash,$totalParaHash,$group_dir,$paraHash)=@_;

	# build %fr2out
	my (%fr2out,%frouthash);
	foreach my $out (keys %{$acceptedOutHash})
	{
		#test01.4.1.2
		my @t=split(/\./,$out);
		my $fr=$t[0];
		$fr2out{$fr}{$out}='';
	}

	# print %fr2out
	if (!(-e "$group_dir/html"))  { system("mkdir $group_dir/html"); }
	if (!(-e "$group_dir/html/fractionHashes"))	{ system("mkdir $group_dir/html/fractionHashes"); }
	foreach my $fr (keys %fr2out)
	{
		foreach my $out (keys %{$fr2out{$fr}})
		{
			$pepxmlParser->add_outfile($totalParaHash, $acceptedOutHash,\%frouthash,$out);
		}
		#$pepxmlParser->printPepXML($totalParaHash,\%frouthash,"$group_dir/fractionHashes/$fr",10);
		store \%frouthash, "$group_dir/html/fractionHashes/$fr.hash";
	}
}

sub buildFeatureHash
{
	my ($peptidehash,$featurehash)=@_;

	while (my ($intpep, $intpephash) = each %$peptidehash)
	{
		$$intpephash{'random'}=1;
		while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}})
		{
			$$featurehash{$outfile}{'XCorr'}=$$outhash{'XCorr'};
                        $$featurehash{$outfile}{'dCn'}=$$outhash{'dCn'};
                        $$featurehash{$outfile}{'ppm'}=$$outhash{'ppm'};
                        #$$featurehash{$outfile}{'ppm_passed'}=$$outhash{'ppm'};
                        $$featurehash{$outfile}{'charge'}=$$outhash{'charge'};
                        $$featurehash{$outfile}{'mis'}=$$outhash{'mis'};
                        $$featurehash{$outfile}{'mod'}=$$outhash{'mod'};
                        $$featurehash{$outfile}{'tryptic'}=$$outhash{'tryptic'};

			my $nomod=$intpep;
			$nomod =~ s/[\@\%\&\~\$\^\#\*\^\~]//g;

                        $$featurehash{$outfile}{'peptide'}=$intpep;
                        $$featurehash{$outfile}{'peptideLength'}=length($nomod);

                        if ($$outhash{'protein'} =~ /Random__/ || $$outhash{'protein'} =~ /Decoy__/)
                        {
                                $$featurehash{$outfile}{'random'} = 1;
                        }
                        else
                        {
                                $$featurehash{$outfile}{'random'} = 0; $$intpephash{'random'}=0;
                        }

		}
	}
}

sub print_DB
{
	my ($proteinhash,$grouphash,$output,$uniqProtein)=@_;

	open(OUT,">$output");
	#foreach my $pro (keys %{$proteinhash})
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
                                                           $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
                                                           $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
                                                } keys %$grouphash)
	{
		if ($uniqProtein and $grouphash->{$protein}->{'order'}!=1) { next; }
		next if ($protein =~ /Decoy/);
		print OUT "\>$protein $$proteinhash{$protein}{annotation}\n$$proteinhash{$protein}{sequence}\n";
	}
	close OUT;
}


sub printPublishedTable
{
	my ($Protein_Hash,$peptidehash,$grouphash,$output,$uniqProtein)=@_;

	#foreach my $pro (keys %$grouphash) { print $pro,',', $grouphash->{$pro}->{group},';'; }

	open(OUT,">$output");
	my ($subgroup_num,$totalProtN)=(0,0);
	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;
	map {$totalProtN++ if($_ =! m/Decoy/)} keys %{$Protein_Hash};

	# unique or total?
	if ($uniqProtein)
	{
		print OUT "Unique proteins identified by mass spectrometry (n = $subgroup_num)\n";
	}
	else
	{
		#print OUT "All proteins identified by mass spectrometry (n = ",scalar(keys %{$Protein_Hash}),")\n";
		print OUT "All proteins identified by mass spectrometry (n = $totalProtN)\n";
	}

	# two search engines
	if ($paramhash{search_engine} eq 'sequest')
	{
		print OUT 'Protein Group#	Protein Accession #	Protein Description	GN	PSM#	Total Peptide#	Unique Peptide#	%Coverage 	Peptide of the Highest Score	Run#	Scan#	m/z	z	ppm	Xcorr	dCn	Q-value(%)';
	}
	elsif ($paramhash{search_engine} eq 'jump')
	{
		print OUT 'Protein Group#	Protein Accession #	Protein Description	GN	PSM#	Total Peptide#	Unique Peptide#	%Coverage 	Peptide of the Highest Score	Run#	Scan#	m/z	z	ppm	Jscore	dJn	Q-value(%)';
	}
	print OUT "$annotationColHeader\n";

	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
                                                           $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
                                                           $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
                                                } keys %$grouphash)
        {
	#print $protein,',', $grouphash->{$protein}->{group},';';
		if ($uniqProtein and $grouphash->{$protein}->{'order'}!=1) { next; }
		#my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
		#next if (!defined($grouphash->{$protein}->{'formal_group'}));
		next if ($protein =~ m/Decoy/);
		if (!defined($grouphash->{$protein}->{'formal_group'})) { print "formal_group not defined: $protein\n";next; }
		my $groupnum = $grouphash->{$protein}->{'formal_group'};
		if (!defined($$Protein_Hash{$protein}{'occurrence'})) {die "zombie protein: $protein";}
                my $description = $$Protein_Hash{$protein}{'annotation'};
                my $SC = $$Protein_Hash{$protein}{'occurrence'};
		my $GN = $$Protein_Hash{$protein}{'GN'};
		my $outfile = $$Protein_Hash{$protein}{best_outfile};
		# test01.2471.1.2.out
		$outfile =~ m/(.*?)\.(\d+)\.(\d+)\.(\d)/; my ($run,$scan,$charge)=($1,$2,$4);
		my $xcorr=$$peptidehash{$$Protein_Hash{$protein}{best_peptide}}{outfiles}{$outfile}{XCorr};
		my $dCn=$$peptidehash{$$Protein_Hash{$protein}{best_peptide}}{outfiles}{$outfile}{dCn};
		my $ppm=$$peptidehash{$$Protein_Hash{$protein}{best_peptide}}{outfiles}{$outfile}{ppm};
		#my $mz=$$peptidehash{$$Protein_Hash{$protein}{best_peptide}}{outfiles}{$outfile}{MH}/$charge;
		my $mz=($$peptidehash{$$Protein_Hash{$protein}{best_peptide}}{outfiles}{$outfile}{MH}-$Hydrogen_mass)/$charge+$Hydrogen_mass;
		
		## Modified by JCho on 02/03/2015 ##

		#printf OUT "$groupnum\t$protein\t$description\t$GN\t$SC\t$$Protein_Hash{$protein}{'total'}\t$$Protein_Hash{$protein}{'unique'}\t%.2f\t",$$Protein_Hash{$protein}{'coverage'};
		#print OUT "$$Protein_Hash{$protein}{best_peptide}\t";
		#printf OUT "$run\t$scan\t%.5f\t$charge\t%.4f\t$xcorr\t$dCn",$mz,$ppm;
		#printf OUT "$groupnum\t$protein\t$description\t$GN\t$SC\t$$Protein_Hash{$protein}{'total'}\t$$Protein_Hash{$protein}{'unique'}\t%.2f\t$$Protein_Hash{$protein}{best_peptide}\t$run\t$scan\t%.5f\t$charge\t%.4f\t$xcorr\t$dCn",$$Protein_Hash{$protein}{'coverage'},$mz,$ppm;
		
		my $nTotPeptide = $$Protein_Hash{$protein}{'total'};
		my $nUniqPeptide = $$Protein_Hash{$protein}{'unique'};
		my $coverage = sprintf("%.2f", $$Protein_Hash{$protein}{'coverage'});
		my $bestPeptide = $$Protein_Hash{$protein}{best_peptide};
		$mz = sprintf("%.5f", $mz);
		$ppm = sprintf("%.4f", $ppm);
		print OUT "$groupnum\t$protein\t$description\t$GN\t$SC\t$nTotPeptide\t$nUniqPeptide\t$coverage\t$bestPeptide\t$run\t$scan\t$mz\t$charge\t$ppm\t$xcorr\t$dCn";

		## End of modification		

		if (defined($$Protein_Hash{$protein}{Qvalue})) 
		{
			printf OUT "\t%.4f",$$Protein_Hash{$protein}{Qvalue};
		}
		else { print OUT "\t",'NA'; }

		print OUT "\t",$grouphash->{$protein}->{'annotationCols'};
		print OUT "\n";
	}

	close OUT;
}

sub printPublishedTable_peptide
{
	my ($Protein_Hash,$peptidehash,$grouphash,$output,$uniqProtein,$modifiedSites)=@_;

	# count unique protein numbers (target proteins, no decoy)
	open(OUT,">$output");
	my ($subgroup_num);
	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;

	# count target peptide numbers
	my $targetPep=0;
	foreach my $pep (keys %{$peptidehash}) { if ($$peptidehash{$pep}{random}==0) {$targetPep++;} }

	if (defined($modifiedSites))
	{
		# unique or total?
		if ($uniqProtein)
		{
			print OUT "Unique modified peptides identified by mass spectrometry\n";
			#print OUT "n = ",scalar(keys %{$peptidehash})," modified peptides\n";
			print OUT "n = ",$targetPep," modified peptides\n";
			print OUT "n = ",$modifiedSites," modified sites\n";
			print OUT "n = ",$subgroup_num," modified proteins\n";
		}
		else
		{
			print OUT "Total modified peptides identified by mass spectrometry\n";
			#print OUT "n = ",scalar(keys %{$peptidehash})," modified peptides\n";
			print OUT "n = ",$targetPep," modified peptides\n";
			print OUT "n = ",$modifiedSites," modified sites\n";
			print OUT "n = ",scalar(keys %{$Protein_Hash})," modified proteins\n";
		}

		# two search engines
		if ($paramhash{search_engine} eq 'sequest')
		{
			print OUT 'Peptides	Protein Group#	Protein Accession #	Protein Description	GN	PSM#	Run#	Scan#	m/z	z	ppm	Xcorr	dCn	Q-value(%)',"\n";
		}
		elsif ($paramhash{search_engine} eq 'jump')
		{
			print OUT 'Peptides	Protein Group#	Protein Accession #	Protein Description	GN	PSM#	Run#	Scan#	m/z	z	ppm	Jscore	dJn	Q-value(%)',"\n";
		}
	}
	else
	{
		# unique or total?
		if ($uniqProtein)
		{
			print OUT "Unique peptides identified by mass spectrometry\n";
			#print OUT "n = ",scalar(keys %{$peptidehash})," peptides\n";
			print OUT "n = ",$targetPep," peptides\n";
			print OUT "n = ",$subgroup_num," proteins\n";
		}
		else
		{
			print OUT "Total peptides identified by mass spectrometry\n";
			#print OUT "n = ",scalar(keys %{$peptidehash})," peptides\n";
			print OUT "n = ",$targetPep," peptides\n";
			print OUT "n = ",scalar(keys %{$Protein_Hash})," proteins\n";
		}

		# two search engines
		if ($paramhash{search_engine} eq 'sequest')
		{
			print OUT 'Peptides	Protein Group#	Protein Accession #	Protein Description	GN	PSM#	Run#	Scan#	m/z	z	ppm	Xcorr	dCn	Q-value(%)',"\n";
		}
		elsif ($paramhash{search_engine} eq 'jump')
		{
			print OUT 'Peptides	Protein Group#	Protein Accession #	Protein Description	GN	PSM#	Run#	Scan#	m/z	z	ppm	Jscore	dJn	Q-value(%)',"\n";
		}
	}

	foreach my $pep (keys %{$peptidehash})
	{
		next if ($$peptidehash{$pep}{random}==1);
		my $SC = scalar(keys %{$$peptidehash{$pep}{outfiles}});
		my $outfile = $$peptidehash{$pep}{best_outfile};
		$outfile =~ m/(.*?)\.(\d+)\.(\d+)\.(\d)/; my ($run,$scan,$charge)=($1,$2,$4);
		my $xcorr=$$peptidehash{$pep}{outfiles}{$outfile}{XCorr};
		my $dCn=$$peptidehash{$pep}{outfiles}{$outfile}{dCn};
		my $ppm=$$peptidehash{$pep}{outfiles}{$outfile}{ppm};
		if ($charge==0) {die "zero charge:$outfile\n";}
		#my $mz=$$peptidehash{$pep}{outfiles}{$outfile}{MH}/$charge;
		my $mz=($$peptidehash{$pep}{outfiles}{$outfile}{MH}-$Hydrogen_mass)/$charge+$Hydrogen_mass;


		$mz = sprintf("%.5f", $mz);
		$ppm = sprintf("%.4f", $ppm);
		$xcorr = sprintf("%.2f", $xcorr);
		$dCn = sprintf("%.2f", $dCn);

		#foreach my $protein (keys %{$$peptidehash{$pep}{proteins}})
		foreach my $protein (sort {$grouphash->{$a}->{'full_group'}<=>$grouphash->{$b}->{'full_group'}} keys %{$$peptidehash{$pep}{proteins}})
		{
			#my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
			my $groupnum = $grouphash->{$protein}->{'formal_group'};
			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $GN = $$Protein_Hash{$protein}{'GN'};
			
			#print OUT "$pep\t$groupnum\t$protein\t$description\t$GN\t$SC\t$run\t$scan\t$mz\t$charge\t$ppm\t$xcorr\t$dCn";
			print OUT "$$peptidehash{$pep}{orig_peptide}\t$groupnum\t$protein\t$description\t$GN\t$SC\t$run\t$scan\t$mz\t$charge\t$ppm\t$xcorr\t$dCn";
			if ( defined($$peptidehash{$pep}{Qvalue}) ) 
			{ 
				#print OUT "\t$$peptidehash{$pep}{Qvalue}"; 
				my $qvalue=$$peptidehash{$pep}{Qvalue};
				$qvalue = sprintf("%.3f", $qvalue); 
				print OUT "\t$qvalue"; 
			}
			else { print OUT "\t",'NA'; }
			if (defined($modifiedSites) and  defined($$peptidehash{$pep}{outfiles}{$outfile}{PTM_sites}))
			{ print OUT "\t$$peptidehash{$pep}{outfiles}{$outfile}{PTM_sites}\t$$peptidehash{$pep}{outfiles}{$outfile}{PTM_score}";  }
			#else { print OUT "\t",'NA'; }
			print OUT "\n";
			if ($uniqProtein ) { last; }
		}
	}
	close OUT;
}

sub pickBestScanForPeptide
{
	my ($peptidehash)=@_;

	foreach my $pep (keys %{$peptidehash})
	{
		my ($maxXCorr,$bestoutfile)=(0,'');
		foreach my $outfile (keys %{$$peptidehash{$pep}{outfiles}})
		{
			if ($maxXCorr<$$peptidehash{$pep}{outfiles}{$outfile}{XCorr})
			{
				$maxXCorr=$$peptidehash{$pep}{outfiles}{$outfile}{XCorr};
				$bestoutfile=$outfile;
			}
		}
		$$peptidehash{$pep}{best_outfile}=$bestoutfile;
	}

}

sub assignBestProtein2Peptide
{
	my ($peptidehash, $grouphash, $decoy_grouphash, $best_proteins)=@_;

	#print scalar(keys %{$best_proteins}),"\n";
	my (%t,%d);
	while(my ($pep, $hash) = each %{$peptidehash})
	{
		my $best_protein='';
		if ( $$hash{random}==0 ) # for targets
		{
			foreach my $protein (sort {$grouphash->{$a}->{'full_group'}<=>$grouphash->{$b}->{'full_group'} } 
						keys %{$$hash{proteins}} )
			{ 
				$best_protein=$protein;  last; 
			}
			$t{$best_protein}='';
		}
		else			# for decoys
		{
			foreach my $protein (sort {$decoy_grouphash->{$a}->{'full_group'}<=>$decoy_grouphash->{$b}->{'full_group'} } 
						keys %{$$hash{proteins}} )
			{ 
				#print "$protein,";
				$best_protein=$protein;  last; 
			}
			$d{$best_protein}=''
		}
		if ($best_protein eq '') { die "Empty best_protein: $pep, $best_protein\n"; }
		$$hash{best_protein}=$best_protein;
		$$best_proteins{$best_protein}='';
	}
	#print scalar(keys %{$best_proteins}),"\n";
	#print scalar(keys %t),',',scalar(keys %d),"\n";
	my $noPit=0; my $a;
	foreach  $a (keys %t) 
	{ 
		if (!defined($grouphash->{$a}) and !defined($decoy_grouphash->{$a})) 
		{ 
			$noPit++;
			#print "Warning: your search database and .pit file could be unmatched (e.g.,$a not found in .pit file)!!!\n"; 
		} 
	}
	if ($noPit) 
	{ 
		print "Warning: $noPit proteins not found in .pit file (e.g.,$a)!!!\n";  
	}
}

sub assignOrder2group
{
	my ($proteinhash, $hash, $best_proteins)=@_;

	my $k=0;
	foreach my $protein (keys %$hash)
	{
		if (defined($$best_proteins{$protein})) { $hash->{$protein}->{order}=1; $k++;}#if ($protein =~ /Decoy/) {print "$protein,";} }
		elsif (defined($$proteinhash{$protein})) { $hash->{$protein}->{order}=2; }
		else { delete $hash->{$protein}; }
	}
	#print "$k\n";
}

sub get_group_num
{
	my ($grouphash)=@_;
        my ($group_num, $subgroup_num, %acceptedGroup);

       	$subgroup_num=0;
       	foreach my $pro (keys %$grouphash) 
	{ 
		if($grouphash->{$pro}->{'order'}==1) 
		{
			$subgroup_num++; 
			my $group=$grouphash->{$pro}->{group};#print "$pro,$group;";
			$acceptedGroup{$group}=''; 
		}
	} 
	$group_num=scalar(keys %acceptedGroup);

	return ($group_num, $subgroup_num);
}

sub update_grouphash
{
	my ($proteinhash, $peptidehash)=@_;

	#my ($grouphash,$decoy_grouphash)=build_grouphash("HUMAN.pit");
	#my ($grouphash,$decoy_grouphash)=build_grouphash("HUMAN_083014.pit");
	if (!defined($paramhash{pit_file})) { die "pit file not defined!!!\nPlease update your parameter file.\n"; }
	my ($grouphash,$decoy_grouphash)=build_grouphash($paramhash{pit_file}); # build from .pit file
        my %best_proteins;
        assignBestProtein2Peptide($peptidehash, $grouphash, $decoy_grouphash, \%best_proteins); # assign Best Protein to Peptide: %peptidehash{$pep}{best_protein}
        assignOrder2group($proteinhash, $grouphash, \%best_proteins); 				# for best protein:  %grouphash->{$protein}->{order}=1
												# for other proteins: %grouphash->{$protein}->{order}=2
        assignOrder2group($proteinhash, $decoy_grouphash, \%best_proteins);			# same for decoy grouphash
        my ($group_num, $subgroup_num)=get_group_num($grouphash); 	# numbers of groups
        my ($decoy_group_num, $decoy_subgroup_num)=get_group_num($decoy_grouphash); # numbers of groups
 
	return ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num);
}

sub build_grouphash
{
	my ($pitfile)=@_;
	# pit initialization:
	my $use_pit=1;
	my ($grouphash,$decoy_grouphash);
	#open(PIT,"HUMAN.pit");
	open(PIT,$pitfile) || die "Cannot open pit file ($pitfile)!!!\n";

	# header
	my $line=<PIT>;
	chomp($line);
	unless ( $line =~ m/^UniprotAC\tSJPGnumber\tGroupName\tAbundance\tResidueLength\tProteinName\tFullDescription(.*?)$/ ) 
	{ die "Unexpected .pit file header!!!\n$line\n"; }
	$annotationColHeader=$1;

	# content
	while(<PIT>)
	{
		chomp;
		#next if (/^UniprotAC/);
		#UniprotAC       SJPGnumber      GroupName       Abundance(emPAI)        ResidueLength   ProteinName     FullDescription TFBSmotif       Kinase  Oncogene        GPCR    EpigeneticRegulator
		#CON_A2A4G1      SJPG000001.001  CON_A2A4G1      -       454     co|CON_A2A4G1|Con       TREMBL:A2A4G1|Contaminant Tax_Id=10090 Gene_Symbol=Krt15 Keratin complex 1, acidic, gene 15     -       -       -       -       -
		unless (/^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)$/)
		{ die "Unexpected .pit file format!!!\n$_\n"; }
		my ($ac,$fmGroup,$gene,$abd,$l,$protein,$desc,$annotationCols)=($1,$2,$3,$4,$5,$6,$7,$8);

		#my ($ac,$fmGroup,$gene,$abd,$l,$protein,$desc,$TFBSmotif,$Kinase,$Oncogene)=split(/\t/,$_);
		#my ($ac,$fmGroup,$gene,$abd,$l,$TFBSmotif,$Kinase,$Oncogene,$protein,$desc)=split(/\t/,$_);

		my ($group,$subgroup)=split(/\./,$fmGroup);
		$subgroup =~ s/^0*//;
		#SJPG18913
		$group =~ m/^[A-Za-z]+(\d+)$/;
		$group = $1;
		$group =~ s/^0*//;


		# prioritize_contaminants = 0: lower contaminant priority
		if ($paramhash{prioritize_contaminants}==0 and $protein =~ m/$CONTAMINAT_PATTERN/)
		{
			$group += $LARGE_GROUP_OFFSET;
		}

		my $full_group="$group\.$subgroup";#print "$protein $full_group;";

		#my $fullgroup=$fmGroup;
		#$fullgroup =~ s///;

		if ($protein =~ m/Decoy__/)
		{
			$protein =~ s/\#\#//;
			#$ac =~ s/\#\#//;
			$decoy_grouphash->{$protein}->{abundance}=$abd;
			$decoy_grouphash->{$protein}->{gene}=$gene;
			$decoy_grouphash->{$protein}->{formal_group}=$fmGroup;
			$decoy_grouphash->{$protein}->{group}=$group;
			$decoy_grouphash->{$protein}->{subgroup}=$subgroup;
			$decoy_grouphash->{$protein}->{full_group}=$full_group;
			$decoy_grouphash->{$protein}->{annotationCols}=$annotationCols;
		}
		else
		{
			$grouphash->{$protein}->{abundance}=$abd;
			$grouphash->{$protein}->{gene}=$gene;
			$grouphash->{$protein}->{formal_group}=$fmGroup;
			$grouphash->{$protein}->{group}=$group;
			$grouphash->{$protein}->{subgroup}=$subgroup;
			$grouphash->{$protein}->{full_group}=$full_group;
			$grouphash->{$protein}->{annotationCols}=$annotationCols;
		}
	}
	close PIT;

	return ($grouphash,$decoy_grouphash,$annotationColHeader);
}

sub UpdateProteinQvalue
{
	my ($proteinhash,$proteinhash1,$qvalue)=@_;

	foreach my $pro (keys %{$proteinhash1})
	{
		if ( !defined($$proteinhash{$pro}) ) { next;}#die "$pro not exist in proteinhash!!! (UpdateProteinQvalue)\n"; }
		if ( !defined($$proteinhash{$pro}{Qvalue}) || $$proteinhash{$pro}{Qvalue}>$qvalue )
		{
			$$proteinhash{$pro}{Qvalue}=$qvalue;
		}
	}
}

sub UpdateProteinQvalue_pephashVersion
{
	my ($proteinhash,$peptidehash,$qvalue)=@_;

	foreach my $pep (keys %{$peptidehash})
	{
		foreach my $pro (keys %{$$peptidehash{$pep}{proteins}})
		{
			if ( !defined($$proteinhash{$pro}) ) { next;}
			if ( !defined($$proteinhash{$pro}{Qvalue}) || $$proteinhash{$pro}{Qvalue}>$qvalue )
			{ $$proteinhash{$pro}{Qvalue}=$qvalue; }
		}
	}
}

sub trimPeptideHash
{
	my ($scanFDRcut1, $peptidehash1, $featurehash1, $dbhash, $proteinhash1)=@_;

	$idutils->find_XCorrdCn_cutoffs($peptidehash1, \%paramhash, $scanFDRcut1, $featurehash1, 0);
	$idutils->create_proteinhash($proteinhash1, $peptidehash1, $dbhash);
}

sub splitpeptidehash
{
	my ($peptidehash,$proteinhash,$pep1hash,$pep2hash,$pep3hash)=@_;

	#print "peptidehash in splitpeptidehash: ",scalar(keys %$peptidehash),"\n";
	foreach my $pep (keys %{$peptidehash})
	{
		my $maxProPep=0;
		foreach my $pro (keys %{$$peptidehash{$pep}{proteins}})
		{
			if ($pro =~ m/\#\#/) { $pro =~ s/\#\#//; }
			if (defined( $$proteinhash{$pro} ))
			{
				if ( $maxProPep < scalar(keys %{$$proteinhash{$pro}{peptides}}) ) 
				{
					$maxProPep = scalar(keys %{$$proteinhash{$pro}{peptides}});
				}
			}
			else { print "not existed protein $pro in proteinhash (splitpeptidehash)!!!\n"; last; }
		}

		#print "$pep\t$maxProPep\n";
		if ( $maxProPep>=3 ) {$$pep3hash{$pep}=clone($$peptidehash{$pep});}
		elsif ( $maxProPep==2 ) {$$pep2hash{$pep}=clone($$peptidehash{$pep});}
		elsif ( $maxProPep==1 ) {$$pep1hash{$pep}=clone($$peptidehash{$pep});}
		else { print "Unexpected maxProPep: $maxProPep (splitpeptidehash)!!!$pep\n"; next;}
	}

	#print "pep1hash: ",scalar(keys %{$pep1hash}),"\n"; 
	my $scFDR1=Calculate_category_specific_FDR($pep1hash,\%paramhash);
	#print "pep2hash: ",scalar(keys %{$pep2hash}),"\n";
	my $scFDR2=Calculate_category_specific_FDR($pep2hash,\%paramhash);
	#print "pep3hash: ",scalar(keys %{$pep3hash}),"\n";
	my $scFDR3=Calculate_category_specific_FDR($pep3hash,\%paramhash);

	return ($scFDR1,$scFDR2,$scFDR3);
}

sub splitpeptidehash_multiPPI
{
	my ($peptidehash,$pep1hash,$pep2hash,$pep0hash,$pep0_1hash,$pep0_2hash)=@_;

	# split into two groups:
	# %pep2hash: 2+ precursor ions
	# %pep1hash: 1 precursor ion
	foreach my $pep (keys %{$peptidehash})
	{
		my %charges;
		#my $maxPPI=0;
		foreach my $out (keys %{$$peptidehash{$pep}{outfiles}})
		{
			#t01.116989.1.2.spout
			my ($fr,$scan,$ppi,$charge,$sf)=split /\./,$out;
			#my $charge="$ppi,$charge";
			$charges{$charge}=1;
		}

		if (scalar(keys %charges)>=2) {$$pep2hash{$pep}=clone($$peptidehash{$pep});}
		elsif (scalar(keys %charges)>=1) {$$pep1hash{$pep}=clone($$peptidehash{$pep});}
		else { print "Zero PPI (splitpeptidehash_multiPPI)!!!$pep\n"; next;}
	}

	# Count modified forms #: build %nomod{$nomod}{$peptide} = 1 #
	my %nomodhash;
	foreach my $pep (keys %{$peptidehash})
	{
		my $nomod = $pep; $nomod =~ s/[\#\*\@\%\&\~\$\^]//g;
		$nomodhash{$nomod}{$pep}=1;
	}

	# split %pep1hash into two groups:
	# %pep1hash: with multiple modified forms
	# %pep0hash: single form
	foreach my $pep (keys %{$pep1hash})
	{
		my $nomod = $pep; $nomod =~ s/[\#\*\@\%\&\~\$\^]//g;
		if ( scalar(keys %{$nomodhash{$nomod}})>=2 ) {}
		else # %pep0hash: single form
		{
			$$pep0hash{$pep}=clone($$pep1hash{$pep});
			delete $$pep1hash{$pep};
		}
	}

	# split %pep0hash into two groups:
	# %pep0_2hash: with multiple scans
	# %pep0_1hash: with single scan
	foreach my $pep (keys %{$pep0hash})
	{
		if (scalar(keys %{$$peptidehash{$pep}{outfiles}})>=2) { $$pep0_2hash{$pep}=clone($$pep0hash{$pep}); }
		else { $$pep0_1hash{$pep}=clone($$pep0hash{$pep}); }
	}
}

sub Calculate_category_specific_FDR 
{
	my ($pep1hash,$paramhash)=@_;
	my ($T,$D,$fdr);
	my $inflateIndex=2;

	# pep1hash
	#print "pep1hash:";
	# scan FDR
	($T,$D)=scan_TD_num($pep1hash);
	#printf "  Outfiles: $T targets and $D decoys (FDR = %.2f%%)\n",$D*100/$T;
	my $scanFDR=$D*100/($T+0.01);
	# peptide FDR
	($T,$D)=peptide_TD_num($pep1hash);
	#printf "  Unique peptides: $T targets and $D decoys (FDR = %.2f%%)\n",$D*100/$T;
	if ( $$paramhash{unique_protein_or_peptide} eq 'peptide' ) { $fdr=$D*100/($T+0.01);  }
	# protein FDR
	($T,$D)=protein_TD_num($pep1hash);
	#printf "  Unique proteins: $T targets and $D decoys (FDR = %.2f%%)\n",$D*100/$T;
	if ( $$paramhash{unique_protein_or_peptide} eq 'protein' ) { $fdr=$D*100/($T+0.01);  }
	# Calculate scan FDR cut for next step
	my $scFDR1=$$paramhash{FDR}*$scanFDR/($fdr+0.001);
	#my $scFDR1=$$paramhash{FDR}*$scanFDR*$inflateIndex/$fdr;

	return $scFDR1;
}

sub scan_TD_num
{
	my ($peptidehash)=@_;
	my ($T,$D)=(0,0);
	foreach my $pep (keys %{$peptidehash})
	{
		if ($$peptidehash{$pep}{random}==1) { $D+=scalar(keys %{$$peptidehash{$pep}{outfiles}}); }
		else { $T+=scalar(keys %{$$peptidehash{$pep}{outfiles}}); }
	}
	return ($T,$D);
}

sub peptide_TD_num
{
	my ($peptidehash)=@_;
	my ($T,$D)=(0,0);
	foreach my $pep (keys %{$peptidehash})
	{
		if ($$peptidehash{$pep}{random}==1) { $D++; }
		else { $T++; }
	}
	return ($T,$D);
}

sub protein_TD_num
{
	my ($peptidehash)=@_;
	my (%proteinhash);
	$idutils->create_proteinhash(\%proteinhash, $peptidehash, \%dbhash);
	my ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num)=update_grouphash(\%proteinhash, $peptidehash);
	return ($subgroup_num,$decoy_subgroup_num);
}

sub mergePeptidehashes
{
	my $peptidehash=shift(@_);
	undef(%{$peptidehash});
	foreach my $hash (@_)
	{
		foreach my $p (keys %{$hash})
		{
			$$peptidehash{$p}=clone($$hash{$p});
		}
	}
}

sub pseudo_protein_Qvalue
{
	my ($proteinhash,$proteinhash1,$peptidehash1,$featurehash1,$paramhash,$scanFDRcut1,$proteinFDR1)=@_;

        $idutils->find_XCorrdCn_cutoffs($peptidehash1, $paramhash, $scanFDRcut1, $featurehash1, 0, 0);
	undef(%$proteinhash1);
        $idutils->create_proteinhash($proteinhash1, $peptidehash1, \%dbhash);
	UpdateProteinQvalue($proteinhash,$proteinhash1,$proteinFDR1);
}

sub print_IDtxt
{
  my ($peptidehash,$proteinhash,$grouphash,$FDR_type,$QvalueCut,$output,$db)=@_;

  open (OUT, ">$output");
	print OUT "Database=$db\n";
  
###### changed by yanji, add rention time, RT, and intensity in the header line of ID.txt
  #print OUT "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;mis_cleavage;pos;precursor_peak_intensity_percentage\n";
  print OUT "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;pos;precursor_peak_intensity_percentage\n";
 # print OUT "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;RT;intensity(totIonCurrent);group;subgroup;unique;tryptic;pos\n";
###### end of change

  for my $peptide (keys %$peptidehash){
    if ( $FDR_type eq 'peptide' and $$peptidehash{$peptide}{Qvalue}>$QvalueCut) {next;}
    for my $outfile (keys %{$$peptidehash{$peptide}{'outfiles'}}){
      for my $protein (keys %{$$peptidehash{$peptide}{'proteins'}}){
				next if (!defined($$proteinhash{$protein}));
	if ( $FDR_type eq 'protein' and $$proteinhash{$protein}{Qvalue}>$QvalueCut) {next;}
				#print Dumper($$proteinhash{$protein}{'peptides'}{$peptide});exit;
				print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'peptide'};";
	  		print OUT "$protein;$$peptidehash{$peptide}{'outfiles'}{$outfile}{'path'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'MH'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'expMH'};";
				if (!defined($$peptidehash{$peptide}{'outfiles'}{$outfile}{'ppm'})){
          print OUT "N/A;";
        } else {
          print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'ppm'};";
          #print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'rawppm'};"; # for ppm debug
        }
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'XCorr'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'dCn'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'ions'};";
	  		print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'red'};";
###### added by ynaji, print rention time and intensity
               #         print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'rt'};";
                #        print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'prec_int'};";
####### end of change
#
	
			#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'mis'};";
	  		#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'mod'};";
	  		#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'tryptic'};";
				#my $nomod = $peptide; $nomod =~ s/[\@\%\&\~\$\^\*\#]//g; 
	  		#printf OUT "%d;", length($nomod);
#	  		print OUT "$$proteinhash{$protein}{'group'};";
#	  		print OUT "$$proteinhash{$protein}{'subgroup'};";
              print OUT "$grouphash->{$protein}->{'group'};";
              print OUT "$grouphash->{$protein}->{'subgroup'};";

				print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'unique'};";
				print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'tryptic'};";
				#print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'mis'};";
				if (!defined($$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'})){
					print "$peptide\n";
					print Dumper($$proteinhash{$protein}{'peptides'}{$peptide});exit;
				}
				#print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'}\n";
				print OUT "$$proteinhash{$protein}{'peptides'}{$peptide}{'peppos'};";
				if (defined($$peptidehash{$peptide}{'outfiles'}{$outfile}{'precursor_peak_intensity_percentage'})) { print OUT "$$peptidehash{$peptide}{'outfiles'}{$outfile}{'precursor_peak_intensity_percentage'}"; }
				else { print OUT "NA"; }
				print OUT "\n";
			}
    }
  }
  close OUT;
}

sub printFilteringResults
{
	my ($peptidehash,$proteinhash,$fprhash,$org_good,$org_bad)=@_;

	my ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num)=update_grouphash($proteinhash, $peptidehash);
	$idutils->count_fpr($proteinhash, $fprhash, $paramhash{'bypass_filtering'},1);
	my ($final_good,$final_bad)=scan_TD_num($peptidehash);
	printf "Accepting %d outfiles: $final_good targets and $final_bad decoys (FDR = %.2f%%)\n",$final_good+$final_bad,$final_bad*100/$final_good;
	printf LOGFILE"Accepting %d outfiles: $final_good targets and $final_bad decoys (FDR = %.2f%%)\n",$final_good+$final_bad,$final_bad*100/$final_good;
	my ($T,$D)=peptide_TD_num($peptidehash);
	printf "  Unique peptides: $T targets and $D decoys (FDR = %.2f%%)\n",$D*100/$T;
	printf LOGFILE"  Unique peptides: $T targets and $D decoys (FDR = %.2f%%)\n",$D*100/$T;
	undef $T; undef $D;
	printf "  Unique proteins: $subgroup_num targets and $decoy_subgroup_num decoys (FDR = %.2f%%)\n",$decoy_subgroup_num*100/$subgroup_num;
	printf LOGFILE "  Unique proteins: $subgroup_num targets and $decoy_subgroup_num decoys (FDR = %.2f%%)\n",$decoy_subgroup_num*100/$subgroup_num;
	my $proteinFDR=$decoy_subgroup_num*100/$subgroup_num;

	printf "  Protein groups (genes): $group_num targets and $decoy_group_num decoys (FDR = %.2f%%)\n", $decoy_group_num*100/$group_num;
	printf LOGFILE "  Protein groups (genes): $group_num targets and $decoy_group_num decoys (FDR = %.2f%%)\n", $decoy_group_num*100/$group_num;

	#Calculate_category_specific_FDR(\%peptidehash,\%paramhash);

	# recovery rate
	#printf "For FDR filtering section:\n  Recovery rate of expected true PSMs after running score filtering = ($final_good-$final_bad)/($org_good-$org_bad)=%.2f%%\n",($final_good-$final_bad)*100/($org_good-$org_bad+0.0001);
	printf "  Recovery rate of expected true PSMs after running score filtering = ($final_good-$final_bad)/($org_good-$org_bad)=%.2f%%\n",($final_good-$final_bad)*100/($org_good-$org_bad+0.0001);
	printf LOGFILE"  Recovery rate of expected true PSMs after running score filtering = ($final_good-$final_bad)/($org_good-$org_bad)=%.2f%%\n",($final_good-$final_bad)*100/($org_good-$org_bad+0.0001);

	return $proteinFDR;
}

sub build_1pct_hash
{
	my ($peptidehash,$proteinhash,$dbhash,$peptidehash_1pct,$proteinhash_1pct)=@_;

	if ($paramhash{unique_protein_or_peptide} eq 'protein')
	{
	foreach my $pep (keys %{$peptidehash})
	{
		my $proteinQ=100;
		foreach my $pro (keys %{$$peptidehash{$pep}{proteins}})
		{
			if (!defined($$proteinhash{$pro})) { next; }
			if (defined($$proteinhash{$pro}{Qvalue}) and  
			$proteinQ>$$proteinhash{$pro}{Qvalue} )
			{ $proteinQ=$$proteinhash{$pro}{Qvalue}; }
		}
		#$peptidehash{$pep}{proteinQ}=$proteinQ;
		if ( $proteinQ<=1 )
		{ $$peptidehash_1pct{$pep}=clone($$peptidehash{$pep}); }

	}
	}
	else
	{
	foreach my $pep (keys %{$peptidehash})
	{
		if ( $$peptidehash{$pep}{Qvalue}<=1 )
		{ $$peptidehash_1pct{$pep}=clone($$peptidehash{$pep}); }
	}
	}

	$idutils->create_proteinhash($proteinhash_1pct, $peptidehash_1pct, $dbhash);
	my ($grouphash, $decoy_grouphash, $group_num, $subgroup_num, $decoy_group_num, $decoy_subgroup_num)=update_grouphash($proteinhash_1pct, $peptidehash_1pct);

	return $grouphash;
}

=head
sub FDR_filtering
{
	while (1)
	{
		$idutils->find_XCorrdCn_cutoffs(\%tmppeptidehash, \%paramhash, $min_fpr, $featurehash, 0, $printResults);
	}
}
=cut
