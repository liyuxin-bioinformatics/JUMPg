#!/usr/bin/perl  -I /data1/pipeline/release/version12.1.0/JUMPsj/

## Release date: 7/01/2016
## Release version: version 12.1.3
## Fixed unsearched dtas 

use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";

use Cwd;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::Params;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::Decharge;
use Spiders::BuildIndex;
use Spiders::Job;
use Spiders::Error;
use Spiders::Digestion;
use Spiders::MakeDB;
use Spiders::MassAccuracy;
use Spiders::PIP;
use Spiders::Path;
use Spiders::RankHits;
use Spiders::SpoutParser;
use Parallel::ForkManager;
use Spiders::MassCorrection;

my $library = $Bin;
my $VERSION = 12.1.3;

my $progname = $0;
# remove path from our name
$progname =~ s@(.*)/@@i;

my ($help,$parameter,$raw_file);
GetOptions('-help|h'=>\$help,
			'-p=s'=>\$parameter,
		);

usage() if ($help || !defined($parameter));
#check_input($raw_file,\$parameter);

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump search using JUMP                  ****     # 
#       ****  Version 12.1.3	                      ****     #
#       ****  Copyright (C) 2012 - 2015               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################

EOF

#print "\n  Initializing jump program\n";
#print "  Features of JUMP search:\n";
#print "  1. ppm can be used for fragment tolerance\n";
#print "  2. ppi can be used to increase ID, but is not recommended for TMT\n";
#print "  3. turning off secord_search can substantially boost searching speed\n";
#print "\n";
 
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 


######## predefine the input #####################
my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();
$params->{'Mn_Mn1'} = 0.5;
$params->{'M_M1'} = 0.3;
$params->{'low_mass'} = 57;
$params->{'high_mass'} = 187;
$params->{'tag_coeff_evalue'} = 1;
$params->{'pho_neutral_loss'} = 0;
#$params->{'cluster'} = 1;
#$params->{'Job_Management_System'} = SGE;
	
database_creation();

## Create the path for multiple raw files
my %rawfile_hash;
print "  Using the following rawfiles:\n";
foreach $arg (sort @ARGV)
{
    my @suffixlist=();
    push @suffixlist,".raw";
    push @suffixlist,".RAW";
    push @suffixlist,".mzXML";

    if($arg=~/.[raw|RAW|mzXML]/)
    {
		print "  $arg","\n";
	}
}

open(LOC,">\.jump_s_tmp");
foreach $arg (sort @ARGV)
{
    my @suffixlist=();
    push @suffixlist,".raw";
    push @suffixlist,".RAW";
    push @suffixlist,".mzXML";
    if($arg=~/.[raw|RAW|mzXML]/)
    {	
		my ($filename, $directory, $suffix) = fileparse($arg,@suffixlist);	
        system(qq(mkdir $directory/$filename >/dev/null 2>&1));
        system(qq(mv $arg $directory/$filename >/dev/null 2>&1));
        my $datafile = "$directory/$filename";
        my $path = new Spiders::Path($datafile);


        my $list = $path->make_directory_list();
        if(@$list)
        {
		if (!defined($params->{'fast_run'}) or $params->{'fast_run'}==0) {
			$newdir = $path->choose_dir_entry($list,"  Choose a .out file directory",$newdir);
		} else {
			$newdir = $path->choose_default_entry($list,"  Choose a .out file directory",$newdir);
		}
        }
        else
        {
			$newdir	= $filename . ".1";
		if (!defined($params->{'fast_run'}) or $params->{'fast_run'}==0) {
			$newdir = $path->ask("  Choose a .out file directory",$newdir);
		}
        }
        print "  Using: $newdir\n";
        $path->add_subdir($newdir);
        my $dir =  $path->basedir() . "/$newdir";
		my $rawfile = "$datafile/$arg";
		$rawfile_hash{$rawfile} = $dir;
 
	system(qq(cp -rf $parameter "$dir/jump.params" >/dev/null 2>&1));
	print LOC "$dir\n";
    }
}
close LOC;
if (!defined($params->{'fast_run'}) or $params->{'fast_run'}==0) {
	system(qq(rm .jump_s_tmp));
}


foreach my $raw_file(sort keys %rawfile_hash)
{		
#	my $raw_file=$ARGV[0];

###### Get working directory #####################
	print "  Searching data: $raw_file\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	print "  Start: ";
	printf "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;	

	
	my $curr_dir = getcwd;
	my $dta_path = $rawfile_hash{$raw_file};
	  
	if($raw_file =~/\.mzXML/)
	{
		$raw_file =~s/\.mzXML//g;
	}

	######### programming starting information #############

	my $proc_raw = new Spiders::ProcessingRAW();
	$proc_raw->set_raw_file($raw_file);
#	my $dta_path = $proc_raw->get_rawfile_basename();

	## Modified by JCho on 01/15/2015
	## For the use of ascore, the generation of sequest.params is re-activated
#	$p->generateSeqeust($params,'sequest.params');
#	system(qq(mv sequest.params $dta_path >/dev/null 2>&1));
	$p->generateSeqeust($params,'sequest.params');
	system(qq(mv sequest.params $dta_path >/dev/null 2>&1));
	## End of modification
	
	###################### window part ##########################
	print "  Converting .raw into .mzXML file\n";

	my $mzXML = $proc_raw->raw2mzXML();
	##################### Linux part ############################
	print "  Extracting peaks from .mzXML\n";
	my $proc_xml = new Spiders::ProcessingMzXML();
	$proc_xml ->set_dta_path($dta_path);

	$proc_xml ->set_mzXML_file($mzXML);

	################### preprocessing #########################
	my (%ms_hash,%msms_hash,@mz_array);


	$proc_xml ->set_parameter($params);
	$proc_xml->generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);
	my $ms1N = scalar(keys %{$ms_hash{'surveyhash'}});
	my $ms2N = scalar(keys %msms_hash)-scalar(keys %{$ms_hash{'surveyhash'}});
	
	printf "\n  There are %d MS and %d MS/MS in the entire run\n", scalar(keys %{$ms_hash{'surveyhash'}}) , scalar(keys %msms_hash)-scalar(keys %{$ms_hash{'surveyhash'}});
	print "\n  Mass correction\n";
	my $masscorr = new Spiders::MassCorrection();
	my ($msms_hash_corrected,$mz_array_corrected) = $masscorr->massCorrection(\%ms_hash, \%msms_hash, \@mz_array, $params);
	%msms_hash = %$msms_hash_corrected;
	@mz_array = @$mz_array_corrected;

	print "\n  Decharging scans\n";
	my $pip = new Spiders::PIP;
	$pip->set_parameter($params);
	$pip->set_origmz_array(\@mz_array);
	$pip->set_origmsms_hash(\%msms_hash);
	$pip->set_isolation_window($params->{'isolation_window'});
	if(!defined($params->{'isolation_window_offset'}))
	{
		$params->{'isolation_window_offset'} = 0;
		print "  Warning: you did not define isolation_window_offset parameter, use default value of 0\n";
	}
	if(!defined($params->{'isolation_window_variation'}))
	{
		$params->{'isolation_window_variation'} = 0;
		print "  Warning: you did not define isolation_window_variation parameter, use default value of 0\n";		
	}
	
	$pip->set_isolation_offset($params->{'isolation_window_offset'});
	$pip->set_isolation_variation($params->{'isolation_window_variation'});
	
	$pip->set_dta_path($dta_path);	
	my $PIPref = $pip->Calculate_PIP();

	my ($charge_dist,$ppi_dist) = $pip->changeMH_folder($PIPref);


	########################## Start Searching #######################################
	print "  Starting database searching\n";

	my $job = new Spiders::Job;
	$job->set_library_path($library);		
	$job->set_dta_path("$dta_path");
	$job->set_pip($PIPref);
	my @file_array = glob("$dta_path/*.dta");
	my $random = int(rand(100));
	if($params->{'second_search'} == 0)
        {
		$job->create_script(0);
	}
	elsif($params->{'second_search'} == 1)
	{
		$job->create_script(1);
	}
	else
	{
		print "  Please specify a right second_search parameter!!\n";
	}
	my $temp_file_array=runjobs(\@file_array,$dta_path,"sch_${random}",$params->{'processors_used'});
	my $rerunN=0;
	my $orig_cluster=$params->{'cluster'};
	while(scalar(@{$temp_file_array})>0 and $rerunN<3)
	{
		$rerunN++;
		print "\n",scalar(@$temp_file_array)," .dta files not finished! Doing re-search (rerunN = $rerunN)\n";
		$params->{'cluster'} = 0;
		my $remaining=runjobs($temp_file_array,$dta_path,"rescue_$rerunN",$params->{'processors_used'});

		$temp_file_array=$remaining;
	}
	$params->{'cluster'} = $orig_cluster;
	if(scalar(@{$temp_file_array})>0)
	{
		#die "\n",scalar(@$temp_file_array)," .dta files not finished!\n";
		print "\nWarning: ",scalar(@$temp_file_array)," .dta files not finished!\n";
	}
		
	my $Rankhit = new Spiders::RankHits();
	my $p = Spiders::Params->new('-path'=>$parameter);
	my $params=$p->parse_param();
	$Rankhit->set_parameter($params);
#		my ($miscleavage_modification_coeff,$modification_coeff,$miscleavage_coeff) = $Rankhit->get_db_mis_mod();
	my $mainhash	= $Rankhit->parse_spOut_files_v5($dta_path);
	my ($sum_target,$sum_decoy,$cutoff_score) = $Rankhit->calculate_FDR($mainhash,0.01);
	print "\n  $sum_target targets and $sum_decoy decoys passed FDR = 1%\n";
=head
	if($params->{'second_search'} == 1)
	{
		print "\n  Research by adjusting mono precursors\n";	
		my ($sum_target,$sum_decoy,$cutoff_score) = $Rankhit->calculate_FDR($mainhash,0);		
		print "  using Jscore $cutoff_score at FDR = 0 for second search\n";
		my $research_scans = $Rankhit->get_score_scans($mainhash,$cutoff_score,$dta_path);
		system(qq(rm -f "$dta_path/runsearch_shell.pl"));
		system(qq(rm -f "$dta_path/runsearch_shell.pl"));
		
		$job->create_script(1);	
		$random = int(rand(100));	
		runjobs($research_scans,$dta_path,"resch_${random}");
		my $mainhash	= $Rankhit->parse_spOut_files_v5($dta_path);
		my ($sum_target,$sum_decoy,$cutoff_score) = $Rankhit->calculate_FDR($mainhash,0.01);
		print "\n  $sum_target targets and $sum_decoy decoys passed FDR = 1% after second search\n";	

		system(qq(rm -f "$summary_out"));
		system(qq(rm -f "$summary_out"));
		system(qq(rm -f "$summary_tag"));
		system(qq(rm -f "$summary_tag"));
	}
=cut
	print "  Generating dtas file\n";
	system(qq(cat $dta_path/*.dtas > $dta_path.dtas));      # make one dtas
	system(qq(rm $dta_path/*.dtas));        # delete job-specific dtas
	print "  Generating tags file\n";
	make_tags($dta_path,"$dta_path\.tags");
	print "  Generating pepXML file\n";
	my $spout=Spiders::SpoutParser->new();

	my $outXML = $dta_path . ".pepXML"; 
	my (%seqParaHash, %outInforHash);
	$spout->readSpsearchParam(\%seqParaHash,"$dta_path/jump.params");

	$spout->parseSpout(\%seqParaHash, \%outInforHash, "$dta_path");
	$spout->printPepXML(\%seqParaHash, \%outInforHash, "$dta_path", $outXML, 5);


	print "\n  Generating summary XLSX report file\n";
	$spout->set_MS1scanNumber($ms1N);
	$spout->set_MS2scanNumber($ms2N);
	$spout->set_totalScan($ms1N+$ms2N);
	$spout->set_chargeDistribution($charge_dist);
	$spout->set_ppiDistribution($ppi_dist);

	my $report_xls = $dta_path . ".xlsx";
	my $topMS2_array=topMS2scans(\%msms_hash);	
	$spout->printXLSX(\%seqParaHash, \%outInforHash, $dta_path,$report_xls,$topMS2_array);

	if(defined($params->{'temp_file_removal'}) and $params->{'temp_file_removal'}==1)
	{
		delete_run_dir($dta_path);
	}
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	print "  Date: ";
	printf "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;	

}

print "  Search finish\n\n";
			
sub Create_Sort_BashFile
{
	
	my ($Params,$dta_path) = @_;	
	my $FileName = "$dta_path/sort_db_".$Params->{"range"}.".sh";	
	my $cmd = join(" ","perl $dta_path/Create_Partial_Idx.pl -m",$Params->{"range"},"-job_num",$Params->{"JobNum"},"-dta_path",$dta_path,"-database_path",$dta_path,"-mass_index",$Params->{"mass_index"},"-peptide_index", $Params->{"peptide_index"},"-protein_index",$Params->{"protein_index"},"-databasename",$Params->{"databasename"},"-num_pro_per_job",$Params->{"num_pro_per_job"},"-prot_job_num",$Params->{"prot_job_num"},"-digestion",$Params->{"digestion"},"\n");
	#my $cmd = "perl Create_Partial_Idx.pl -m ".$Params->{"range"}." -job_num ".$Params->{"JobNum"}." -dta_path ".$Params->{"dta_path"}." -curr_dir ".$Params->{"curr_dir"}." -mass_index ".$Params->{"mass_index"}." -peptide_index ".$Params->{"peptide_index"}." -protein_index ".$Params->{"protein_index"}."\n";
	
	LuchParallelJob($FileName,$cmd,$Params->{"GridType"},"sort_db_".$Params->{"range"},$dta_path);

}

sub runjobs
{
	my ($file_array,$dta_path,$job_name,$MAX_PROCESSES) = @_;
	my $curr_dir = getcwd;
	#my $MAX_PROCESSES = 32;	
	$MAX_PROCESSES = defined($MAX_PROCESSES)?$MAX_PROCESSES:4;

	my $dta_num_per_file = 10;
	my $job_num = int($#$file_array / $dta_num_per_file) + 1;
	## Set the maximum number of jobs to 4000
	if ($job_num > 4000) {
		$job_num = 4000;
		$dta_num_per_file = int($#$file_array / $job_num) + 1;
	}

	if($params->{'cluster'} eq '0')
	{
		$job_num = $MAX_PROCESSES;
		$dta_num_per_file = int($#$file_array / $job_num) + 1;
	}
	#my $dta_num_per_file = int($#$file_array / $job_num) + 1;
	#my $dta_num_per_file = 200/($num_dynamic_mod*2);
	#my $job_num=int($#file_array/$dta_num_per_file)+1;
	#$job_num = $#$file_array+1 if($#$file_array<$job_num);
 
	for($i = 0; $i < $job_num; $i++)
	{	
		if (($i * $dta_num_per_file) > $#$file_array) {
			$job_num = $i;
			last;
		}
		open(JOB,">$dta_path/${job_name}_${i}.sh") || die "can not open the job files\n";
		my $dta_file_temp="";
	#	my $dta_file_temp2="";
		my @dta_file_arrays=();
		my $multiple_jobs_num  = 0;
		for($j = 0; $j < $dta_num_per_file; $j++)
		{
			if(($i*$dta_num_per_file+$j)<=$#$file_array)
			{
				$dta_file_temp .= " $$file_array[$i*$dta_num_per_file+$j]";
				push (@dta_file_arrays, $$file_array[$i*$dta_num_per_file+$j])
	#			$dta_file_temp .= " $$file_array[$i+$job_num*$j]";
	#			$dta_file_temp2 .= " $$file_array[$i+$job_num*$j]";
				
				#$multiple_jobs_num++;
				#if($multiple_jobs_num>100)
				#{
				#	push (@dta_file_arrays,$$file_array[$i+$job_num*$j]);
				#	$dta_file_temp2="";
				#	$multiple_jobs_num=0;
				#}				
			}
		}
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q normal\n";
			print JOB "#BSUB -M 2000\n";
			print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
			print JOB "#BSUB -eo $dta_path/${job_name}_${i}.e\n";
			print JOB "#BSUB -oo $dta_path/${job_name}_${i}.o\n";
			print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $curr_dir/$parameter -dta_path $dta_path $dta_file_temp\n";		
		}
		elsif($params->{'Job_Management_System'} eq 'SGE')
		{
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N ${job_name}_${i}\n";
			print JOB "#\$ -e $dta_path/${job_name}_${i}.e\n";
			print JOB "#\$ -o $dta_path/${job_name}_${i}.o\n";
			
			
	#		print JOB "perl $dta_path/runsearch_shell.pl -param $curr_dir/$parameter -dta_path $dta_path $dta_file_temp\n";	
			foreach (@dta_file_arrays)
			{
				print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $curr_dir/$parameter -dta_path $dta_path $_\n";	
			}
		}
		else
		{
			print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $curr_dir/$parameter -dta_path $dta_path $dta_file_temp\n";	
		}
		
		close(JOB);
	}

	######### running jobs ######################### 
	my $job_list;
	if($params->{'cluster'} eq '1')
	{
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			for(my $i=0;$i<$job_num;$i++)
			{
				$command_line = qq(cd $dta_path && bsub <${job_name}_${i}.sh);
				my $job=qx[$command_line];
				chomp $job;
				my $job_id=0;
				if($job=~/Job \<(\d*)\> is/)
				{
					$job_id=$1;
				}
				$job_list->{$job_id}=1;
			}
		}
		elsif($params->{'Job_Management_System'} eq 'SGE')
		{

			for(my $i=0;$i<$job_num;$i++)
			{
				my $job_name = "${job_name}_${i}.sh";
				$command_line = qq(cd $dta_path && qsub -cwd $job_name);
				my $job=qx[$command_line];
				chomp $job;
				my $job_id=0;
				if($job=~/$job_name \<(\d*)\> is/)
				{
					$job_id=$1;
				}
				$job_list->{$job_id}=1;
				my $count = $i+1;
				print "\r  $count jobs were submitted";				
			}	
		}
		print "\n  You submitted $job_num jobs for database search\n";
		Check_Job_stat("${job_name}_",$job_num,$dta_path);		
	}
	elsif($params->{'cluster'} eq '0')
	{
	print "Only use headmaster to run jobs (MAX_PROCESSES = $MAX_PROCESSES)\n";
        my $pm = new Parallel::ForkManager($MAX_PROCESSES);
        for my $i ( 0 .. $MAX_PROCESSES )
        {
            $pm->start and next;
			my $job_name = "${job_name}_${i}.sh";			
            system("cd $dta_path && sh $job_name >/dev/null 2>&1");
			print "\r  $i jobs were submitted";	
			Check_Job_stat("${job_name}_",$job_num,$dta_path);				
			
            $pm->finish; # Terminates the child process
        }
	
        $pm->wait_all_children;		
	}
#=head
	#### checking whether all dta files has been searched  
	my $temp_file_array;
	for(my $k=0;$k<=$#$file_array;$k++)
	{
		my $data_file = $file_array->[$k];
		my $out_file = $data_file;
		#$out_file =~ s/\.dta$/\.out/;
		$out_file =~ s/\.dta$/\.spout/;
		if(!-e $out_file)
		{
			push (@{$temp_file_array},$data_file);
		}
	}
	return $temp_file_array;
=head
	if(scalar(@temp_file_array)>0)
	{
		print "\n",scalar(@temp_file_array)," .dta files not finished! Doing re-search\n";
		runjobs(\@temp_file_array,$dta_path,"rescue");
		undef @temp_file_array;
	        for(my $k=0;$k<=$#$file_array;$k++)
        	{
               	 	my $data_file = $file_array->[$k];
                	my $out_file = $data_file;
                	#$out_file =~ s/\.dta$/\.out/;
                	$out_file =~ s/\.dta$/\.spout/;
                	if(!-e $out_file)
                	{
                       		 push (@temp_file_array,$data_file);
                	}
        	}
#=head
		if(scalar(@temp_file_array)>0)
		{		
			die scalar(@temp_file_array)," .dta files not finished!\nPlease search again!!!\n";
		}
#=cut
	}
=cut
#=cut	
#	print "\n  $job_name finished\n\n";
}

sub LuchParallelJob{
	
	my($FileName,$cmd,$GridType,$outputName,$dta_path)= @_;
	
	open(JOB,">$FileName") || die "can not open $FileName\n";
	if($GridType eq 'LSF')
	{	
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q normal\n";
			print JOB "#BSUB -M 20000\n";
			print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
			
			print JOB "#BSUB -eo $dta_path/$outputName.e\n";
			print JOB "#BSUB -oo $dta_path/$outputName.o\n";
			print JOB $cmd;		
			close(JOB);
			system(qq(bsub <$FileName >/dev/null 2>&1));	
	}
	if($GridType eq 'SGE')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#\$ \-S /bin/bash\n";  #In our cluster this line is esential for executing some bash commands such as for
		print JOB "#\$ \-N $outputName\n";
		print JOB "#\$ \-e $dta_path/$outputName.e\n";
		print JOB "#\$ \-o $dta_path/$outputName.o\n";
		print JOB $cmd;
		close(JOB);
		if (defined($params->{'digestion'}) and $params->{'digestion'} eq 'partial')
		{
			system(qq(qsub -cwd -cwd -pe mpi 64 -l mem_free=100G,h_vmem=100G $FileName >/dev/null 2>&1));
		}
		else
		{
			system(qq(qsub -cwd -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G $FileName >/dev/null 2>&1));
		}
		#system(qq(qsub -cwd -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G $FileName >/dev/null 2>&1));
	}
	if($GridType eq 'PBS')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#PBS -N $outputName\n";
		print JOB "#PBS -e $dta_path/$outputName.e\n"; 
		print JOB "#PBS -o $dta_path/$outputName.o"; 			
		print JOB $cmd;
		close(JOB);
		system(qq(qsub -cwd $FileName >/dev/null 2>&1));
	}
	close(JOB);
}

sub MergeFile
{
	my $nbRanges = shift;
	my $cmd = "for i in {0..$nbRanges} do; cat "
}
	
sub Check_Job_stat
{
# test whether the job is finished
	my ($jobs_prefix,$job_num,$dta_path) = @_;
	$job_info=1;
    my ($username) = getpwuid($<);
	my $command_line="";
	my $dot = ".";
	while($job_info)
	{
		if($params->{'cluster'} eq 0)
		{
			if($jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/)
			{
			
				my @outfile = glob("$dta_path\/\*.spout");
				my @dtafile = glob("$dta_path\/\*.dta");

				my $outfile_num=scalar @outfile;
				my $dtafile_num=scalar @dtafile;
				print "\r  $outfile_num files have done         ";				
				if($outfile_num == $dtafile_num)
				{
					$job_info = 0;
				}
				else
				{
					$job_info = 1;
				}
				#chomp $outfile_num;
				
				#$command_line =  "qstat -u $username | wc -l";	
				#my $job_status=qx[$command_line];
				#my @job_status_array=split(/\n/,$job_status);
				#my $job_number = $job_num - scalar (@job_status_array);
				
				#my $count = $#file_array+1;
				#my $outfile_num = $count - ($job_number - 2) * $dta_num_per_file;

				sleep(30);
			}
		}
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";
		}
		elsif($params->{'Job_Management_System'} eq 'SGE')
		{
			$command_line =  "qstat -u $username";
		}
		my $job_status=qx[$command_line];
		#print "$job_status\n";
## end of change
		my $job_status=qx{$command_line 2>&1};
		
		my @job_status_array=split(/\n/,$job_status);
	
		#Consider only the one that we submitted
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";	

			my $job_status=qx[$command_line];
			my @job_status_array=split(/\n/,$job_status);
			my $job_number = $job_num - scalar (@job_status_array);
			if(scalar (@job_status_array) == 0)
			{
				print "\r  $job_num jobs finished          ";
			}
			else
			{
				print "\r  $job_number jobs finished          ";
				sleep(5);
			}
			if(scalar(@job_status_array)>0)
			{
				$job_info=1;				
			}
			else
			{
				$job_info=0;		
			}			
		}
		elsif($params->{'Job_Management_System'} eq 'SGE')
		{
		
			@job_status_array = grep(/$jobs_prefix/,@job_status_array);
			if($job_status=~/No unfinished job found/)
			{
				$job_info=0;
				print "  \n";
			}
			elsif((scalar(@job_status_array))==0)
			{
				$job_info=0;
			}
			elsif($job_status_array[1]=~/PEND/)
			{
				print "\r cluster is busy, please be patient!          ";
				sleep(100);
			}
			elsif($jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/)
			{
				my $check_command = "ls -f $dta_path\/\*.spout \| wc -l";
				my @outfile = glob("$dta_path\/\*.spout");

				my $outfile_num=scalar @outfile;
				#chomp $outfile_num;
				
				#$command_line =  "qstat -u $username | wc -l";	
				#my $job_status=qx[$command_line];
				#my @job_status_array=split(/\n/,$job_status);
				#my $job_number = $job_num - scalar (@job_status_array);
				
				#my $count = $#file_array+1;
				#my $outfile_num = $count - ($job_number - 2) * $dta_num_per_file;
				print "\r  $outfile_num files have done         ";
				sleep(30);
				
=head				
				print "\r  $job_number jobs have finished         ";

				$dot .= ".";
				if(length($dot) == 25)
				{
					$dot = ".";
				}
				print "\r  Searching $dot          ";
=cut				
			}
			elsif($jobs_prefix eq "job_db_" || $jobs_prefix eq "sort_" || $jobs_prefix eq "merge_"  || $jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/)
			{
				if($params->{'Job_Management_System'} eq 'LSF')
				{	
					$command_line =  "bjobs -u $username";	
				}
				elsif($params->{'Job_Management_System'} eq 'SGE')
				{
					$command_line =  "qstat -u $username";			
				}			
				my $job_status=qx[$command_line];
				my @job_status_array=split(/\n/,$job_status);
				my $job_number = $job_num - scalar (@job_status_array) + 2;
				if(scalar (@job_status_array) == 0)
				{
					print "\r  $job_num jobs finished          ";
				}
				else
				{
					print "\r  $job_number jobs finished          ";
				}
				if(scalar(@job_status_array)>0)
				{
					$job_info=1;				
				}
				else
				{
					print "\r  $job_num jobs finished          ";				
					$job_info=0;		
				}
			
			}
			
		}

	}
}


sub check_input
{
	my ($raw_file,$sim_path,$parameter)=@_;

	my $err = new Spiders::Error;
	if(!defined($raw_file))
	{
		$err->no_rawfile_error();
	}
	if(!defined($parameter))
	{
		$parameter = "spSearch.params";
	}
}

########################################################################


sub database_creation
{
	#my $fasta_file = $params->{'database_name'};
	my $databasename = $params->{'database_name'};
	my $database_path = dirname($databasename); 
	my $database_basename = basename($databasename);
	my $digestion = $params->{'digestion'};
	my $num_dynamic_mod = 1;
	my $random = int(rand(100));
	my $tmp_database_path = $database_path . "/.tmp$random/";
	foreach my $key (keys %$params)
	{
		if($key=~/dynamic_/)
		{
			$num_dynamic_mod++;
		}
	}

	if ($params->{search_engine} eq 'SEQUEST')
	{
	}
	elsif(!(-e ($databasename) and $databasename=~/.mdx/))
	{
######## Needs further test
=head	
		if($param->{'inclusion_decoy'}==1)
		{
			print "  Generating reverse database\n";
			my $makedb = new Spiders::MakeDB('db_file'=>$databasename);
			$makedb->set_parameter($params);
			$makedb->make_db_with_decoy();
		}
=cut	
	
		print "  Creating database\n";
		my $mass_index = $databasename . ".mdx";

		my $peptide_index = $databasename . ".pdx";
		my $protein_index = $databasename . ".prdx";
		if(-e($mass_index))
		{
			print "  Do you want to remove the old database with same name? (y/n): ";
			chomp(my $choice = <STDIN>);
			if($choice eq "yes" || $choice eq "y")
			{
				print "  Removing old database\n";
				system(qq(rm -f $mass_index >/dev/null 2>&1));
				system(qq(rm -f $peptide_index >/dev/null 2>&1));
				system(qq(rm -f $protein_index >/dev/null 2>&1));
			}	
		}
		$params->{'database_name'} = $mass_index;
		################# create database using cluster system ########################
		my $total_protein_number=0;
		open(FASTA, $databasename) || die "can not open the database\n";
		while(<FASTA>)
		{
			$total_protein_number++ if($_=~/^\>/);
		}		
		close(FASTA);
		
		my $num_protein_per_job;
		#my $num_protein_per_job = int($total_protein_number/(200*$num_dynamic_mod))+1;
		if ($total_protein_number>2000) {
			$num_protein_per_job = int($total_protein_number/(200*$num_dynamic_mod))+1; # v12 statement (for regular database, e.g., human proteome)
		} else {
			$num_protein_per_job = int($total_protein_number/(5*$num_dynamic_mod))+1; # for small test
		}

		my $protein_num=0;
		my $k=0;
		
		open(FASTA, $databasename) || die "can not open the database\n";
		system(qq(mkdir $tmp_database_path >/dev/null 2>&1));
		while(<FASTA>)
		{
			$protein_num++ if($_=~/^\>/);
			if(($protein_num % $num_protein_per_job)==1)
			{
				if($_=~/^\>/)
				{
					$k++;
					print "\r  Generating $k temporary files";
					open(FASTATEMP,">$tmp_database_path/temp_${k}_${database_basename}");
				}
				print FASTATEMP "$_";

			}
			else
			{
				print FASTATEMP "$_";		
			}
		}
		close(FASTA);
		print "\n";
		
	########### make job file ##############	
		my $job = new Spiders::Job();
		my $abs_parameter = abs_path ($parameter);
		$job->set_library_path($library);			
	#	my $num_mass_region = 20 * $num_dynamic_mod;
		## Modified by JCho on 12/02/2014
		my $num_mass_region = 40;

		$job->make_createdb_script($tmp_database_path,$abs_parameter,$num_mass_region,$num_protein_per_job);

	######### submit job file ##############	
		my $job_num=int($protein_num/$num_protein_per_job)+1;	
		for($i=1;$i<=$job_num;$i++)
		{
			open(JOB,">$tmp_database_path/job_db_$i.sh") || die "can not open the job db files\n";
			if($params->{'cluster'} eq '1') {
			if($params->{'Job_Management_System'} eq 'LSF')
			{
				print JOB "#BSUB -P prot\n";
				print JOB "#BSUB -q normal\n";
				print JOB "#BSUB -M 20000\n";
				print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
				print JOB "#BSUB -eo $tmp_database_path/$i.e\n";
				print JOB "#BSUB -oo $tmp_database_path/$i.o\n";
				print JOB "perl $tmp_database_path/create_db.pl $tmp_database_path/temp_${i}_${database_basename}\n";
				print JOB "rm -f $tmp_database_path/temp_${i}_${database_basename}";
				close(JOB);
				#system(qq(cd $tmp_database_path && bsub <job_db_$i.sh >/dev/null 2>&1));	
				
			}
			if($params->{'Job_Management_System'} eq 'SGE')
			{
				print JOB "#!/bin/bash\n";
				print JOB "#\$ \-N job_db_$i\n";
				print JOB "#\$ \-e $tmp_database_path/$i.e\n";
				print JOB "#\$ \-o $tmp_database_path/$i.o\n";
				print JOB "perl $tmp_database_path/create_db.pl $tmp_database_path/temp_${i}_${database_basename}\n";
				print JOB "rm -rf $tmp_database_path/temp_${i}_${database_basename}\n";
				close(JOB);
				#system(qq(qsub -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G "$tmp_database_path/job_db_$i.sh" >/dev/null 2>&1));
			}
			} else {
				print JOB "perl $tmp_database_path/create_db.pl $tmp_database_path/temp_${i}_${database_basename} >$tmp_database_path/$i.o 2>$tmp_database_path/$i.e\n";
				print JOB "rm -rf $tmp_database_path/temp_${i}_${database_basename}\n";
			}
			close(JOB);
		}

		my ($MAX_PROCESSES);

		# submit jobs
		if($params->{'cluster'} eq '1') {
			if($params->{'Job_Management_System'} eq 'LSF') {
			} elsif($params->{'Job_Management_System'} eq 'SGE') {
				for($i=1;$i<=$job_num;$i++) {
					system(qq(qsub -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G "$tmp_database_path/job_db_$i.sh" >/dev/null 2>&1));
				}
			}
		
		print "  You submit $job_num jobs for creating index files\n";
		
		Check_Job_stat("job_db_",$job_num);
		} elsif($params->{'cluster'} eq '0') {
			$MAX_PROCESSES=defined($params->{'processors_used'})?$params->{'processors_used'}:4;
			my $pm = new Parallel::ForkManager($MAX_PROCESSES);
			for my $i ( 1 .. $job_num ) {
				$pm->start and next;
				system(qq(sh $tmp_database_path/job_db_$i.sh >/dev/null 2>&1));
				print "\r  $i jobs were submitted";
				$pm->finish;
			}
			$pm->wait_all_children;
		}

	###### Merge all temp_ files into big data file #######################
		
		my $j=0;
		my $l=0;
		my $prev_run = 0;
		print "\n  Sorting indexes\n";

		$job->make_partialidx_script($tmp_database_path);
		
		my $prot_job_num = int($num_mass_region/2);

#		my $inclusion_decoy = $params->{'inclusion_decoy'};

		my $tmp_mass_index = "$tmp_database_path/$database_basename" . ".mdx";
		my $tmp_peptide_index = "$tmp_database_path/$database_basename" . ".pdx";
		my $tmp_protein_index = "$tmp_database_path/$database_basename" . ".prdx";
		for($m=0;$m<$num_mass_region;$m++)
		{ 		
			
			my ($mass_hash,$peptide_hash,$protein_hash);	
			#create bash files
			my $Params = {'range'=> $m,				
							  'GridType'=> $params->{'Job_Management_System'},
							  'JobNum'=>$job_num,
							  'database_path'=> $tmp_database_path,
							  'dta_path'=> $tmp_database_path,
							  'mass_index'=>  $tmp_mass_index.".$m",
							  'peptide_index'=> $tmp_peptide_index.".$m",
							  'protein_index'=> $protein_index,
							   'databasename'=> $database_basename,
							   'num_pro_per_job'=>$num_protein_per_job,
							   'prot_job_num'=>$prot_job_num,
							   'digestion'=>$digestion,
							   };
							   

										
			my $dta_path=$tmp_database_path;
			my $FileName = "$dta_path/sort_db_".$Params->{"range"}.".sh"; #print "$FileName\n";
			my $cmd = join(" ","perl $dta_path/Create_Partial_Idx.pl -m",$Params->{"range"},"-job_num",$Params->{"JobNum"},"-dta_path",$dta_path,"-database_path",$dta_path,"-mass_index",$Params->{"mass_index"},"-peptide_index", $Params->{"peptide_index"},"-protein_index",$Params->{"protein_index"},"-databasename",$Params->{"databasename"},"-num_pro_per_job",$Params->{"num_pro_per_job"},"-prot_job_num",$Params->{"prot_job_num"},"-digestion",$Params->{"digestion"},"\n");
			
			my $outputName="sort_db_".$Params->{"range"};

			open(JOB,">$FileName") || die "can not open $FileName\n";

			if ($params->{'cluster'} eq '1') {
				if($params->{'Job_Management_System'} eq 'LSF') {
				} elsif($params->{'Job_Management_System'} eq 'SGE') {
					print JOB "#!/bin/bash\n";
					print JOB "#\$ \-S /bin/bash\n";
					print JOB "#\$ \-N $outputName\n";
					print JOB "#\$ \-e $dta_path/$outputName.e\n";
					print JOB "#\$ \-o $dta_path/$outputName.o\n";
					print JOB $cmd;
				}
			} elsif($params->{'cluster'} eq '0') {
				print JOB "#!/bin/bash\n";
				print JOB "$cmd >$dta_path/$outputName.o 2>$dta_path/$outputName.e";
			}

			close JOB;
		}

		if ($params->{'cluster'} eq '1') {

			if($params->{'Job_Management_System'} eq 'LSF') {
			} elsif($params->{'Job_Management_System'} eq 'SGE') {
				for($m=0;$m<$num_mass_region;$m++) {
					my $dta_path=$tmp_database_path;
					my $FileName = "$dta_path/sort_db_".$m.".sh";
					if (defined($params->{'digestion'}) and $params->{'digestion'} eq 'partial') {
						system(qq(qsub -cwd -cwd -pe mpi 64 -l mem_free=100G,h_vmem=100G $FileName >/dev/null 2>&1));
					} else {
						system(qq(qsub -cwd -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G $FileName >/dev/null 2>&1));
					}
				}
			}

			print "  You submit $num_mass_region jobs for sorting index files\n";
			
			Check_Job_stat("sort_",$num_mass_region);
		} elsif($params->{'cluster'} eq '0') {
			my $pm;
			if (defined($params->{'digestion'}) and $params->{'digestion'} eq 'partial') {
				my $k=int($MAX_PROCESSES/10)+1;
			} else {
				my $k=int($MAX_PROCESSES/2)+1;
			}
			$pm = new Parallel::ForkManager($k); # the sorting usually requires large memory
			for($m=0;$m<$num_mass_region;$m++) {
				my $dta_path=$tmp_database_path;
				my $FileName = "$dta_path/sort_db_".$m.".sh";
				$pm->start and next;
				system(qq(sh $FileName >/dev/null 2>&1));
				$pm->finish;
			}
			$pm->wait_all_children;
		}
		print "  $num_mass_region jobs finished";


		#Merge all the files
		print "\n  Mergering files\n";
		my $merge_job_num = $num_mass_region-1;

		if ($params->{'cluster'} eq '1') {
			my $cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_mass_index.".'$i'." >> $mass_index\n done\n";
			my $FileName = "$tmp_database_path/merge_mass_index.sh";
			LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_mass_index",$tmp_database_path);		
			$cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_peptide_index.".'$i'." >> $peptide_index\n done\n";
			$FileName = "$tmp_database_path/merge_peptide_index.sh";
			LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_peptide_index",$tmp_database_path);
		} elsif($params->{'cluster'} eq '0') {
			my $cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_mass_index.".'$i'." >> $mass_index\n done\n";	
			my $FileName = "$tmp_database_path/merge_mass_index.sh";
			open(JOB,">$FileName") || die "can not open $FileName\n";
			print JOB "$cmd >$dta_path/merge_mass_index.o 2>$dta_path/merge_mass_index.e";
			close(JOB);
			system(qq(sh $FileName >/dev/null 2>&1));

			my $cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_peptide_index.".'$i'." >> $peptide_index\n done\n";
			$FileName = "$tmp_database_path/merge_peptide_index.sh";
			open(JOB,">$FileName") || die "can not open $FileName\n";
			print JOB "$cmd  >$dta_path/merge_peptide_index.o 2>$dta_path/merge_peptide_index.e";
			close(JOB);
			system(qq(sh $FileName >/dev/null 2>&1));

			
		}

# to implement standalone version, old code deactiavted by Yuxin on 9/6/16
=head		
		for($m=0;$m<$num_mass_region;$m++)
		{ 		
			#print $m,"\n";
			my ($mass_hash,$peptide_hash,$protein_hash);	
			#create bash files
			my $parameters = {'range'=> $m,				
							  'GridType'=> $params->{'Job_Management_System'},
							  'JobNum'=>$job_num,
							  'database_path'=> $tmp_database_path,
							  'dta_path'=> $tmp_database_path,
							  'mass_index'=>  $tmp_mass_index.".$m",
							  'peptide_index'=> $tmp_peptide_index.".$m",
							  'protein_index'=> $protein_index,
							   'databasename'=> $database_basename,
							   'num_pro_per_job'=>$num_protein_per_job,
							   'prot_job_num'=>$prot_job_num,
							   'digestion'=>$digestion,
							   };
							   #'inclusion_decoy'=>$inclusion_decoy,

			Create_Sort_BashFile($parameters,$tmp_database_path)							
		}
		print "  You submit $num_mass_region jobs for sorting index files\n";
		#Waint For the jobs until they get finished
		Check_Job_stat("sort_",$num_mass_region);
		print "  $num_mass_region jobs finished";		


		#Merge all the files		
		print "\n  Mergering files\n";
		my $merge_job_num = $num_mass_region-1;
		
		my $cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_mass_index.".'$i'." >> $mass_index\n done\n";
		my $FileName = "$tmp_database_path/merge_mass_index.sh";
		LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_mass_index",$tmp_database_path);		
		$cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_peptide_index.".'$i'." >> $peptide_index\n done\n";
		$FileName = "$tmp_database_path/merge_peptide_index.sh";
		LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_peptide_index",$tmp_database_path);
=cut
		
	### create a script to sum the numbers #############	
		my @summaryfiles=();
		for($i=0;$i<$merge_job_num;$i++)
		{
			my $sumfile = "$tmp_mass_index.${i}.summary";
			push(@summaryfiles,$sumfile);
		}
	## generate a database summary file#####	
		my $outfile = $mass_index;
		$outfile =~ s/mdx/sdx/;
		summary_db(\@summaryfiles,$outfile);

	####################################################	
		print "  You submit 2 jobs for merging index files\n";
		Check_Job_stat("merge_","2");		
		print "  2 jobs finished";			
		#Clean Temporary files
		print "\n  Removing temporary files\n";
	#	system(qq(mv $tmp_database_path/$mass_index $database_path/$mass_index >/dev/null 2>&1));
	#	system(qq(mv $tmp_database_path/$peptide_index $database_path/$peptide_index >/dev/null 2>&1));
	#	system(qq(mv $tmp_database_path/$protein_index $database_path/$protein_index >/dev/null 2>&1));		
	#	system(qq(mv $tmp_database_path/$outfile $database_path/$outfile >/dev/null 2>&1));
	#	system(qq(rm -rf $tmp_database_path >/dev/null 2>&1));
		if (defined($params->{'temp_file_removal'}) and $params->{'temp_file_removal'}==1) {
			system(qq(rm -rf $tmp_database_path >/dev/null 2>&1)); # deactivated for debug
		}
				
		print "  Database creation completed\n";
	}		
}

sub summary_db
{
	my ($files,$outfile)=@_;

	my $fasta_file = $params->{'database_name'};
	my $enzyme_info = $params->{'enzyme_info'};
	my $digestion = $params->{'digestion'};	
	my $mis_cleavage = $params->{'max_mis_cleavage'};
	my $min_mass = $params->{'min_peptide_mass'};	
	my $max_mass = $params->{'max_peptide_mass'};	
	
	my %hash;
	foreach my $file(@$files)
	{
		open(FILE,$file);
		while(<FILE>)
		{
			my ($key,$value) = split(/\:/,$_);
			if(!defined($hash{$key}))
			{
				$hash{$key}=$value;
			}
			else
			{
				$hash{$key} += $value;
			}
		}
	}
	close(FILE);
	open(OUTPUT,">$outfile");
	print OUTPUT "JUMP Database index file\n";
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon,$mday,$hour,$min,$sec;	
	print OUTPUT "Database_name = ", $fasta_file,"\n";
	print OUTPUT "enzyme_info = ", $enzyme_info,"\n";
	print OUTPUT "digestion = ", $digestion,"\n";
	print OUTPUT "max_mis_cleavage = ", $mis_cleavage,"\n";
	print OUTPUT "min_peptide_mass = ", $min_mass,"\n";
	print OUTPUT "max_peptide_mass = ", $max_mass,"\n";
	
	foreach (sort keys %hash)
	{
		print OUTPUT $_," = ", $hash{$_},"\n";
	}
	close(OUTPUT);	
}

sub topMS2scans
{
        my ($msmshash)=@_;

        my (%hash,@array,$preNumber);

        foreach my $number (keys %{$msmshash})
        {
                $preNumber=$$msmshash{$number}{'survey'};#print "$preNumber\n";
                if (defined($preNumber))
                {
                #if (exists($hash{$preNumber})) { $hash{$preNumber}++;  print "mark1\n";}
                if (defined($hash{$preNumber})) { $hash{$preNumber}++; }#print "mark1\n";}
                else { $hash{$preNumber}=1; }#print "mark2\n";}
                }
        }

        for (my $i=0; $i<=10; $i++) { $array[$i]=0; }
        foreach my $number (keys %hash)
        {
                for (my $i=0; $i<$hash{$number}; $i++) { $array[$i]++; }
        }

        return \@array;
}

sub delete_run_dir
{
	my ($run)=@_;
	my $tmp = int(rand(100000));
	my $tmp = "." . $tmp;
	system(qq(mkdir $tmp;cp $run\/*.params $tmp;cp $run\/*.pl $tmp));
	system(qq(rm -rf $run > /dev/null 2>&1));
	system(qq(mkdir $run > /dev/null 2>&1));
	system(qq(mv $tmp/* $run/));
	system(qq(rm -rf $tmp > /dev/null 2>&1));
}


sub make_tags
{
	my ($dta_path,$output)=@_;
	my @tagfile = glob("$dta_path\/\*.tag");
	open(OUT,">$output");
        my $k=0;
        my $tt=scalar(@tagfile);
        foreach my $tag (@tagfile)
        {
                $k++;
                print "\r  collecting $k out of $tt tag files";

                my @t=split /\//,$tag;
                $tag=$t[$#t];

                print OUT "$tag\n";

                my $fullpath="$dta_path\/$tag";
                open(IN,$fullpath) or die "cannot open $fullpath";
                my @lines=<IN>;
                close IN;

                for (my $i=0; $i<=$#lines;$i++) { print OUT "$lines[$i]"; }
        }
        close OUT;
        print "\n";
}


sub usage {

print <<"EOF";
	
################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump search using JUMP                  ****     # 
#       ****  Version 12.1.3                          ****     #
#       ****  Copyright (C) 2012 - 2015               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################

Usage: $progname -p parameterfile rawfile.raw 
	or
       $progname -p parameterfile rawfile.mzXML
	

EOF
exit 1;
}
