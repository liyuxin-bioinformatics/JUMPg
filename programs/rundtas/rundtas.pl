#!/usr/bin/perl -w -I /home/yli4/development/JUMPg/JUMPg_v2.3.6/programs/s 

use strict;
#use outFileParser;
use Spiders::SpoutParser;
use Parallel::ForkManager;
use Cwd;
use Spiders::Path;

if (scalar(@ARGV)!=2)
{
	die "Usage: perl rundtas rundtas.params qc_MSMS_input.txt\n";
}

# program path
#my $currectDir=getcwd;
#$runsearch_shell="$currectDir/runsearch_shell.pl";
my $runsearch_shell="/home/yli4/development/JUMPg/JUMPg_v2.3.6/programs/rundtas/runsearch_shell.pl";

#initialization
my (%parahash);
parse_params($ARGV[0],\%parahash);
parse_inputs($ARGV[1],\%parahash);
makeOutputPaths(\%parahash);

# manipulation and search for each fraction
foreach my $run (sort keys %{$parahash{'runs'}})
{
	print "\n  Running search for $parahash{'runs'}{$run}\n";

	# time:
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	print "  Search starts: ";
	printf "%4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

	# generate dta (and tag if available) files 
	print "\n  Generating .dta files:\n";
	next unless (generate_dta_files(\%parahash,$run));
	print "\n  Generating .tag files:\n";
	generate_tag_files(\%parahash,$run);

	# run search
	if ( $parahash{'search_engine'} eq 'JUMP' )
	{
		# normal run jobs
		print "\n  Running jump search\n";
		my $dta_path=$parahash{'outputPaths'}{$run};
		my @file_array = glob("$dta_path/*.dta");
		my $random = int(rand(100));
		#my $temp_file_array=runjobs(\@file_array,$dta_path,"sch_${random}");
		my ($done,$temp_file_array)=runjobs(\@file_array,$dta_path,"sch_${random}");
		#$done=1; # some scans will never generate spout if no candidate peptide exist (?)

		# re-run for 'missed' dta files (e.g., due to dead nodes)
		my $rerunN=0;
		while($done==0 and $rerunN<3)
		#while($done==0 and $rerunN<10)
		{
			$rerunN++;
			print "\n  ",scalar(@$temp_file_array)," .dta files not finished! Doing re-search (rerunN = $rerunN)\n";
			($done,my $remaining)=runjobs($temp_file_array,$dta_path,"rescue_$rerunN");
			$temp_file_array=$remaining;
		}
		if($done==0)
		{
			print "\n  Warning: ",scalar(@$temp_file_array)," .dta files not finished!\n";
		}

		# generate .tags file
		#if (!defined($parahash{'use_existed_tag'}) or $parahash{'use_existed_tag'}==0)
		if (1)
		{
			print "  Generating tags file\n";
			make_tags(\%parahash,$run);
		}

		# generate pepXML
		#if (defined($parahash{generate_pepXML}) && $parahash{generate_pepXML})
		if (1)
		{
			print "  Generating pepXML file\n";
			my $spout=Spiders::SpoutParser->new();

			my $outXML = $dta_path . ".pepXML";
			my (%seqParaHash, %outInforHash);
			$spout->readSpsearchParam(\%seqParaHash,"$dta_path/jump.params");

			$spout->parseSpout(\%seqParaHash, \%outInforHash, "$dta_path");
			$spout->printPepXML(\%seqParaHash, \%outInforHash, "$dta_path", $outXML, 10);
		}
	}
	elsif ( $parahash{'search_engine'} eq 'SEQUEST' )
	{
	}

	# rm small files
	if (defined($parahash{'temp_file_removal'}) and $parahash{'temp_file_removal'}==1)
	{
		print "  Removing small files\n";
		delete_run_dir($parahash{'outputPaths'}{$run});
	}

	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	print "\n  End: ";
	printf "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
}

=head
#run search
if ( $parahash{search_engine} eq 'sequest' )
{
	foreach my $run (keys %{$parahash{'runs'}})
	{
		print "\nRunning sequest for $run\n";
		chomp(my $pwd = `pwd`);
		$run="$pwd/$run";
		ex_sequest_jobs($run,\%parahash);

		# check if dta and out flies are equal
	        my ($ndta,$nout)=count_dta_out($run);
        	print "  There are $ndta dta files and $nout out files generated: finished \% = ";
	        if ($ndta) {printf "%.6f", $nout*100/$ndta;} else {print 'NA';} print "\n";
        	if ($nout < $ndta) { print "Search not complete!!!\n"; }
	 	else{
			# generate pepXML
			if (defined($parahash{generate_pepXML}) && $parahash{generate_pepXML})	{ generate_pepXML($run); }
		}
	}
}
else
{
}
=cut
# time:

print "Done!\n";

#------------------------------------------------------
sub delete_run_dir
{
        my ($run)=@_;
        my $tmp = int(rand(100000));
        $tmp = "." . $tmp;
        system(qq(mkdir $tmp;cp $run\/*.params $tmp;cp $run\/*.pl $tmp));
                system(qq(rm -r $run)); system(qq(mv $tmp $run));
}

sub makeOutputPaths
{
	my ($parahash)=@_;

	# list all inputs
	print "  Using the following input paths:\n";
	foreach my $run (sort values %{$$parahash{'runs'}}) { print "  $run\n"; }

	open(LOC,">\.jump_s_tmp");
	# Create the path for multiple runs
	foreach my $runSymbel (sort keys %{$$parahash{'runs'}})
	{
		my $run=$$parahash{'runs'}{$runSymbel};

		# get fraction name
		# /home/yli4/development/spectrumQC/051515/qc_10fr_test1/ad_pl07
		my @t=split "\/",$run;
		my $fr=$t[$#t];
		my $directory=getcwd;

		# mkdir: 1st layer
		system(qq(mkdir $directory/$fr >/dev/null 2>&1));
		system(qq(touch $directory/$fr/$fr\.raw >/dev/null 2>&1));

		# mkdir: 2nd layer
		my $datafile = "$directory/$fr";
		my $path = new Spiders::Path($datafile);
		my $list = $path->make_directory_list();
		my $newdir;
		if(@$list)
		{
			if (!defined($$parahash{'fast_run'}) or $$parahash{'fast_run'}==0)
			{ $newdir = $path->choose_dir_entry($list,"  Choose a .out file directory",$newdir);}
			else
			{ $newdir = $path->choose_default_entry($list,"  Choose a .out file directory",$newdir); }
		}
		else
		{
			$newdir = $fr . ".1";
			if (!defined($$parahash{'fast_run'}) or $$parahash{'fast_run'}==0)
			{ $newdir = $path->ask("  Choose a .out file directory",$newdir);}
		}
		print "  Using: $newdir\n";
		$path->add_subdir($newdir);

		# store $dir as .outputPath. in %paramhash
		my $dir =  $path->basedir() . "/$newdir";
		$$parahash{'outputPaths'}{$runSymbel}=$dir;
		print LOC "$dir\n";

		# cp jump.params to $dir
		system(qq(cp -rf $ARGV[0] "$dir/jump.params" >/dev/null 2>&1));
		$$parahash{$runSymbel}{search_params}="$dir/jump.params";
	}
	close LOC;
}

sub count_dta_out
{
        my ($run)=@_;
        my ($a,$b,%ha,%hb);

        opendir(DIR,$run);
        while(my $file=readdir(DIR))
        {
                if ($file =~ /\.(\d+)\.(\d+)\.(\d)\.dta$/) { $ha{"$1\.$2\.$3"}=1; }
                elsif ($file =~ /\.(\d+)\.(\d+)\.(\d)\.out$/) { $hb{"$1\.$2\.$3"}=1; }
        }
        $a=scalar(keys %ha);
        $b=scalar(keys %hb);
        closedir DIR;

        return ($a,$b);
}



sub generate_pepXML
{
	my ($dta_path)=@_;
	my $pOut=outFileParser->new();
	my $outXML = $dta_path;
	my $run=$dta_path;
	my $topHit=30;
	my (%seqParaHash, %outInforHash);

	$pOut->readSeqParam(\%seqParaHash,$run);
	$pOut->readOutFiles(\%seqParaHash, \%outInforHash, $run);
	$pOut->printPepXML(\%seqParaHash, \%outInforHash, $run,$outXML,$topHit);
	print "\n  PepXML transformation Done!\n";
}

sub make_tags
{
	my ($parahash,$run)=@_;

	my $dta_path=$$parahash{'outputPaths'}{$run};

	my $output=$$parahash{'outputPaths'}{$run};
	$output.='.tags';

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

sub generate_tag_files
{
	my ($parahash,$run)=@_;

	# locate tags file
	my $tags=$$parahash{'runs'}{$run};
	if (1)
	{
		my @t=split /\//,$tags;
		$tags.="\/$t[$#t]\.filtered\.tags";
		if (!-e $tags) 
		{ 
			print "WARNING: $tags not found for run $$parahash{'runs'}{$run}! De novo tagging is now activated\n"; 
			$$parahash{$run}{'use_existed_tag'}=0;
			return 0;
		}
		$$parahash{$run}{'use_existed_tag'}=1;
		modify_jump_params($parahash,$run);
		my $output_tags=$$parahash{'outputPaths'}{$run};
		$output_tags.='.tags';
		#system(qq(ln -s $tags $output_tags));
	}

	open(IN,$tags)  or die "Cannot open $tags!!!\n";
	my $k=0;
	my $tagfile='';
	my @lines;
	while(<IN>)
	{
		chomp;
		if (/\.tag/)
		{
			if ($tagfile ne '')
			{
				# t20.10000.1.2.tag
				my ($runname,$scan,$ppi,$z,$suff)=split /\./,$tagfile;
				if ( $scan>=$$parahash{first_scan_extraction} and $scan<=$$parahash{last_scan_extraction} )
				{
					open(OUT,">$$parahash{'outputPaths'}{$run}\/$tagfile");
					for (my $i=0; $i<scalar(@lines); $i++) { print OUT "$lines[$i]\n";}
					close OUT;
					$k++; print "  $k tag files generated\r";
				}
			}

			undef @lines;
			$tagfile=$_;
		}
		else
		{
			push @lines, $_;
		}
	}
	close IN;
	return 1;
}

sub generate_dta_files
{
	my ($parahash,$run)=@_;

	# locate dtas file
	my $dtas=$$parahash{'runs'}{$run};
	if (1)
	{
		my @t=split /\//,$dtas;
		$dtas.="\/$t[$#t]\.filtered\.dtas";
		if (!-e $dtas) 
		{ 
			print "WARNING: $dtas not found for run $$parahash{'runs'}{$run}! This run is kipped for search.\n"; 
			return 0;
		}
		my $output_dtas=$$parahash{'outputPaths'}{$run};
		$output_dtas.='.dtas';
		system(qq(ln -s $dtas $output_dtas));
	}

	open(IN,$dtas)  or die "Cannot open $dtas!!!\n";
	my $k=0;
	while(<IN>)
	{
		my @t=split(/\s/,$_);
		my $line;
		$line=<IN>; my @mz=split(/\s/,$line);
		$line=<IN>; my @int=split(/\s/,$line);

		# t20.10000.1.2.dta
		$t[0] =~ /^(.*?)\.(\d+)\.(\d+)\.\d\.dta$/;
		my ($runname,$scan,$ppi)=($1,$2,$3);
		#my $run=$1; my $ppi=$2; 

		# scan range control
		next if ( $scan<$$parahash{first_scan_extraction} or $scan>$$parahash{last_scan_extraction} );

		# ppi control
		if ( defined($$parahash{first_ppi_only}) and $$parahash{first_ppi_only}==1
		and $ppi>1 ) { next; }

		#open(OUT,">$run\/$t[0]");
		open(OUT,">$$parahash{'outputPaths'}{$run}\/$t[0]");
		print OUT "$t[1] $t[2]\n";
		for (my $i=0; $i<scalar(@int); $i++) { print OUT "$mz[$i] $int[$i]\n"; }
		close OUT;
		$k++; print "  $k dta files generated\r";
	}
	close IN;

	return 1;
}

sub modify_seqeust_params
{
	my ($parahash)=@_;

	open(IN,$$parahash{search_params}) or die "Cannot open $$parahash{search_params}!!!\n";
	my @lines=<IN>;
	close IN;

	open(OUT,">$$parahash{search_params}");
	for (my $i=0; $i<=$#lines;$i++)
	{
		my @t=split(' = ',$lines[$i]);
		if ($t[0] eq 'first_database_name')
		{
			print OUT "first_database_name = $$parahash{database}\n";
		}
		else {print OUT $lines[$i];}
	}
	close OUT;
}

sub modify_jump_params
{
	my ($parahash,$run)=@_;

	open(IN,$$parahash{$run}{search_params}) or die "Cannot open $$parahash{$run}{search_params}!!!\n";

	my @lines=<IN>;
	close IN;

	open(OUT,">$$parahash{$run}{search_params}");
	for (my $i=0; $i<=$#lines;$i++)
	{
		my @t=split(' = ',$lines[$i]);
		if ($t[0] eq 'tag_generation')
		{
			print OUT $lines[$i];
			if ( defined($$parahash{$run}{use_existed_tag}) and $$parahash{$run}{use_existed_tag}==1 )
			{
				print OUT "use_existed_tag = $$parahash{$run}{use_existed_tag}\n";
			}
		}
		else {print OUT $lines[$i];}
	}
	close OUT;
}

sub parse_params
{
	my ($par,$parahash)=@_;
	open(IN,$par);
	my $k=0;
	while (<IN>)
	{
		s/^\s+//; # rm space at front
		next if (/^#/); # skip comment lines
		chomp;

		if (/ = /) # normal parameters
		{

			s/\s*([;\#].*)$//; # delete comments

			my ($key,$value) = split(' = ',$_);
			$$parahash{$key}=$value;
		}
=head
		elsif (/\//) # potential input paths
		{
			s/\s+$//; # rm space at end
			if (-e $_)
			{
				$k++;
				$$parahash{'runs'}{"run$k"}=$_;
			}
			else { die "Folder not found in $_!!!\nPlease specify correct paths for input folders.\n"; }
		}
=cut
	}
	close IN;
}

sub parse_inputs
{
        my ($inputlist,$parahash)=@_;
        open(IN,$inputlist) or die "cannot open $inputlist";
        my $k=0;
        while (<IN>)
        {
                s/^\s+//; # rm space at front
                next if (/^#/); # skip comment lines
                s/\s+$//; # rm space at end
                if (-e $_)
                {
                        $k++;
                        $$parahash{'runs'}{"run$k"}=$_;
                }
                else { die "Folder not found in $_!!!\nPlease specify correct paths for input folders.\n"; }
        }
        close IN;
}

sub MS1_shift
{
	my ($dta,$out)=@_;

	open(IN,$dta);
	my @inarray = <IN>;
	close IN;

	open (OUT, ">$out");
	my $tmp=$inarray[0];
	my @a=split(/ /,$tmp);
	$a[0]+=100;
	$tmp="$a[0] $a[1]";
	$inarray[0]=$tmp;
	print OUT @inarray;
	close OUT;
}

sub runjobs
{
        my ($file_array,$dta_path,$job_name) = @_;
	#system("cp /home/yli4/development/rundtas/JUMPversion/runsearch_shell.pl $dta_path");
	#system("cp /home/yli4/bin/runsearch_shell.pl $dta_path");
	system("cp $runsearch_shell $dta_path");
        #my $curr_dir = getcwd;
	chomp(my $curr_dir = `pwd`);
        #my $MAX_PROCESSES = 10;
        my $MAX_PROCESSES = $parahash{'processors_used'};

	my ($job_num,$dta_num_per_file);
	#if ( defined($parahash{use_existed_tag}) and $parahash{use_existed_tag}==1 )
	#{
	#	$job_num = $parahash{job_num};
	#	$dta_num_per_file = int($#$file_array / $job_num) + 1;
	#}
	#else
	#{
	        $dta_num_per_file = 10;
        	$job_num = int($#$file_array / $dta_num_per_file) + 1;
	        ## Set the maximum number of jobs to 4000
        	if ($job_num > 4000) 
		{
	                $job_num = 4000;
                	$dta_num_per_file = int($#$file_array / $job_num) + 1;
        	}
	#}

	for(my $i = 0; $i < $job_num; $i++)
        {
                if (($i * $dta_num_per_file) > $#$file_array) {
                        $job_num = $i;
                        last;
                }
                open(JOB,">$dta_path/${job_name}_${i}.sh") || die "can not open the job files\n";
                my $dta_file_temp="";
                my @dta_file_arrays=();
                my $multiple_jobs_num  = 0;
                for(my $j = 0; $j < $dta_num_per_file; $j++)
                {
                        if(($i*$dta_num_per_file+$j)<=$#$file_array)
                        {
                                $dta_file_temp .= " $$file_array[$i*$dta_num_per_file+$j]";
                                push (@dta_file_arrays, $$file_array[$i*$dta_num_per_file+$j])
                        }
                }
                if($parahash{'Job_Management_System'} eq 'LSF')
                {
                }
                elsif($parahash{'Job_Management_System'} eq 'SGE')
                {
                        print JOB "#!/bin/bash\n";
                        print JOB "#\$ -N ${job_name}_${i}\n";
                        print JOB "#\$ -e $dta_path/${job_name}_${i}.e\n";
                        print JOB "#\$ -o $dta_path/${job_name}_${i}.o\n";
                        foreach (@dta_file_arrays)
                        {
				my @t=split /\//,$_;
                                #print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $dta_path/$parahash{search_params} -dta_path $dta_path $t[$#t]\n";
                                print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $dta_path/jump.params -dta_path $dta_path $t[$#t]\n";
                        }
                }

                close(JOB);
        }
        my $job_list;
        if($parahash{'cluster'} eq '1')
        {
                if($parahash{'Job_Management_System'} eq 'LSF')
                {
                }
                elsif($parahash{'Job_Management_System'} eq 'SGE')
                {

                        for(my $i=0;$i<$job_num;$i++)
                        {
                                my $job_name = "${job_name}_${i}.sh";
                                my $command_line = qq(cd $dta_path && qsub -cwd $job_name);
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
                Check_Job_stat_jump("${job_name}_",$job_num,$dta_path);
        }
	elsif($parahash{'cluster'} eq '0')
	{
		my $pm = new Parallel::ForkManager($MAX_PROCESSES);
		#for my $i ( 0 .. $MAX_PROCESSES )
		for my $i ( 0 .. ($job_num-1) )
		{
			$pm->start and next;
			my $job_name = "${job_name}_${i}.sh";
			system("cd $dta_path && sh $job_name >/dev/null 2>&1");
			$pm->finish; # Terminates the child process
			print "\r  $i jobs were submitted";
			
		}
		#Check_Job_stat_jump("${job_name}_",$job_num,$dta_path);
		$pm->wait_all_children;
	}
	#print "Check_Job_stat_jump passed\n";
	#### checking whether all dta files has been searched
	my $done=1;
	my $temp_file_array;
        for(my $k=0;$k<=$#$file_array;$k++)
        {
                my $data_file = $file_array->[$k];
                my $out_file = $data_file;
                $out_file =~ s/\.dta$/\.spout/;
                if(!-e $out_file)
                {
                        push (@{$temp_file_array},$data_file);
			$done=0;
                }
        }
        return ($done,$temp_file_array);

}

sub Check_Job_stat_jump
{
	# test whether the job is finished for JUMP search
	my ($jobs_prefix,$job_num,$dta_path) = @_;
        my $job_info=1;
    my ($username) = getpwuid($<);
        my $command_line="";
        my $dot = ".";
        while($job_info)
        {
                if($parahash{'Job_Management_System'} eq 'LSF')
                {
                        $command_line =  "bjobs -u $username";
                }
                elsif($parahash{'Job_Management_System'} eq 'SGE')
                {
                        $command_line =  "qstat -u $username";
                }
                my $job_status=qx[$command_line];
                $job_status=qx{$command_line 2>&1};

                my @job_status_array=split(/\n/,$job_status);

                #Consider only the one that we submitted
                if($parahash{'Job_Management_System'} eq 'LSF')
                {
                }
                elsif($parahash{'Job_Management_System'} eq 'SGE')
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
                                #my $check_command = "ls -f $dta_path\/\*.spout \| wc -l";
                                my @outfile = glob("$dta_path\/\*.spout");

                                my $outfile_num=scalar @outfile;
                                print "\r  $outfile_num files have done         ";


				#print "\r  ",scalar @job_status_array," jobs still running";

                                #sleep(30);
                                #sleep(300);
                                sleep(60);
                        }
                        elsif($jobs_prefix eq "job_db_" || $jobs_prefix eq "sort_" || $jobs_prefix eq "merge_"  || $jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/)
                        {
                                if($parahash{'Job_Management_System'} eq 'LSF')
                                {
                                        $command_line =  "bjobs -u $username";
                                }
                                elsif($parahash{'Job_Management_System'} eq 'SGE')
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

#}

sub ex_sequest_jobs
{
#######################################################
#my $sequest28single='/usr/local/bin/sequest28single';
#my $sequest28single='/state/partition1/spiders/bin/sequest28single';
#my $sequest28single='/home/xwang4/scripts/sequest28single';
my $sequest28single='/data/bin/sequest28single';
#######################################################

        my ($run,$parahash)=@_;

        my $job_num=200;
	if (defined($$parahash{job_num})) { $job_num=$$parahash{job_num}; }
        my @array = glob("$run/*.dta");
	my $dta_num_per_file = int($#array / $job_num) + 1;

        #for(my $i=0;$i<=int(scalar(@array)/$job_num);$i++)
	print "Job_Management_System: ",$$parahash{Job_Management_System},"\n";
        for(my $i=0;$i<$job_num;$i++)
        {
                open(JOB,">$run/job_$i.sh") || die "can not open the job files\n";
=head
                my $dta_file_temp="";
                for(my $j=0;$j<$dta_num_per_file;$j++)
                {
                        if(($i*$dta_num_per_file+$j)<=$#array)
                         {
                               $dta_file_temp .= " $array[$i*$dta_num_per_file+$j]";
                          }
                }
=cut
		if ($$parahash{Job_Management_System} eq 'LSF')
		{
                        print JOB "#BSUB -P prot\n";
                        print JOB "#BSUB -q normal\n";
        #               print JOB "#BSUB -M 2000\n";
        #               print JOB "#BSUB -R \"rusage[mem=20000]\"\n";
                        print JOB "#BSUB -eo $run/$i.e\n";
                        print JOB "#BSUB -oo $run/$i.o\n";
 		}
		elsif ($$parahash{Job_Management_System} eq 'SGE')
		{
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N job_$i\n";
			print JOB "#\$ -e $run/$i.e\n";
			print JOB "#\$ -o $run/$i.o\n";
		}
                for(my $j=0;$j<$dta_num_per_file;$j++)
                {
                        if(($i*$dta_num_per_file+$j)<=$#array)
                         {
				print JOB "$sequest28single $array[$i*$dta_num_per_file+$j]\n";
				#print JOB "$sequest28single $array[$i*$dta_num_per_file+$j]; rm $array[$i*$dta_num_per_file+$j]\n";
                                #$dta_file_temp .= " $array[$i*$dta_num_per_file+$j]";
                          }
                }
                #print JOB "$sequest28single $dta_file_temp\n";
                #print JOB "rm $run/job_$i.sh\n";
		close JOB;
        }

	my $job_list;
	if ($$parahash{Job_Management_System} eq 'LSF')
	{
	}
	elsif ($$parahash{Job_Management_System} eq 'SGE')
	{
		#	print "aaaaa\n";
	        for(my $i=0;$i<$job_num;$i++)
        	{
			#my $command_line = "cd $run && qsub -cwd job_" . "$i.sh ";
			my $command_line = "cd $run && qsub -cwd job_$i.sh ";
			my $job=qx[$command_line];
        	        chomp $job;
			my $job_id=0;
	                if($job=~/Job \<(\d*)\> is/)
        	        {
                	        $job_id=$1;
	                }
        	        $job_list->{$job_id}=1;
                	my $count = $i+1;
	                print "\r  $count jobs were submitted";

		}
	}
	print "\n  You submitted $job_num jobs for database search\n";
	Check_Job_stat("job_",$job_num);
	print "\n  Searching finished\n\n";
=head
        my $MAX_PROCESSES = 10;
        my $pm = new Parallel::ForkManager($MAX_PROCESSES);

                for my $i ( 0 .. int(scalar(@array)/$job_num) )
                {
                        $pm->start and next;
                        #system("cd $run && sh job_$i.sh >/dev/null 2>&1");
                        system("cd $run && sh job_$i.sh >/dev/null 2>&1");
                        $pm->finish; # Terminates the child process
                }
                $pm->wait_all_children;
=cut
}

sub Check_Job_stat
{
# test whether the job is finished
	my ($jobs_prefix,$job_num) = @_;
	my $job_info=1;
    my ($username) = getpwuid($<);
	my $command_line="";
	my $dot = ".";
	while($job_info)
	{
		sleep(10);

		#if($params->{'Job_Management_System'} eq 'LSF')
		#{
		#	$command_line =  "bjobs -u $username";
		#}
		#elsif($params->{'Job_Management_System'} eq 'SGE')
		#{
			$command_line =  "qstat -u $username";
		#}
		my $job_status=qx[$command_line];
		#print "$job_status\n";
## end of change
		#my $job_status=qx{$command_line 2>&1};
		
		my @job_status_array=split(/\n/,$job_status);
	
		#Consider only the one that we submitted
=head
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";	

			my $job_status=qx[$command_line];
			my @job_status_array=split(/\n/,$job_status);
			my $job_number = $job_num - scalar (@job_status_array) + 1;
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
=cut		
			@job_status_array = grep(/$jobs_prefix/,@job_status_array); 
			if($jobs_prefix eq "job_")
			{
		#		my $check_command = "ls -f $dta_path\/\*.spout \| wc -l";
		#		my @outfile = glob("$dta_path\/\*.spout");

		#		my $outfile_num=scalar @outfile;
				#chomp $outfile_num;
=head				
				$command_line =  "qstat -u $username | wc -l";	
				my $job_status=qx[$command_line];
				my @job_status_array=split(/\n/,$job_status);
				my $job_number = $job_num - scalar (@job_status_array) + 2;
				
				my $count = $#file_array+1;
				my $outfile_num = $count - ($job_number - 2) * $dta_num_per_file;
				print "\r  Searching $outfile_num out of $count dta files          ";
=cut
				$dot .= ".";
				if(length($dot) == 25)
				{
					$dot = ".";
				}
				print "\r  Searching $dot          ";
			}
			elsif($jobs_prefix eq "job_db_" || $jobs_prefix eq "sort_" || $jobs_prefix eq "merge_")
			{
=head
				if($params->{'Job_Management_System'} eq 'LSF')
				{	
					$command_line =  "bjobs -u $username";	
				}
				elsif($params->{'Job_Management_System'} eq 'SGE')
				{
=cut
					$command_line =  "qstat -u $username";			
#				}			
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
					$job_info=0;		
				}
			}		
#		}
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
			sleep(5);
		}
	}
}
