#!/usr/bin/perl  -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/s
use Getopt::Long;

#use Cwd 'abs_path';
#use Storable;
#use File::Basename;
use Spiders::Params;
#use Spiders::Deisotope;
#use Spiders::Consolidation;
#use Spiders::Simulation;
#use Spiders::ProcessingRAW;
#use Spiders::ProcessingMzXML;
#use Spiders::Decharge;
#use Spiders::Sequest;
use Spiders::Tag;
#use Spiders::Error;
use Spiders::Search;
use Spiders::Dta;
#use Spiders::Dtas;

my ($help,$parameter,$sim_path);
GetOptions('-help|h'=>\$help,
		'-job_num=s'=>\$job_num,
		'-param=s'=>\$parameter,
		'-dta_path=s'=>\$dta_path,
		);
my @dtafiles = @ARGV;		
my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();	
$params->{'Mn_Mn1'} = 0.5;
$params->{'M_M1'} = 0.3;
$params->{'low_mass'} = 57;
$params->{'high_mass'} = 187;
$params->{'tag_coeff_evalue'} = 1;
$params->{'pho_neutral_loss'} = 0;

my $mass_tolerance = $params->{'peptide_tolerance'};

my $mass_tolerance_units = $params->{'peptide_tolerance_units'};
my $frag_tolerance = $params->{'frag_mass_tolerance'};
my $tag_search_method = $params->{'tag_search_method'};
my $max_number_tags_for_search = $params->{'max_number_tags_for_search'};
my $min_tag_length = $params->{'min_tag_length'};
my $databasename = $params->{'database_name'};

#my $dynamic_mass_tolerance_hash = retrieve("$dta_path/.dynamic_mass_tolerance") if($params->{'vary_tolerance'});
#my $pip = retrieve("$dta_path/.pip_hash");
#my $dtas = Spiders::Dtas->new(); # turn off dtas

$params->{'MS2_deisotope'}=0;
$params->{'ms2_consolidation'}=100;

foreach my $dta_file (@dtafiles)
{	
	next if($dta_file eq ".");
	next if($dta_file eq "..");
	next if($dta_file !~/.dta/);
	my $prec_peak_int;
	if($dta_file =~ /(([A-Za-z0-9\_\-]+)\.(\d+).(\d+).(\d+).dta)/)
	{
		my $scan = $3;
		my $order = $4;
		#$prec_peak_int = $$pip{$scan}{$order}{'pip'};	
		$prec_peak_int = 0;	
	}
#	my $decharge = new Spiders::Decharge();
#	$decharge->set_parameter($params);
#	$decharge->set_dta_path($dta_path);	
#	my $dta_hash = $decharge->create_dtahash($dta_file,\%msms_hash);
#	my $dta_file = $decharge->decharge($dta_file, $dta_hash, \%msms_hash, \%ms_hash, \@mz_array, \%realcharge);

	my $dtafile = "$dta_file";
	my $ms2_signal_noise_ratio = 0;
	my $prec_mass = 0;

#	if(1)
#	{
#		my $dta = new Spiders::Dta();
#		$dta->set_dta_file($dta_file);
#		$dta->process_dtafile();
#		$dta->parse_dtafile_name();	
#		my $charge = $dta->get_charge();
#		$prec_mass = $dta->get_prec_mz();
#		my $exp_mz_int = $dta->get_mz_int_hash();
#		undef $dta;
#		my $deisotope = new Spiders::Deisotope();
#		$deisotope->set_parameter($params);	
#		$ms2_signal_noise_ratio = $deisotope->calculate_signal_noise_ratio($prec_mass,$charge,$exp_mz_int);
#		undef $deisotope;
#	}
#	else
#	{
=head
		if($params->{'MS2_deisotope'}==1)
		{
			my $deisotope = new Spiders::Deisotope();
			$deisotope->set_parameter($params);
		
	#	my $dtafile = abs_path($dtafile);
		
			$deisotope->set_dta($dtafile);
			($ms2_signal_noise_ratio,$pho_neutral_loss) = $deisotope->MS2_deisotope();
	#		$deisotope->print_mass_error("$dtafile.mserr");
			undef $deisotope;
		}	
		my $consolidation = new Spiders::Consolidation('-dta_file'=>$dtafile,'-keepnum'=>$params->{'ms2_consolidation'});
		$consolidation->set_parameter($params);
		$consolidation->Consolidation();

		undef $consolidation;
=cut
#	}
=head	
	if($params->{'simulation'}==1)
	{	
		$sim = new Spiders::Simulation();
		$sim->set_dta_file($dtafile);
		$sim->set_param($params);
		$sim->simulation();
		undef $sim;
	}
=cut
=head
	if($params->{'search_engine'} eq "SEQUEST")
	{
		if($params->{'tag_generation'} == 1)
		{
			my $tag = new Spiders::Tag();
			$tag->set_parameter($params);
			$tag->set_dta($dtafile);
			my $msms_hash_dta=$tag->get_msms_dta($dtafile);
			$prec_mass = $tag->get_precursor_mz();
			my $cand_tag = $tag->derive_tag($dtafile,$msms_hash_dta);
			my $tag_rank_num = 0;

			my $tagfile = $searchfile = $dtafile;
			$tagfile =~ s/dta/tag/;
			$searchfile =~ s/dta/spout/;	
			
			my @searchtagarray;	
			foreach $sel_tag (@$cand_tag)
			{	
				next if (length($sel_tag->{'tag'})<$min_tag_length);	
				$tag_rank_num++;
				my $search_tag = $tag->construct_search_tag($tagfile,$sel_tag);
				$search_tag->{'tag_rank_num'}=$tag_rank_num;
				
				push (@searchtagarray,$search_tag);		
			}
		}		
		my $sequest= new Spiders::Sequest();
		$sequest->set_sequest_path("/data1/pipeline/release/version11.1.1");
		$sequest->Run_Sequest_standalone($dtafile);
		
	}
=cut
#	else # JUMP
#	{
		my @searchtagarray;
		if( defined($params->{'use_existed_tag'}) and $params->{'use_existed_tag'}==1 )
		{
			my $tag = new Spiders::Tag();
			$tag->set_parameter($params);
			$tag->set_dta($dtafile);
			my $msms_hash_dta=$tag->get_msms_dta($dtafile);
			$prec_mass = $tag->get_precursor_mz();
			#my $cand_tag = $tag->derive_tag($dtafile,$msms_hash_dta);

			my $tagfile = $searchfile = $dtafile;
			$tagfile =~ s/dta/tag/;	
			$searchfile =~ s/dta/spout/;			
			$searchtagarray_ref = $tag->construct_search_tag_from_file($tagfile);
			@searchtagarray = @$searchtagarray_ref;
			undef $tag;
		}
		else
		{
			my $tag = new Spiders::Tag();
			$tag->set_parameter($params);
			$tag->set_dta($dtafile);
			my $msms_hash_dta=$tag->get_msms_dta($dtafile);
			$prec_mass = $tag->get_precursor_mz();
			my $cand_tag = $tag->derive_tag($dtafile,$msms_hash_dta);
			my $tag_rank_num = 0;

			my $tagfile = $searchfile = $dtafile;
			$tagfile =~ s/dta/tag/;
			$searchfile =~ s/dta/spout/;	
			
	
			foreach $sel_tag (@$cand_tag)
			{	
				next if (length($sel_tag->{'tag'})<$min_tag_length);	
				$tag_rank_num++;
				my $search_tag = $tag->construct_search_tag($tagfile,$sel_tag);
				$search_tag->{'tag_rank_num'}=$tag_rank_num;
				
				push (@searchtagarray,$search_tag);		
			}
			undef $tag;			
		}
		
		my $search = new Spiders::Search();
		$search->set_parameter($params);
		$search->set_scanNum($scan);
		$search->set_dtafile($dta_file);
		
		$search->set_database("$databasename");
		$search->set_precMass($prec_mass);
		
		$search->set_mass_tolerance_units($mass_tolerance_units);

		$search->set_pho_neutral_loss($pho_neutral_loss);
		
		my $dta = new Spiders::Dta();
		$dta->set_dta_file($dta_file);
		$dta->process_dtafile();
		$dta->parse_dtafile_name();	
		my $charge = $dta->get_charge();
		my $exp_mz = $dta->get_mz_array();
		my $exp_mz_int = $dta->get_mz_int_hash();
		undef $dta;
		
		$search->set_exp_mz($exp_mz);
		$search->set_exp_mz_int($exp_mz_int);
		
		$search->set_precCharge($charge);
		
		if($params->{'vary_tolerance'})
		{
			my $dynamic_mass_tolerance = $dynamic_mass_tolerance_hash->{$dta->get_scan_num()};
			$search->set_mass_tolerance($dynamic_mass_tolerance);	
		}
		else
		{
			$search->set_mass_tolerance($mass_tolerance);
		}
		
		$search->set_frag_tolerance($frag_tolerance);	
		#my $searchmassresults = $search->SearchMass(1);
		my $searchmassresults = $search->SearchMass($params->{second_search});
		my $updatedresults;
		if($tag_search_method == 1)
		{
			foreach $search_tag (@searchtagarray)
			{
				$search->set_tag($search_tag);
				my $searchresults = $search->SearchTag($searchmassresults);
				if($searchresults)
				{
					$updatedresults = $search->mergeresults($updatedresults,$searchresults,$search_tag->{'tagSeq'});
					last;
				}		
			}	
		}	
		elsif($tag_search_method == 2)
		{
			my $number4search = $max_number_tags_for_search <= $#searchtagarray ? $max_number_tags_for_search : scalar @searchtagarray;
	#		my %matched_tag=();


		
			for(my $i=0;$i<$number4search;$i++)
			{
	############## 10/14/2013 ###################		
				if($charge == 1)
				{
					next if ($i>5);
				}
	#############################################			
	#			next if ($matched_tag{$searchtagarray[$i]->{'tagSeq'}});
				$search->set_tag($searchtagarray[$i]);
				my $searchresults = $search->SearchTag($searchmassresults);
				if($searchresults)
				{
	#				$matched_tag{$searchtagarray[$i]->{'tagSeq'}}=1;
					$updatedresults = $search->mergeresults($updatedresults,$searchresults,$searchtagarray[$i]->{'tagSeq'});
				}		
			}
	#		undef %matched_tag;
		}
			
		$search->set_tag_number(scalar (@searchtagarray));	

		
		if(scalar( keys %$updatedresults)<1)
		{
			
			my $updatedresults = $search->SearchWithoutTag($searchmassresults);
			$search->WriteResults4NoTags($searchfile,$updatedresults,$tag_rank_num,$prec_peak_int,$ms2_signal_noise_ratio,\@searchtagarray);
		}
		else
		{
			$search->WriteResults($searchfile,$updatedresults,$tag_rank_num,$prec_peak_int,$ms2_signal_noise_ratio);
		}
		undef $updatedresults;
		
		undef $search;
		print "$dta_file searched
";
		#$dtas->add_dta($dtafile);		 # turn off dtas
	#	system(qq(rm -rf $dta_file));		
	#	system(qq(rm -rf $dta_file)) if(1);
#	}

}
#$dtas->print_dtas("job_$job_num.dtas");	 # turn off dtas
