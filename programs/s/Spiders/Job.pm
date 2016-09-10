#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Job
######### Job #################################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
package Spiders::Job;


use strict;
use warnings;
use File::Basename;
use Storable;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.03;
######### Version 12.1.0 #######################
# change the code to allow multiple tags for searching 
######## Version 12.1.0 (2/1/2013) ########################
# add two functions: make_createdb_script and make_partialidx_script
######## Version 12.1.0 (2/5/2013 #########################
# allow specify min_tag_length for search
############################################ 
##### Version 12.1.0 (2/20/2013) ###########
# added the simulation function ###########


@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _dta_path => undef,
        _sim_path  => undef,
    };
    bless $self, $class;
    return $self;
}

sub set_dta_path
{
	my ($self,$dta_path)=@_;
	$self->{_dta_path}=$dta_path;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{_dta_path};
}

sub set_pip
{
	my ($self,$pip)=@_;
	$self->{_pip}=$pip;	
}

sub get_pip
{
	my $self=shift;
	return $self->{_pip};
}

sub set_library_path
{
	my ($self,$lib)=@_;
	$self->{_lib}=$lib;	
}

sub get_library_path
{
	my ($self)=@_;
	return $self->{_lib};	
}

sub create_script
{
	my ($self,$search_flag)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $pip = $self->get_pip();
	if(!defined($search_flag))
	{
		$search_flag = 0;
	}
	
	store(\%$pip,"$dir/.pip_hash");
	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;	
#!/usr/bin/perl  -I $lib
use Getopt::Long;

use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::Params;
use Spiders::Deisotope;
use Spiders::Consolidation;
use Spiders::Simulation;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::Decharge;
use Spiders::Sequest;
use Spiders::Tag;
use Spiders::Error;
use Spiders::Search;
use Spiders::Dta;
use Spiders::Dtas;
use Storable;

my (\$help,\$parameter,\$sim_path);
GetOptions('-help|h'=>\\\$help,
		'-job_num=s'=>\\\$job_num,
		'-param=s'=>\\\$parameter,
		'-dta_path=s'=>\\\$dta_path,
		);
my \@dtafiles = \@ARGV;		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	
\$params->{'Mn_Mn1'} = 0.5;
\$params->{'M_M1'} = 0.3;
\$params->{'low_mass'} = 57;
\$params->{'high_mass'} = 187;
\$params->{'tag_coeff_evalue'} = 1;
\$params->{'pho_neutral_loss'} = 0;

my \$mass_tolerance = \$params->{'peptide_tolerance'};

my \$mass_tolerance_units = \$params->{'peptide_tolerance_units'};
my \$frag_tolerance = \$params->{'frag_mass_tolerance'};
my \$tag_search_method = \$params->{'tag_search_method'};
my \$max_number_tags_for_search = \$params->{'max_number_tags_for_search'};
my \$min_tag_length = \$params->{'min_tag_length'};
my \$databasename = \$params->{'database_name'};

my \$dynamic_mass_tolerance_hash = retrieve("\$dta_path\/\.dynamic_mass_tolerance") if(\$params->{'vary_tolerance'});
my \$pip = retrieve("\$dta_path\/\.pip_hash");
my \$dtas = Spiders::Dtas->new();

foreach my \$dta_file (\@dtafiles)
{	
	next if(\$dta_file eq "\.");
	next if(\$dta_file eq "\.\.");
	next if(\$dta_file !~/\.dta/);
	my \$prec_peak_int;
	if(\$dta_file =~ /(([A-Za-z0-9\\_\\-]+)\\.(\\d+)\.(\\d+)\.(\\d+)\.dta)/)
	{
		my \$scan = \$3;
		my \$order = \$4;
		\$prec_peak_int = \$\$pip{\$scan}{\$order}{'pip'};	
	}
#	my \$decharge = new Spiders::Decharge();
#	\$decharge->set_parameter(\$params);
#	\$decharge->set_dta_path(\$dta_path);	
#	my \$dta_hash = \$decharge->create_dtahash(\$dta_file,\\\%msms_hash);
#	my \$dta_file = \$decharge->decharge(\$dta_file, \$dta_hash, \\\%msms_hash, \\\%ms_hash, \\\@mz_array, \\\%realcharge);

	my \$dtafile = "\$dta_file";
	my \$ms2_signal_noise_ratio = 0;
	my \$prec_mass = 0;

#	if($search_flag)
#	{
#		my \$dta = new Spiders::Dta();
#		\$dta->set_dta_file(\$dta_file);
#		\$dta->process_dtafile();
#		\$dta->parse_dtafile_name();	
#		my \$charge = \$dta->get_charge();
#		\$prec_mass = \$dta->get_prec_mz();
#		my \$exp_mz_int = \$dta->get_mz_int_hash();
#		undef \$dta;
#		my \$deisotope = new Spiders::Deisotope();
#		\$deisotope->set_parameter(\$params);	
#		\$ms2_signal_noise_ratio = \$deisotope->calculate_signal_noise_ratio(\$prec_mass,\$charge,\$exp_mz_int);
#		undef \$deisotope;
#	}
#	else
#	{
		my \$consolidation = new Spiders::Consolidation('-dta_file'=>\$dtafile,'-keepnum'=>\$params->{'ms2_consolidation'});
		if(\$params->{'TMT_Cterm_ion_removal'})
		{
			\$minPeak=\$params->{'TMT_Cterm_ion_removal'};
			\$dis="0.99703,1.00335";
			\$peakTol = 0.01;
			\$consolidation->remove_TMTc_peaks(\$dtafile,\$minPeak,\$dis,\$peakTol);
		}
		if(\$params->{'MS2_deisotope'}==1)
		{
			my \$deisotope = new Spiders::Deisotope();
			\$deisotope->set_parameter(\$params);
		
	#	my \$dtafile = abs_path(\$dtafile);
		
			\$deisotope->set_dta(\$dtafile);
			(\$ms2_signal_noise_ratio,\$pho_neutral_loss) = \$deisotope->MS2_deisotope();
	#		\$deisotope->print_mass_error("\$dtafile.mserr");
			undef \$deisotope;
		}	

		\$consolidation->set_parameter(\$params);
		\$consolidation->Consolidation();

		undef \$consolidation;
#	}
	
	if(\$params->{'simulation'}==1)
	{	
		\$sim = new Spiders::Simulation();
		\$sim->set_dta_file(\$dtafile);
		\$sim->set_param(\$params);
		\$sim->simulation();
		undef \$sim;
	}
	if(\$params->{'search_engine'} eq "SEQUEST")
	{
		if(\$params->{'tag_generation'} == 1)
		{
			my \$tag = new Spiders::Tag();
			\$tag->set_parameter(\$params);
			\$tag->set_dta(\$dtafile);
			my \$msms_hash_dta=\$tag->get_msms_dta(\$dtafile);
			\$prec_mass = \$tag->get_precursor_mz();
			my \$cand_tag = \$tag->derive_tag(\$dtafile,\$msms_hash_dta);
			my \$tag_rank_num = 0;

			my \$tagfile = \$searchfile = \$dtafile;
			\$tagfile =~ s\/dta\/tag\/;
			\$searchfile =~ s\/dta\/spout\/;	
			
			my \@searchtagarray;	
			foreach \$sel_tag (\@\$cand_tag)
			{	
				next if (length(\$sel_tag->{'tag'})<\$min_tag_length);	
				\$tag_rank_num++;
				my \$search_tag = \$tag->construct_search_tag(\$tagfile,\$sel_tag);
				\$search_tag->{'tag_rank_num'}=\$tag_rank_num;
				
				push (\@searchtagarray,\$search_tag);		
			}
		}		
		my \$sequest= new Spiders::Sequest();
		\$sequest->set_sequest_path("/data1/pipeline/release/version12.1.0");
		\$sequest->Run_Sequest_standalone(\$dtafile);
		
	}
	else
	{
		my \@searchtagarray;
#		if($search_flag)
#		{
#			my \$tag = new Spiders::Tag();
#			my \$tagfile = \$searchfile = \$dtafile;
#			\$tagfile =~ s\/dta\/tag\/;	
#			\$searchfile =~ s\/dta\/spout\/;			
#			\$searchtagarray_ref = \$tag->construct_search_tag_from_file(\$tagfile);
#			\@searchtagarray = \@\$searchtagarray_ref;
#			undef \$tag;
#		}
#		else
#		{
			my \$tag = new Spiders::Tag();
			\$tag->set_parameter(\$params);
			\$tag->set_dta(\$dtafile);
			my \$msms_hash_dta=\$tag->get_msms_dta(\$dtafile);
			\$prec_mass = \$tag->get_precursor_mz();
			my \$cand_tag = \$tag->derive_tag(\$dtafile,\$msms_hash_dta);
			my \$tag_rank_num = 0;

			my \$tagfile = \$searchfile = \$dtafile;
			\$tagfile =~ s\/dta\/tag\/;
			\$searchfile =~ s\/dta\/spout\/;	
			
	
			foreach \$sel_tag (\@\$cand_tag)
			{	
				next if (length(\$sel_tag->{'tag'})<\$min_tag_length);	
				\$tag_rank_num++;
				my \$search_tag = \$tag->construct_search_tag(\$tagfile,\$sel_tag);
				\$search_tag->{'tag_rank_num'}=\$tag_rank_num;
				
				push (\@searchtagarray,\$search_tag);		
			}
			undef \$tag;			
#		}
		
		my \$search = new Spiders::Search();
		\$search->set_parameter(\$params);
		\$search->set_scanNum(\$scan);
		\$search->set_dtafile(\$dta_file);
		
		\$search->set_database("\$databasename");
		\$search->set_precMass(\$prec_mass);
		
		\$search->set_mass_tolerance_units(\$mass_tolerance_units);

		\$search->set_pho_neutral_loss(\$pho_neutral_loss);
		
		my \$dta = new Spiders::Dta();
		\$dta->set_dta_file(\$dta_file);
		\$dta->process_dtafile();
		\$dta->parse_dtafile_name();	
		my \$charge = \$dta->get_charge();
		my \$exp_mz = \$dta->get_mz_array();
		my \$exp_mz_int = \$dta->get_mz_int_hash();
		undef \$dta;
		
		\$search->set_exp_mz(\$exp_mz);
		\$search->set_exp_mz_int(\$exp_mz_int);
		
		\$search->set_precCharge(\$charge);
		
		if(\$params->{'vary_tolerance'})
		{
			my \$dynamic_mass_tolerance = \$dynamic_mass_tolerance_hash->{\$dta->get_scan_num()};
			\$search->set_mass_tolerance(\$dynamic_mass_tolerance);	
		}
		else
		{
			\$search->set_mass_tolerance(\$mass_tolerance);
		}
		
		\$search->set_frag_tolerance(\$frag_tolerance);	
		my \$searchmassresults = \$search->SearchMass($search_flag);
		my \$updatedresults;
		if(\$tag_search_method == 1)
		{
			foreach \$search_tag (\@searchtagarray)
			{
				\$search->set_tag(\$search_tag);
				my \$searchresults = \$search->SearchTag(\$searchmassresults);
				if(\$searchresults)
				{
					\$updatedresults = \$search->mergeresults(\$updatedresults,\$searchresults,\$search_tag->{'tagSeq'});
					last;
				}		
			}	
		}	
		elsif(\$tag_search_method == 2)
		{
			my \$number4search = \$max_number_tags_for_search <= \$\#searchtagarray \? \$max_number_tags_for_search \: scalar \@searchtagarray;
	#		my \%matched_tag=();


		
			for(my \$i=0;\$i<\$number4search;\$i++)
			{
	############## 10/14/2013 ###################		
				if(\$charge == 1)
				{
					next if (\$i>5);
				}
	#############################################			
	#			next if (\$matched_tag{\$searchtagarray[\$i]->{'tagSeq'}});
				\$search->set_tag(\$searchtagarray[\$i]);
				my \$searchresults = \$search->SearchTag(\$searchmassresults);
				if(\$searchresults)
				{
	#				\$matched_tag{\$searchtagarray[\$i]->{'tagSeq'}}=1;
					\$updatedresults = \$search->mergeresults(\$updatedresults,\$searchresults,\$searchtagarray[\$i]->{'tagSeq'});
				}		
			}
	#		undef \%matched_tag;
		}
			
		\$search->set_tag_number(scalar (\@searchtagarray));	

		
		if(scalar( keys \%\$updatedresults)<1)
		{
			
			my \$updatedresults = \$search->SearchWithoutTag(\$searchmassresults);
			\$search->WriteResults4NoTags(\$searchfile,\$updatedresults,\$tag_rank_num,\$prec_peak_int,\$ms2_signal_noise_ratio,\\\@searchtagarray);
		}
		else
		{
			\$search->WriteResults(\$searchfile,\$updatedresults,\$tag_rank_num,\$prec_peak_int,\$ms2_signal_noise_ratio);
		}
		undef \$updatedresults;
		
		undef \$search;
		print "\$dta_file searched\n";
		\$dtas->add_dta(\$dtafile);		
	#	system(qq(rm -rf \$dta_file));		
	#	system(qq(rm -rf \$dta_file)) if($search_flag);
	}

}
\$dtas->print_dtas("job_\$job_num.dtas");	
EOF
	close(RUNSHELL);
}

sub make_createdb_script
{

	my ($self,$dta_path,$parameter,$num_mass_region,$num_protein_per_job)=@_;
	my $lib = $self->get_library_path();	
	open(DB,">$dta_path/create_db.pl") || die "can not open the create_db.pl file:\n";

	print DB "#!/usr/bin/perl -I $lib\n";
	print DB "my \$fasta_path = \$ARGV[0];\n";
	print DB "use Spiders::BuildIndex;\n\n";
	print DB "use Spiders::Params;\n";
	print DB "use Storable;\n";

	
	print DB "my \$p = Spiders::Params\-\>new('-path'=>\"$parameter\");\n";
	print DB "my \$params=\$p->parse_param();\n";			
	print DB "my \$index = new Spiders\:\:BuildIndex();\n";
	print DB "\$index->set_dta_path(\"$dta_path\");\n";
	print DB "\$index->set_parameter(\$params);\n";
	print DB "my (\$masshash, \$peptidehash, \$proteinhash)=\$index->create_index(\$fasta_path);\n";

	print DB "my \$max_mass = \$params->{'max_peptide_mass'};\n";
	print DB "my \$min_mass = \$params->{'min_peptide_mass'};\n";
	
	print DB "my \$mass_range = int((\$max_mass-\$min_mass)/$num_mass_region);\n";
	print DB "my \$masshash_temp;\n";
	print DB "my \$peptidehash_temp;\n";

	print DB "for(my \$i=0;\$i<$num_mass_region;\$i++)\n";
	print DB "{\n";

##### Removed it because the range has been controlled by the max and min 8/1/2013 ############	
#	print DB "   next if (\$mass_range\*(\$i+1)<\$params->{'min_peptide_mass'});\n";          # 
###############################################################################################
		
	print DB "   my \$j=0;\n";	
	print DB "	foreach my \$mass (keys(\%\{\$masshash\}))\n";
	print DB "	{\n";

	print DB "		if(\$mass <(\$mass_range\*(\$i+1)+\$min_mass)\*1000 and \$mass > (\$mass_range\*\$i+\$min_mass)\*1000)\n";
	print DB "		{\n";
	print DB "			\$masshash_temp->{\$i}->{\$mass} =  \$masshash->{\$mass};\n";			
	print DB "			my \@pepIds = \@{\$masshash->{\$mass}};\n";
	print DB "			foreach my \$peptide (\@pepIds){\n";
	print DB "				\$j++;\n";
	print DB "				\$peptidehash_temp->{\$i}->{\$peptide}->{'seq'}=\$peptidehash->{\$peptide}->{'seq'};\n";
	print DB "				\$peptidehash_temp->{\$i}->{\$peptide}->{'proteinID'}=\$peptidehash->{\$peptide}->{'proteinID'};\n";
	print DB "			}\n";				
	print DB "		}\n";
	print DB "	}\n";
	print DB "next if (!defined(\$masshash_temp->{\$i}));\n";
	print DB "my \$nbMass = scalar keys %{\$masshash_temp->{\$i}}; \n";
	print DB 'my $mass_out = $fasta_path ."_$nbMass"."_$j"."_mass_$i";'."\n";
	print DB "	store(\$masshash_temp->{\$i},\$mass_out);\n";
	print DB "	my \$pep_out = \$fasta_path . \"_pep_\$i\";\n";
	print DB "	store(\$peptidehash_temp->{\$i},\$pep_out);\n";
	print DB "}\n";
	print DB "my \$pro_out = \$fasta_path . \"_pro\";\n";
	print DB "store(\$proteinhash,\$pro_out);\n";	
	close(DB);	
	
}

sub make_partialidx_script
{
	my ($self,$dir) = @_;
	my $lib = $self->get_library_path();
	open(PIDX,">$dir/Create_Partial_Idx.pl");
print PIDX <<EOF;
#!/usr/bin/perl -I $lib

use Getopt::Long;	

use Cwd;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::Params;
use Spiders::BuildIndex;
	
my (\$m,\$job_num,\$dta_path,\$database_path,\$mass_index,\$peptide_index,\$protein_index,\$databasename,\$num_pro_per_job,\$prot_job_num,\$inclusion_decoy,\$digestion);
GetOptions('-m=i'=>\\\$m,
		'-job_num=i'=>\\\$job_num,
		'-dta_path=s'=>\\\$dta_path,
		'-database_path=s'=>\\\$database_path,
		'-mass_index=s'=>\\\$mass_index,
		'-peptide_index=s'=>\\\$peptide_index,
		'-protein_index=s'=>\\\$protein_index,
		'-databasename=s'=>\\\$databasename,
		'-num_pro_per_job=i'=>\\\$num_pro_per_job,
		'-prot_job_num=i'=>\\\$prot_job_num,
		'-inclusion_decoy=i'=>\\\$inclusion_decoy,
		'-digestion=s'=>\\\$digestion,	
		);	

	if(\$inclusion_decoy)
	{
		\$num_pro_per_job = +2*\$num_pro_per_job;
	}
				
	#Get the stating indexes
	my (\$prev,\$j) = GetPreviousFilesMass(\$m,\$job_num); 	
	my \$l=0;	
		

	my (\$mass_hash,\$peptide_hash,\$protein_hash);
	my \%temp_peptideseq={};
	for(\$i=1;\$i<=\$job_num;\$i++)
	{
		## Get the mass file name
		my \$filename = "\${database_path}/temp_\${i}_\${databasename}_*_*_mass_\$m";
		my \$cmd = qx(ls -1 \$filename);
		next if(\$cmd eq \"\");
		my \@cmd_array = split("\\n",\$cmd);
		#my \$mass_out = "\${database_path}/temp_\${i}_\${databasename}_mass_\$m";
		my \$mass_out = \$cmd_array[0];
		my \$mass_retrieve = retrieve(\$mass_out);
		my \$pep_out = "\${database_path}/temp_\${i}_\${databasename}_pep_\$m";
		my \$pep_retrieve = retrieve(\$pep_out);
				
		print "\${database_path}/temp_\${i}_\${databasename}_mass_\$m\n";

		
		######### temporary peptide seq hash ######

		foreach my \$mass_mw (keys \%\$mass_retrieve)
		{	
			foreach \$k (\@{\$mass_retrieve->{\$mass_mw}})
			{
			
    #   reduce the database size #######################################
				my \$key = \$temp_peptideseq{\$pep_retrieve->{\$k}->{'seq'}};
				if (defined(\$key) and \$key=~/\\d+/)
				{
					push (\@{\$peptide_hash->{\$key}->{'proteinID'}},\$pep_retrieve->{\$k}->{'proteinID'}+(\$i-1)*\$num_pro_per_job);											
					next;
				}
				\$j++;				
				\$temp_peptideseq{\$pep_retrieve->{\$k}->{'seq'}}=\$j;
	##########################################################################################			

				push \@{\$mass_hash->{\$mass_mw}},(\$j);
				\$peptide_hash->{\$j}->{'seq'}=\$pep_retrieve->{\$k}->{'seq'};
				push (\@{\$peptide_hash->{\$j}->{'proteinID'}},\$pep_retrieve->{\$k}->{'proteinID'}+(\$i-1)*\$num_pro_per_job);							

####### the protein ID has to be adjusted  					
####### *2 is because ID comprise of decoy (generated by makeDB module) and target 
			}
		}
		
		if(\$m==\$prot_job_num)
		{
			my \$pro_out = "\${database_path}/temp_\${i}_\${databasename}_pro";
			my \$pro_retrieve = retrieve(\$pro_out);		
			\$last_pro_num = \$l;
			foreach my \$proteinid (keys \%\$pro_retrieve)
			{
				\$l++;
				\$protein_hash->{\$proteinid+\$last_pro_num}=\$pro_retrieve->{\$proteinid};
			}
		}		
	}
	my \$index = new Spiders::BuildIndex();
	\$index->set_dta_path(\$dta_path);
	\$index->set_digestion(\$digestion);
#	\$index->create_mass_index_file(\$mass_hash,\$mass_index,\$prev_run);
	
	\$index->create_masspep_index_file(\$mass_hash,\$peptide_hash,\$mass_index,\$peptide_index,\$prev);
	
	undef \$mass_hash;
	undef \$peptide_hash;
	if(\$m==\$prot_job_num)
	{	
		\$index->create_protein_index_file(\$protein_hash,\$protein_index);
		undef \$protein_hash;
	}		

sub GetPreviousFilesMass{
	my (\$m,\$job_num) = \@_;
	
	my \$prevn =0;
	my \$prevpep=0;
	my \@files = glob("\${database_path}/temp_\*_\${databasename}_\*_\*_mass_\*");
	
	for(my \$j=0;\$j<\$m;\$j++)
	{
		for(my \$i=1;\$i<\$job_num;\$i++)
		{
			for(my \$k=1;\$k<\$#files;\$k++)
			{
				
		#	my \$filename = "\${database_path}/temp_\${i}_\${databasename}_\*_\*_mass_\$j";
		#	my \$cmd = qx(ls -1 \$filename);
		#	next if (\$cmd eq \"\");
		#	my \@cmd_array = split("\\n",\$cmd);
		#	\$cmd = \$cmd_array[0];
				if(\$files[\$k] =~ /.*temp_\${i}_\${databasename}_(\\d+)_(\\d+)_mass_\$j/)
				{			
					\$prevn += \$1;		
					\$prevpep+= \$2;
				}
			}
		}
	}	
	
	return(\$prevn,\$prevpep);
}
EOF

}

sub create_job_files
{
	my ($self,$dta_path)=@_;
	
#	open(DECHARGE,">job.pl");
}

1;
