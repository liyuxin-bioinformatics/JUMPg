#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Search

package Spiders::Search;

######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Search database	                      ****    #
#       ****					                      ****    #
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #
#       ****all rights reserved.		              ****    #
#       ****xusheng.wang@stjude.org		              ****    #
#       ****					                      ****    #
#       ****					                      ****    #
#       **************************************************    #
###############################################################

use Spiders::Params;
use Spiders::MassUtils;
use Spiders::Hypergeometric;
use Spiders::BuildIndex;
use Spiders::MathUtils;


require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(findTag readTags exportMatches);
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 2.01;

######### version 1.02 ############
# To speed the search algorithm, it is changed to use the index file to search the database 
######### version 1.04 ############
# change the search module to allow search for low resolution data 
# 1. added the fragmented ion mass tolerance
######## Version 12.1.0 ##############
#add the amino acid to both sides of peptide  
######################################
######## Version 12.1.0 ##########
# Fix the hyper p value issue 

sub new{

	my ($class, %arg) = @_;
	my $self = {};
	bless $self, $class;
	return $self
}

sub set_H_value
{
	my ($self,$c_value)=@_;
	$self->{_H_value}=$c_value;	
}

sub get_H_value
{
	my $self=shift;
	if(!defined($self->{_H_value}))
	{
		$self->{_H_value}=1.007276466812;
	}
	return $self->{_H_value};
}

sub set_C_value
{
	my ($self,$c_value)=@_;
	$self->{_C_value}=$c_value;	
}

sub get_C_value
{
	my $self=shift;
	if(!defined($self->{_C_value}))
	{
		$self->{_C_value}=1.0033548378;
	}
	return $self->{_C_value};
}

sub set_pho_neutral_loss
{
	my ($self,$pho_neutral_loss)=@_;
	$self->{_pho_neutral_loss}=$pho_neutral_loss;	
}

sub get_pho_neutral_loss
{
	my $self=shift;
	if(!defined($self->{_pho_neutral_loss}))
	{
		$self->{_pho_neutral_loss}=0;
	}
	return $self->{_pho_neutral_loss};
}

sub readTags{

	my ($self, $f) = @_;
	
	my $tagHash = {};
	open F, $f or die;	
	while(<F>){
		if ($. > 1){
			chomp;
			# 6 columns: (1) Scan number, (2) precursor mass; (3) Tag sequences; (4) Side mass; (5) Rank p value; (6) Hyper p value 
			my @cols = split(/\t/,$_);
			my $h = {'scanNum'=>$cols[0],
					'precMass'=>$cols[1],
					'tagSeq'=>$cols[2],
					'sideMass'=>$cols[3],
					'rankP'=>$cols[4],
					'hyperP'=>$cols[5],
					'combP'=>$cols[6],
					};
			$tagHash->{$cols[0]} = $h;
		}		
	}
	close F;

	return $tagHash;	
}

sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}

sub set_tag
{
	my ($self,$tag)=@_;
	$self->{'tag'}=$tag;
}

sub get_tag
{
	my ($self)=@_;
	return $self->{'tag'};
}

sub set_tag_number
{
	my ($self,$tag_num)=@_;
	$self->{'tag_num'}=$tag_num;
}

sub get_tag_number
{
	my ($self)=@_;
	return $self->{'tag_num'};
}

sub set_precMass
{
	my ($self,$precmass) = @_;
	$self->{'precMass'} = $precmass;
}

sub get_precMass
{
	my ($self) = @_;
	return $self->{'precMass'};
}

sub set_precCharge
{
	my ($self,$prec_charge) = @_;
	$self->{'precCharge'} = $prec_charge;
}

sub get_precCharge
{
	my ($self) = @_;
	return $self->{'precCharge'};
}

sub set_scanNum
{
	my ($self,$scan) = @_;
	$self->{'scanNum'} = $scan;		
}

sub get_scanNum
{
	my $self=shift;
	return $self->{'scanNum'};
}

sub set_dtafile
{
	my ($self,$scan) = @_;
	$self->{'dtafile'} = $scan;		
}

sub get_dtafile
{
	my $self=shift;
	return $self->{'dtafile'};
}

sub get_tagSeq
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'tagSeq'} = $tag->{'tagSeq'};	
	return $self->{'tagSeq'};
}

sub get_sideMass
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'sideMass'} = $tag->{'sideMass'};
	return $self->{'sideMass'};
}

sub set_database
{
	my ($self,$database) = @_;
	$database=~s/\.mdx//;
	$self->{'database'} = $database;
}

sub get_database
{
	my ($self) = @_;
	return $self->{'database'};
}

sub set_masshash
{
	my ($self,$masshash) = @_;
	$self->{'masshash'} = $masshash;
}

sub get_masshash
{
	my ($self) = @_;
	return $self->{'masshash'};
}

sub set_peptidehash
{
	my ($self,$peptidehash) = @_;
	$self->{'peptidehash'} = $peptidehash;
}

sub get_peptidehash
{
	my ($self) = @_;
	return $self->{'peptidehash'};
}

sub set_proteinhash
{
	my ($self,$proteinhash) = @_;
	$self->{'proteinhash'} = $proteinhash;
}

sub get_proteinhash
{
	my ($self) = @_;
	return $self->{'proteinhash'};
}


sub get_tagPvalue
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'rankP'} = $tag->{'rankP'};
	return $self->{'rankP'};
}

sub set_mass_tolerance
{
	my ($self,$mass_tolerance)=@_;
	$self->{'mass_tolerance'} = $mass_tolerance;
}

sub get_mass_tolerance
{
	my $self=shift;
	return $self->{'mass_tolerance'};
}

sub set_mass_tolerance_units
{
	my ($self,$mass_tolerance_units)=@_;
	$self->{'mass_tolerance_units'} = $mass_tolerance_units;
}

sub get_mass_tolerance_units
{
	my $self=shift;
	return $self->{'mass_tolerance_units'};
}



sub set_frag_tolerance
{
	my ($self,$mass_tolerance)=@_;
	$self->{'frag_tolerance'} = $mass_tolerance;
}

sub get_frag_tolerance
{
	my $self=shift;
	return $self->{'frag_tolerance'};
}

sub get_AA_mass
{
	my $self=shift;
	
	my $massutil = new Spiders::MassUtils();
	my $params = $self->get_parameter();	
	$massutil->set_parameter($params);
	
	my $AA_mass = $massutil->get_AA_mass();
	return $AA_mass;
}


	

sub SearchMass{

	my ($self,$research) = @_;
	my %MatchedMass;
	
#	my %proteinid;
#	my %precursor;	
#	my %tagGroups;
#	my %threotical_MH;
	
	my $index = new Spiders::BuildIndex();
	my $database = $self->get_database();
	
	my $mass_index_file = $database . ".mdx";
	my $peptide_index_file = $database . ".pdx";	
	my $protein_index_file = $database . ".prdx";	

	my $error = $self->get_mass_tolerance();
	my $error_unit = $self->get_mass_tolerance_units();
	my $ms1_pho_loss_num = $self->get_pho_neutral_loss();
	my $precmass = $self->get_precMass();
	my $C = $self->get_C_value();	
	
#	my $max_modif = 0;
######### version 1.1 using index to search the value #####################	
# get precursor mass 
# if there is no modification 
####### ----------------------o-------------------------o------------
####### ---------------------start---------------------end------------
#############################=============selected=======############
##### Version 12.1.0 both ppm and th can be used ################
	if($error_unit==2)
	{
		$error = $error*$precmass/1000000;
	}

	my $search_loop = 0;
	if($research)
	{
		$search_loop = $self->get_isotopic_distribution($precmass);
		$search_loop = 1 if ($search_loop > 1);
	}
	 
	for(my $loop = -$search_loop; $loop<=$search_loop;$loop++)
	{
	
		my $prec_start = int(($precmass-$error - $C*$loop)*1000);
		my $prec_end = int(($precmass+$error - $C*$loop)*1000);

	#	$self->search_core(\%proteinid,\%tagGroups,\%threotical_MH,$prec_start,$prec_end);

		for(my $prec_mass = $prec_start; $prec_mass <= $prec_end; $prec_mass++)
		{
			my $massId = $index->readMassIndex($mass_index_file,$prec_mass);

			if($massId!=-1)
			{
				my $peptidehash =$index->getAssociatedPeptides($massId,$peptide_index_file);

				foreach $PeptideID (keys %{$peptidehash})
				{	
							
					my $mod;

					my %mod_pos=();
					
					my $nomod_pep = $peptidehash->{$PeptideID}->{'seq'};

					my @seq_array = split(/\./,$nomod_pep);
					$nomod_pep = $seq_array[1];

					if($ms1_pho_loss_num)
					{
						my $num_S = 0;
						while($nomod_pep =~ /S\W+/g)
						{
							$num_S++;
						}
						my $num_T = 0;
						while($nomod_pep =~ /T\W+/g)
						{
							$num_T++;
						}
						my $num_ST = $num_S + $num_T;	

						next if($num_ST<$ms1_pho_loss_num);
					}						
			
					while($nomod_pep =~ /[^a-zA-Z]+/g)
					{
						$mod_pos{$-[0]}=1;
						$nomod_pep =~ s/[^a-zA-Z]+//;
					}				
					if((scalar keys %mod_pos)==0)
					{
	############### 8/16/2013 fixed a length issue 				
						my $length=length($peptidehash->{$PeptideID}->{'seq'});
						$length=length($nomod_pep);
						
						#$nomod_pep=$peptidehash->{$PeptideID}->{'seq'};
						$mod = ":" x ($length+1);
					}
					else
					{
	############## if the peptide has modification, get the modification info ##################	
	# Fix a bug: K.QIVWKYCGR.M@			
	#					my $nomod_pep = $peptidehash->{$PeptideID}->{'seq'};
	#					my @seq_array = split(/\./,$nomod_pep);
	#					$nomod_pep = $seq_array[1];			
						my @nomodseq = split(//,$nomod_pep);
					
						my $length=length($nomod_pep);

						for(my $i=0;$i<=$length;$i++)
						{
							if($mod_pos{$i+1})
							{
								$mod .= ":";
								$mod .= $nomodseq[$i];
							}
							else
							{
								$mod .= ":";
							}
						}

					}

	#				print $peptidehash->{$PeptideID}->{'seq'},"\t",$nomod_pep,"\t",$mod,"\n";
	#				my $peptide_mass = 0;				
					
	#				$peptide_mass = $massutil->get_peptide_mass($nomod_pep,$mod,"N");
					
		#			print $prec_mass,"\t",$peptide_mass,"\n";
		
	########### In version2.0.1, we merged all individual hashes into one MatchedMass hash 	
					$MatchedMass{$PeptideID}{$mod}{'proteinID'} = $peptidehash->{$PeptideID}->{'proteinID'};
					$MatchedMass{$PeptideID}{$mod}{'origseq'} = $peptidehash->{$PeptideID}->{'seq'};
					$MatchedMass{$PeptideID}{$mod}{'seq'} = $seq_array[1];
					$MatchedMass{$PeptideID}{$mod}{'nomod_seq'} = $nomod_pep;
					$MatchedMass{$PeptideID}{$mod}{'theoretical_MH'} = $prec_mass;	

	###################################################				
	#				$proteinid{$PeptideID}{$mod} = $peptidehash->{$PeptideID}->{'proteinID'};			
	#				$tagGroups_orig{$PeptideID}{$mod}=$peptidehash->{$PeptideID}->{'seq'};
	#				$tagGroups{$PeptideID}{$mod}=$nomod_pep;
	########### in version 2: we used $prec_mass as the theorectical mass even it is a rough value ############				
	#				$threotical_MH{$PeptideID}{$mod}=$peptide_mass;
	#				$threotical_MH{$PeptideID}{$mod} = $prec_mass;
				}			
			}
		}
	}	
	$self->{'matched_prec_peptide_num'} = scalar keys %MatchedMass;
	return \%MatchedMass;
}

####### In version 2.0.1, we change "findTags" function into two independent functions: SearchMass and SearchTag. 
####### The advantage of sperating into two functions is only search once for mass no matter how many tags will be used for searching, where is especially useful for low resolution.
 

sub SearchTag
{	
	my ($self,$MatchedMassRef) = @_;

	my $matchResults;
	
	my $index = new Spiders::BuildIndex();
	my $database = $self->get_database();	
	my $protein_index_file = $database . ".prdx";	
	
	my $massutil = new Spiders::MassUtils();
	my $params = $self->get_parameter();	
	$massutil->set_parameter($params);
	my $tags = $self->get_tag();
	my $frag_tolerance = $self->get_frag_tolerance();
	
	my $charge = $self->get_precCharge();
	
	my $AA_mass = $self->get_AA_mass();	
	my $parameter = new Spiders::Params();
	
	my ($dynamic_mod,$largest_mass) = $parameter->get_dynamic_modifications($params);

	my $static_mod = $parameter->get_static_modification($params);	
	
	my %MatchedMass = %$MatchedMassRef;
	
	foreach my $pep (keys %MatchedMass)
	{
		foreach $mod (keys %{$MatchedMass{$pep}})
		{
########################## in order to include the modified tag, change the nomod_seq to orig_seq
					my $pepseq = $MatchedMass{$pep}{$mod}{'nomod_seq'};
#	my $pepseq_core = $MatchedMass{$pep}{$mod}{'nomod_seq'};
					my $pepseq_core = $MatchedMass{$pep}{$mod}{'seq'};
										
############################################################################
				# discard those peptides with charge = 1 and length <=5

					next if($charge==1 and length($tags->{'tagSeq'})<=3 and length($pepseq_core)<=8);
		#			next if(length($pepseq_core)<=5);					
######################################################################
#			my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$mod,$tags,$error);

			my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$pepseq_core,$mod,$tags,$frag_tolerance,$AA_mass,$dynamic_mod,$static_mod);
			
			if ($sideMassMatchResult->{'lSide_match'} || $sideMassMatchResult->{'rSide_match'})
			{

				my $protein_id = $MatchedMass{$pep}{$mod}{'proteinID'};
				
				#my ($proteinName,$proteinDesc) = $index->getProtein($protein_id,$protein_index_file);
							
				$matchResults->{$pep}->{$mod}->{'scanNum'} = $tags->{'scanNum'};
				$matchResults->{$pep}->{$mod}->{'precMass'} = $tags->{'precMass'};
				$matchResults->{$pep}->{$mod}->{'theoretical_MH'} = $MatchedMass{$pep}{$mod}{'theoretical_MH'};
				$matchResults->{$pep}->{$mod}->{'tagSeq'} = $tags->{'tagSeq'};

				$matchResults->{$pep}->{$mod}->{'tag_rank_num'} = $tags->{'tag_rank_num'};	
				
				$matchResults->{$pep}->{$mod}->{'sideMass'} = $tags->{'sideMass'};			
				$matchResults->{$pep}->{$mod}->{'pepseq'} = $pepseq;
				
				$matchResults->{$pep}->{$mod}->{'proteinid'} = $protein_id;

												
				$matchResults->{$pep}->{$mod}->{'tagSeqMatch'} = $tagSeqMatchResult;
				$matchResults->{$pep}->{$mod}->{'lsideMassMatch'} = $sideMassMatchResult->{'lSide_match'};
				$matchResults->{$pep}->{$mod}->{'rsideMassMatch'} = $sideMassMatchResult->{'rSide_match'};
				$matchResults->{$pep}->{$mod}->{'lthreot_Mass'} = sprintf("%.4f",$sideMassMatchResult->{'lthreot_Mass'});
				$matchResults->{$pep}->{$mod}->{'rthreot_Mass'} = sprintf("%.4f",$sideMassMatchResult->{'rthreot_Mass'});

				$matchResults->{$pep}->{$mod}->{'matchPeptide'} = $pepseq;

					
				$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'origseq'};

				$matchResults->{$pep}->{$mod}->{'rank_p'} = $tags->{'rankP'};
				$matchResults->{$pep}->{$mod}->{'hyper_p'} = $tags->{'hyperP'};
				$matchResults->{$pep}->{$mod}->{'comb_p'} = $tags->{'combP'};					
			}
			
	
		}
	}
	my $num_peptides = scalar keys %$matchResults;

	
	if($num_peptides > 0)
	{
			return $matchResults;
	}
	else
	{
		return 0;
	}
}





sub matchtwoseq
{
	my ($self,$seq1,$seq2)=@_;
	my $matched=0;
	my $pseq = $seq1;
	$pseq =~ s/I/J/g;
	$pseq =~ s/L/J/g;

	my $tseq = $seq2;
	$tseq =~ s/I/J/g;
	$tseq =~ s/L/J/g;
	my $rev_tseq=reverse($tseq);
	if($pseq=~/$tseq/)
	{
		$matched=1;
	}
	elsif($pseq=~/$rev_tseq/)
	{
		$matched=1;
	}
	else
	{
		$matched=0;
	}
	return $matched;
}

sub matchTagSeq{

	my ($self,  $pepseq, $pepseq_mod, $mod, $tags, $error,$AA_mass,$dynamic_mod,$static_mod) = @_;


	my $pseq = $pepseq_mod;
#	my $pseq = $pepseq_mod;
	$pseq =~ s/I/J/g;
	$pseq =~ s/L/J/g;
	
#	$pepseq =~ s/I/J/g;
#	$pepseq =~ s/L/J/g;
	
	
	my $tseq = $tags->{'tagSeq'};
	$tseq =~ s/I/J/g;
	$tseq =~ s/L/J/g;
	
	my $tagSeqMatch = 0;
	my $sideMassMatch;
	my $pepsideMass = 0;
	my $whole = -1;
	my $lPart = -1;
	my $rPart = -1;
	my $removed_seq = "";

	while (1) 
	{	

	#	$whole = index($pseq,$tseq,$whole+1);
		$whole = index($pseq,$tseq,$whole+1);

		if($whole < 0)
		{
			$sideMassMatch->{'lSide_match'} = 0;
			$sideMassMatch->{'rSide_match'} = 0;
			last;			
		}
		$tagSeqMatch = 1;
		
		$sideMassMatch = $self->checkSideMass($whole,$tags,$pseq,$mod,$error,$removed_seq,$AA_mass,$dynamic_mod,$static_mod);
		if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0)
		{
			return ($tagSeqMatch,$sideMassMatch);
			last;
		}		
	}

	$whole = -1;
	if($sideMassMatch->{'lSide_match'}==0 and $sideMassMatch->{'rSide_match'}==0)
	{
		my $rev_tseq = reverse ($tseq);
#		$rev_tseq =~ s/([^\w])([\w])/$2$1/g;
	
		while (1) {	
			$whole = index($pseq,$rev_tseq,$whole+1);
	
			last if($whole < 0);
			$sideMassMatch = $self->checkSideMass_rev($whole,$tags,$pseq,$mod,$error,$removed_seq,$AA_mass,$dynamic_mod,$static_mod);

			if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0)
			{
				$tagSeqMatch = 1;

				return ($tagSeqMatch,$sideMassMatch); 			
			}			
		}
	}
	

	
	if(!defined($sideMassMatch->{'lSide_match'}))
	{
		$sideMassMatch->{'lSide_match'}=0;
	}
	if(!defined($sideMassMatch->{'rSide_match'}))
	{
		$sideMassMatch->{'rSide_match'}=0;
	}
	
	return ($tagSeqMatch,$sideMassMatch);
	
}

sub checkSideMass{

	my ($self, $ind, $tag, $pseq, $mod, $error,$removed_seq,$AA_mass,$dynamic_mod,$static_mod) = @_;
	my $massutil = new Spiders::MassUtils();
	
	my $parameter = $self->get_parameter();
	$massutil->set_parameter($parameter);
	
######## N-term mass has to add H2O + H
	my $result;
	my $pepSideMass;

	my $N_term_Mass = 19.017806;
	my $lSide = substr($pseq, 0, $ind);
	my $rSide = substr($pseq, $ind);
	
	my $lSide_mod = substr($mod,0,$ind);
	
	$lSide_mod = $lSide_mod . ":";
	my $rSide_mod = substr($mod, $ind);

	my $tolerance_unit = 1;
	if(defined($parameter->{'frag_mass_tolerance_unit'}))
	{
		$tolerance_unit = $parameter->{'frag_mass_tolerance_unit'};
	}

	
	if($removed_seq ne "")
	{
	##### need check ??????##
	#	$tag->{'sideMass'} = $tag->{'sideMass'} - $massutil->get_peptide_mass($removed_seq,"","C") + 1.007825;
	}
####### because the tag was reverse, see tag.pm module, so the N and C terminus was reversed
	$lSide =~ s/[^a-zA-Z]+//g;
	$rSide =~ s/[^a-zA-Z]+//g;

	$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"N",$AA_mass,$dynamic_mod);

	$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"C",$AA_mass,$dynamic_mod);

	#print $pseq,"\t",$tag->{"tagSeq"},"\t",$mod,"\t",$lSide,"\t",$rSide,"\t",$lSide_mod,"\t",$rSide_mod,"\t",$tag->{'sideMass'},"\t",$result->{'lthreot_Mass'},"\t",$result->{'rthreot_Mass'},"\t",$error,"\n"; 
	#	print abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}),"\n";
### commented by xusheng; only error without multiple with $tag->{'sideMass'}
### 6/28/2013 initially change the $error into 2*$error, but changed back in this version 	

	my $error_Da = $error;
	if($tolerance_unit == 2)
	{
		$error_Da = $error * $result->{'lthreot_Mass'} / 1000000;
	}	


	if (abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ) <= $error_Da )
	{
		$result->{'lSide_match'} =  1;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};
	}
	else
	{
#		$result->{'lthreot_Mass'} = $self->calMW($lSide);

		$result->{'lSide_match'} = 0;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};	
	}
	if($tolerance_unit == 2)
	{
		$error_Da = $error * $result->{'rthreot_Mass'} / 1000000;
	}	
	if (abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}) <= $error_Da ){
		$result->{'rSide_match'} = 1 ;

		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}
	else{
		$result->{'rSide_match'} = 0 ;
		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}

	
	return $result;
}

sub checkSideMass_rev{

	my ($self, $ind, $tag, $pseq, $mod, $error,$removed_seq,$AA_mass,$dynamic_mod,$static_mod) = @_;
	my $massutil = new Spiders::MassUtils();
	
	my $parameter = $self->get_parameter();
	$massutil->set_parameter($parameter);
	
######## N-term mass has to add H2O + H
	my $result;
	my $pepSideMass;

## to match the side Mass, needs to add the length of tag because N and C terminus shift;  	
	$ind += length($tag->{'tagSeq'});
#####	
	my $lSide = substr($pseq, 0, $ind);
	my $rSide = substr($pseq, $ind);

	$lSide =~ s/[^a-zA-Z]+//g;
	$rSide =~ s/[^a-zA-Z]+//g;	
	my $lSide_mod = substr($mod,0,$ind);
	$lSide_mod = $lSide_mod . ":";	
	my $rSide_mod = substr($mod, $ind);
	if($removed_seq ne "")
	{
#		$tag->{'sideMass'} = $tag->{'sideMass'} - $massutil->get_peptide_mass($removed_seq,"","C") + 1.007825;
	}
####### because the tag was reverse, see tag.pm module, so the N and C terminus was reversed
########### ????????????????????????????????
#	$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"N") - 19.017806 + 1.007825;

#	$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"C") + 19.017806 - 1.007825;

	$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"N",$AA_mass,$dynamic_mod);

	$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"C",$AA_mass,$dynamic_mod);
	
### commented by xusheng; only error without multiple with $tag->{'sideMass'}


	my $tolerance_unit = 1;
	if(defined($parameter->{'frag_mass_tolerance_unit'}))
	{
		$tolerance_unit = $parameter->{'frag_mass_tolerance_unit'};
	}

	if($tolerance_unit == 2)
	{
		$error = $error * $tag->{'sideMass'} / 1000000;
	}	

#	print abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ),"\t",$error,"\n";
	if (abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ) <= $error ){

		$result->{'lSide_match'} =  1;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};
	}
	else
	{
#		$result->{'lthreot_Mass'} = $self->calMW($lSide);
		$result->{'lSide_match'} = 0;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};	
	}

	if (abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}) <= $error ){
		$result->{'rSide_match'} = 1 ;

		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}
	else{
		$result->{'rSide_match'} = 0 ;
		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}

	return $result;
}

############## Version 12.1.0 ######################
## Add the function allowing to search without Tag matches 

sub SearchWithoutTag
{	
	my ($self,$MatchedMassRef) = @_;

	my $matchResults;
	
	my $index = new Spiders::BuildIndex();
	my $database = $self->get_database();	
	my $protein_index_file = $database . ".prdx";	
	
	my $massutil = new Spiders::MassUtils();
	my $params = $self->get_parameter();	
	$massutil->set_parameter($params);
	my $tags = $self->get_tag();

	my $charge = $self->get_precCharge();	
	
	my %MatchedMass = %$MatchedMassRef;
	
	foreach my $pep (keys %MatchedMass)
	{
		foreach $mod (keys %{$MatchedMass{$pep}})
		{
			my $pepseq = $MatchedMass{$pep}{$mod}{'nomod_seq'};
	
############## 10/14/2013 #######################
#			Discard those peptides with charge = 1, No tag, and peptide length < 8
			my $length_peptide = length($pepseq);

			next if($charge==1 and $length_peptide <= 7);
# 			Discard those peptides with length <= 5
			next if($length_peptide<=6);
# 			Discard those peptides with charge = 2, No tag, and miscleavage > 1
	
###################################################
	

			
			my $protein_id = $MatchedMass{$pep}{$mod}{'proteinID'};

		
			$matchResults->{$pep}->{$mod}->{'scanNum'} = $tags->{'scanNum'};
			$matchResults->{$pep}->{$mod}->{'precMass'} = $tags->{'precMass'};
			$matchResults->{$pep}->{$mod}->{'theoretical_MH'} = $MatchedMass{$pep}{$mod}{'theoretical_MH'};
			$matchResults->{$pep}->{$mod}->{'tagSeq'} = "N/A";			
			$matchResults->{$pep}->{$mod}->{'sideMass'} = "";			
			$matchResults->{$pep}->{$mod}->{'pepseq'} = $pepseq;
			$matchResults->{$pep}->{$mod}->{'proteinid'} = $protein_id;

												
			$matchResults->{$pep}->{$mod}->{'tagSeqMatch'} = "";
			$matchResults->{$pep}->{$mod}->{'lsideMassMatch'} = "";
			$matchResults->{$pep}->{$mod}->{'rsideMassMatch'} = "";
			$matchResults->{$pep}->{$mod}->{'lthreot_Mass'} = "";
			$matchResults->{$pep}->{$mod}->{'rthreot_Mass'} = "";

			$matchResults->{$pep}->{$mod}->{'matchPeptide'} = $pepseq;

			$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'origseq'};
			$matchResults->{$pep}->{$mod}->{'rank_p'} = "";
			$matchResults->{$pep}->{$mod}->{'hyper_p'} = "";
			$matchResults->{$pep}->{$mod}->{'comb_p'} = "";
		}
	}
	my $num_peptides = scalar keys %$matchResults;
	if($num_peptides > 0)
	{
		return $matchResults;
	}
	else
	{
		return 0;
	}	
}

sub generate_peptide_theoritical_mass
{
## version 11.0.2
#	my ($self,$peptide,$modif)=@_;
	my ($self,$peptide,$modif,$peptidemass,$AA_mass)=@_;
#######
	
	my $parameter = $self->get_parameter();
	my $massutils = new Spiders::MassUtils();
	$massutils->set_parameter($parameter);

##########################################################	
##########Version 12.1.0########################	
# add the ion series function 	
	my $ion_series = $parameter->{'ion_series'};

	
	my @ions_array = split(/\s+/,$ion_series);
	if(scalar(@ions_array) != 9)
	{
		print "please put the right ion_series parameter!\n";
		exit;
	}
	my @ion_series_name = ('a', 'b', 'c', 'd', 'v', 'w', 'x', 'y', 'z','a++', 'b++', 'c++', 'd++', 'v++', 'w++', 'x++', 'y++', 'z++');
	my @ion_series_used;
	for (my $i=0;$i<$#ions_array;$i++)
	{
		if($ions_array[$i]==1)
		{
			push (@ion_series_used,$ion_series_name[$i]);
		}
	}
# add the charge state 
	my $charge = $self->get_precCharge();
	if($charge>=2)
	{
		for (my $i=0;$i<$#ions_array;$i++)
		{	
			if($ions_array[$i]==1)
			{
				push (@ion_series_used,$ion_series_name[$i+9]);
			}
		}
	}

####### 11/3/2014 #########################
	my $ion_loss = 0;
	if(defined($parameter->{'ion_losses_MS2'}))
	{
		$ion_loss = $parameter->{'ion_losses_MS2'};

		my @ions_loss_array = split(/\s+/,$ion_loss);
		my @temp_series = @ion_series_used;

		
		if($ions_loss_array[0]==1)
		{
			foreach (@temp_series)
			{
				my $new_loss = $_ . '-H2O';
				push (@ion_series_used,$new_loss);	
			}
		}
		if($ions_loss_array[1]==1)
		{
			my $num_S = $peptide =~ tr/S//;
			my $num_T = $peptide =~ tr/T//;
			
			next if(($num_S+$num_T/3)<0.6);		
			foreach (@temp_series)
			{
				my $new_loss = $_ . '-HPO3';
				push (@ion_series_used,$new_loss);	
			}
		}		
		if($ions_loss_array[2]==1)
		{
			my $num_S = $peptide =~ tr/S//;
			my $num_T = $peptide =~ tr/T//;
			
			next if(($num_S+$num_T/3)<0.6);			
			foreach (@temp_series)
			{
				my $new_loss = $_ . '-H3PO4';
				push (@ion_series_used,$new_loss);	
			}
		}
		if($ions_loss_array[3]==1)
		{
			foreach (@temp_series)
			{
				my $new_loss = $_ . '-NH3';
				push (@ion_series_used,$new_loss);	
			}
		}		
	}	

#	push (@ion_series_used,'immo');
# version 11.0.2
	$massutils->getFragmentMasses($AA_mass,$fragType,$dynamic_mod,$series,$loss,pept=>$peptide, modif=>$modif, fragTypes=>[@ion_series_used], spectrum=>\%spectrum);

#	$massutils->getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>[@ion_series_used], spectrum=>\%spectrum, peptidemass=>$peptidemass,AA_mass=>$AA_mass);
#	$massutils->getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>['b','y','b++','y++','immo'], spectrum=>\%spectrum);
####################end ####################################


## For different ions
	my $exp_mz = $self->get_exp_mz();
########################
	my $ion_scoring = 1;
	$ion_scoring = $parameter->{'ion_scoring'};
	my (@aa, @bb);
	my $match_ion;
	foreach (keys %{$spectrum{'mass'}{'term'}}){
	#	print $peptide,"\t",$modif,"\t",$_,"\n";
		push (@aa, @{$spectrum{'mass'}{'term'}{$_}});
		#foreach my $ions (@{$spectrum{'mass'}{'term'}{$_}})
		#{
		#	print $ions," ";
		#}
		#print "\n";
#########################################
		if($ion_scoring==2)
		{
			my ($matched, $matched_peaks) = $self->compare_theoritical_experiment($exp_mz,\@{$spectrum{'mass'}{'term'}{$_}});
			$match_ion->{$_} = 	$matched . "|" . $#{$spectrum{'mass'}{'term'}{$_}};
		}
#		print $_,"\t",$match_ion->{$_},"\n";
##############################################	
		
	}

	

#	foreach (keys %{$spectrum{'mass'}{'intern'}}){
#		@bb =  @{$_->{'masses'}};
#		print $_,"\t",@bb,"\n";
#	}

	#my @big = (@aa,@bb);

		
	return (\@aa,$match_ion);
	
}

sub compare_theoritical_experiment
{
	my ($self,$exp_mz,$theor_mz) = @_;
	my $parameter = $self->get_parameter();
	my $matched = 0;
	my $tolerance_input = $self->get_frag_tolerance();
	my $tolerance_unit = 1;
	if(defined($parameter->{'frag_mass_tolerance_unit'}))
	{                        
		$tolerance_unit = $parameter->{'frag_mass_tolerance_unit'};
	}

	
#### 6/28/2013 Fixed a bug: count multiple times because of similar theoretical mass generated from 1 charge and 2 charge	
	my %matched_hash={};

	$tolerance = $tolerance_input;
### version 11.0.2 to boot the speed

	for (my $i=0;$i<=$#{$theor_mz};$i++)
	{	
	
		my $theor_mz_value = $theor_mz->[$i];
		my $theor_mz_value_int = int($theor_mz_value);

		if($tolerance_unit == 2)
		{
			$tolerance = $tolerance_input * $theor_mz_value_int / 1000000;
		}
		my $theor_mz_value_min = int($theor_mz_value-$tolerance);
		my $theor_mz_value_max = int($theor_mz_value+$tolerance)+1;
	
		for (my $theor_mz_value_i = $theor_mz_value_min; $theor_mz_value_i < $theor_mz_value_max; $theor_mz_value_i++)
		{
			if(defined($exp_mz->{$theor_mz_value_i}))
			{
				foreach my $exp_mz_value (keys %{$exp_mz->{$theor_mz_value_i}})
				{
				
					if(($exp_mz_value+$tolerance)>$theor_mz_value and ($exp_mz_value-$tolerance)<$theor_mz_value)
					{
						next if($matched_hash{$theor_mz_value});

						$matched++;

						$matched_hash{$theor_mz_value}=1;
					}

				}
			}
		}			
	}

	return ($matched,\%matched_hash);
}

sub set_exp_mz
{
	my ($self,$exp_mz)=@_;
	$self->{'_exp_mz'} = $exp_mz;
}

sub get_exp_mz
{
	my $self=shift;
	return $self->{'_exp_mz'};
}

sub get_exp_mz_hash
{
	my $self=shift;
	my $exp_mz = $self->{'_exp_mz'};
	my %exp_hash=();
	foreach (@$exp_mz)
	{
		my $int_mz = int($_);
		$exp_hash{$int_mz}{$_}=1;
	}
	return \%exp_hash;
#	return $self->{'_exp_mz'};
}

sub set_exp_mz_int
{
	my ($self,$exp_mz_int)=@_;
	$self->{'_exp_mz_int'} = $exp_mz_int;	
}

sub get_exp_mz_int
{
	my $self=shift;
	return $self->{'_exp_mz_int'};
}

#### Version 12.1.0 in order to increase the speed #####################

sub get_peptide_matched
{
## version 11.0.2
#	my ($self,$peptide,$modif)=@_;
	my ($self,$peptide,$modif,$peptidemass,$AA_mass)=@_;

#	my ($theor_mz,$match_ion) = $self->generate_peptide_theoritical_mass($peptide,$modif);
	my ($theor_mz,$match_ion) = $self->generate_peptide_theoritical_mass($peptide,$modif,$peptidemass,$AA_mass);
########

	my $exp_mz = $self->get_exp_mz_hash();
	my ($matched, $matched_peaks) = $self->compare_theoritical_experiment($exp_mz,$theor_mz);

############ Version 12.1.0 ########################### 	
# add return theoretical ions and matched ions
	$self->{'Theor_ion_num'} = $#$theor_mz + 1; 
	$self->{'Matched_ion_num'} = $matched; 
	$self->{'Matched_peaks'} = $matched_peaks;
################### end ##############################	

	return ($matched,$match_ion);
}

sub get_peptide_pvalue
{
	my ($self,$matched,$total)=@_;
######### use different method for matching score #############
	my $parameter = $self->get_parameter();
	
###### Version 12.1.0 fix a bug: using frag tolerance instead of mass tolerance	
	my $mass_tolerance = $self->get_frag_tolerance();
	
	my $tolerance_unit = 1;
	if(defined($parameter->{'frag_mass_tolerance_unit'}))
	{
		$tolerance_unit = $parameter->{'frag_mass_tolerance_unit'};
	}
	
	my $exp_mz = $self->get_exp_mz();
	my $exp_mz_int = $self->get_exp_mz_int();
######### Version 12.1.0 fix a bug: for low resolution, using exp mass range rather than theoretical mass range	
#	my $total_location = $total/($mass_tolerance*2);
	
	if($tolerance_unit == 2)
	{
### use the average mass to convert the unit	
		$mass_tolerance = $mass_tolerance *($exp_mz->[$#$exp_mz] + $exp_mz->[0]) / 2000000;
	}			
	
######### Trim the data 
#########  exp -----------+-----------------------------+-------------
#########  the --------------+---------------------+------------------                	
	my $total_number = int($exp_mz->[$#$exp_mz] - $exp_mz->[0])/($mass_tolerance*2);
#	my $total_number = $total_location > $total_location_alt ? $total_location : $total_location_alt; 

	my $matched_peaks = $self->{'Matched_peaks'};

	my $exp_mass_num = scalar (@$exp_mz);
	
	my $hyper = new Spiders::Hypergeometric();
	my $wilcox_test = Spiders::WilcoxonRankSum->new();	
	my $log_peptide_pvalue = 1;
	
	if($parameter->{'matching_method'} eq "hyper_p")
	{
#		my $peptide_pvalue=$hyper->hypergeom($total_number,$total,$exp_mass_num,$matched);
	

	my $peptide_pvalue=$hyper->Hypergeometric($total_number,$exp_mass_num,$total,$matched);	
		$log_peptide_pvalue = sprintf("%.6f",-log($peptide_pvalue)/log(10));
	}
	elsif($parameter->{'matching_method'} eq "rank_p")
	{
		my @sel_array;
		my @unsel_array;

		foreach my $mz (keys %$exp_mz_int)
		{
			if($matched_peaks->{$mz})
			{
				push (@sel_array,$exp_mz_int->{$mz});
			}
			else
			{
				push(@unsel_array,$exp_mz_int->{$mz});
			}
		}

		$wilcox_test->load_data(\@sel_array, \@unsel_array);
		my $prob = $wilcox_test->probability();
		$log_peptide_pvalue = sprintf("%.6f",-log($prob)/log(10));
	
	}
	elsif($parameter->{'matching_method'} eq "comb_p")
	{
		my @sel_array;
		my @unsel_array;

		foreach my $mz (keys %$exp_mz_int)
		{
			if($matched_peaks->{$mz})
			{
				push (@sel_array,$exp_mz_int->{$mz});
			}
			else
			{
				push(@unsel_array,$exp_mz_int->{$mz});
			}
		}


		$wilcox_test->load_data(\@sel_array, \@unsel_array);
		my $prob = $wilcox_test->probability();

		my $log_peptide_rankp = sprintf("%.6f",-log($prob)/log(10));
		
#		my $peptide_pvalue=$hyper->hypergeom($total_number,$total,$exp_mass_num,$matched);

	my $peptide_pvalue=$hyper->Hypergeometric($total_number,$exp_mass_num,$total,$matched);	
		my $log_peptide_hyperp = sprintf("%.6f",-log($peptide_pvalue)/log(10));		

		$log_peptide_pvalue = ($log_peptide_rankp + $log_peptide_hyperp);
	}
	else
	{
		print "please use the right method for calculating matching score!\n";
	}
	
	return ($log_peptide_pvalue);
}

sub mergeresults
{
	my ($self,$result1,$result2,$tag) = @_;
	foreach (keys %$result2)
	{

		$result1->{$tag}->{$_} = $result2->{$_};
	}
	return $result1;
}

sub writingHeader
{
	my ($self,$outfile) = @_;
	my $parameter = $self->get_parameter();
## multiple tag #######	
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";

	print OUTPUT "\nspTag version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	my $dtafile = $self->get_scanNum();
	$dtafile =~ s/tag/dta/;
	print OUTPUT "DTA file = ", $dtafile,"\n";
	print OUTPUT "Database =", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass =",$self->get_precMass(),"\n\n";
	close(OUTPUT);
	
}

sub WriteResults
{
	my ($self,$outfile,$matches,$tag_rank_num,$prec_peak_int,$ms2_signal_noise_ratio) = @_;

	my $database = $self->get_database();
	my $AA_mass = $self->get_AA_mass();
	
	my $protein_index_file = $database . ".prdx";	
	my $index = new Spiders::BuildIndex();
	my $massutils = new Spiders::MassUtils();
	my $params = $self->get_parameter();
	$massutils->set_parameter($params);
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";
	my $parameter = $self->get_parameter();
#	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");	
	my %mod_symbol = ('M'=>"@",'S'=>"#",'T'=>"%",'Y'=>"*",'G'=>"^",'K'=>"&",'D'=>'?','A'=>'~','Q'=>'!','P'=>"(",'E'=>")",'E'=>"{",'V'=>"}","V"=>"[","H"=>"]","C"=>":","F"=>",","I"=>';',"L"=>',',"R"=>"<","N"=>">","W"=>"'");
	
## multiple tag #######	


	print OUTPUT "\nJUMP version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon,$mday,$hour,$min,$sec;	
	my $dtafile = $self->get_dtafile();
#	$dtafile =~ s/tag/dta/;
	print OUTPUT "DTA file = ", $dtafile,"\n";
	
	print OUTPUT "Database = ", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass = ",$self->get_precMass(),"   ", "Percentage of precursor peak intensity = ", sprintf("%0.2f",$prec_peak_int*100),"\%\n";
	print OUTPUT "Precursor matched peptides = ", $self->{'matched_prec_peptide_num'},"\, ","Peptides mass tolerance = ", $self->get_mass_tolerance(), "  Peptides mass tolerance units = $parameter->{'peptide_tolerance_units'}\n";
	print OUTPUT "MS2 signal noise ratio = ",sprintf("%0.2f",$ms2_signal_noise_ratio),"\, ","Tag number = ",$self->get_tag_number(),"\,  ","Tag ranking method = ", $parameter->{'tag_select_method'},"\n"; 
	
	print OUTPUT "ion series ABCDVWXYZ: ", $parameter->{'ion_series'},"\n";
	print OUTPUT "ion losses H2O H3PO4 Ammonia: ", $parameter->{'ion_losses_MS2'},"\n";
	
	our ($AA_mass,$dynamic_mod,$fragType,$series,$loss) = $massutils->get_fragment_parameter();
###### output modification information ###############################
	my $i=0;
	my $mod_symbol_used="";
	foreach (keys %$parameter)
	{
		if($_ =~/add_/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,"=",$parameter->{$_},"  ";
		}
		if($_ =~/dynamic_(\w+)/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,$mod_symbol{$1},"=",$parameter->{$_},"  ";
			$i++;
		}

#		if($_ =~/dynamic_/ and $parameter->{$_} != 0)
#		{
#			print OUTPUT $_,$mod_symbol[$i],"=",$parameter->{$_},"  ";
#			$i++;
#		}
	}

	print OUTPUT "\n\n";
#	print OUTPUT "Tag = ", $self->get_tagSeq(),"\n";
#	my $tag_e_value = sprintf("%.2f",$self->get_tagPvalue());
#	print OUTPUT "Tag E_value = ", $tag_e_value,"\n","Tag Rank Number = ",$tag_rank_num,"\n","Tag Side Mass = ",$self->get_sideMass(),"\n"; 
	my $i=0;
	my %sort_results;
	my %comb_results;
	my %peptide_matched_ratio;
	my %ions;
	my %matched_ions;
## if it is the same pep, same matched and thereotical ions, skip	
	my %matched_pep_ions;
	
	foreach my $tag (keys %{$matches})
	{
	
		foreach my $pep (sort keys %{$matches->{$tag}})
		{	
			my $j=1;
			foreach my $mod (keys %{$matches->{$tag}->{$pep}})
			{	
					next if ($matches->{$tag}->{$pep}->{$mod}->{'tagSeqMatch'} < 1);
		#		print $matches->{$pep}->{$mod}->{'lsideMassMatch'},"\t",$matches->{$pep}->{$mod}->{'rsideMassMatch'},"\n";
	#			if ((abs($matches->{$pep}->{$mod}->{'lthreot_Mass'} - $self->get_sideMass()) < 0.02) || (abs($matches->{$pep}->{$mod}->{'rthreot_Mass'} - $self->get_sideMass()) < 0.02))
	#			{
					my $peptideseq = $matches->{$tag}->{$pep}->{$mod}->{'pepseq'};
					
					my $pep_pvalue = 0;
				#	if(!defined($peptide_matched_ratio{$pep}))
				#	{					
#						my $pep_pvalue = $self->get_peptide_pvalue($peptideseq,$mod) + $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}};
####### Version 12.1.0; increase search speed 

##################### version 11.0.2 
# in order to increase the speed
# pass the peptide mass to the subroutine
#
#						my ($matched,$match_ion) = $self->get_peptide_matched($peptideseq,$mod);

						my ($matched,$match_ion) = $self->get_peptide_matched($peptideseq,$mod,$matches->{$tag}->{$pep}->{$mod}->{'theoretical_MH'}/1000,$AA_mass);

#######################
						$ions{$pep} = "$self->{'Matched_ion_num'}\/$self->{'Theor_ion_num'}";

						
						if(!defined($matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}}))
						{
							$matched_ions{$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'}} = $match_ion;
							next if ($matched<3);
													
	############### seperately scoring for different ion types	##################					
							my $ion_scoring = 1;
							$ion_scoring = $parameter->{'ion_scoring'};

							if($ion_scoring == 1)
							{
								$pep_pvalue = sprintf("%0.2f",$self->get_peptide_pvalue($self->{'Matched_ion_num'}, $self->{'Theor_ion_num'}));
							}
							elsif($ion_scoring == 2)
							{
								my %merge_ion_losses;
								foreach my $ion_type (keys %$match_ion)
								{
									my ($ion,$loss)=split(/-/,$ion_type);
									my ($matched_num,$theor_num) = split(/\|/,$match_ion->{$ion_type}); 
									$merge_ion_losses{$ion}{'matched_ion_num'} += $matched_num;
									$merge_ion_losses{$ion}{'Theor_ion_num'} += $theor_num;								
								}
								foreach my $merged_ion (keys %merge_ion_losses)
								{
									next if($merge_ion_losses{$merged_ion}{'matched_ion_num'}==0);
							
									$pep_pvalue += sprintf("%0.2f",$self->get_peptide_pvalue($merge_ion_losses{$merged_ion}{'matched_ion_num'}, $merge_ion_losses{$merged_ion}{'Theor_ion_num'}));
								}
							}
							else
							{
								print "please input the right ion scoring methods\n";
							}
							$matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}}=$pep_pvalue;
							
						}
						else
						{
							$pep_pvalue = $matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}};
							$peptide_matched_ratio{$pep} =  $matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}};	
						}
##################### 6/12/2013 Add the modification weight for the peptide with 				
						#$pep_pvalue = $self->weight_p_value_notag($pep,$pep_pvalue,$mod_symbol_used);
						next if (!defined($pep_pvalue));
						next if ($pep_pvalue==0);				
						next if($pep_pvalue=~/nan/);
						next if($pep_pvalue !~ /\d+/);
###########################################################################		

						my $pep_matched_ratio = $self->{'Matched_ion_num'} / $self->{'Theor_ion_num'};

#						$peptide_pvalue{$pep} = $pep_pvalue;
						$peptide_matched_ratio{$pep} =  $pep_pvalue;	
#						$peptide_matched_ratio{$pep} =  $pep_matched_ratio;	
						
						
						

		#			}
	#				print $peptide_matched_ratio{$pep},"\t",$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'},"\t",$matches->{$tag}->{$pep}->{$mod}->{'pepseq'},"\n";
	# version 1.02 $matches->{$tag}->{$pep}->{$mod}->{'rank_p'},
	# version 1.03 $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}}

					my $results = sprintf("%-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-0.2f   %-3s    %-0.4f  %-10s  %-30s    %-20s", 
					$matches->{$tag}->{$pep}->{$mod}->{'theoretical_MH'}/1000,$matches->{$tag}->{$pep}->{$mod}->{'lthreot_Mass'},$matches->{$tag}->{$pep}->{$mod}->{'rthreot_Mass'},
					$peptide_matched_ratio{$pep},$matches->{$tag}->{$pep}->{$mod}->{'tagSeq'},$matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}},$matches->{$tag}->{$pep}->{$mod}->{'tag_rank_num'},$matches->{$tag}->{$pep}->{$mod}->{'sideMass'},$ions{$pep},$matches->{$tag}->{$pep}->{$mod}->{'proteinid'},
					$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'});
					$sort_results{$peptide_matched_ratio{$pep}}{$tag}{$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'}}{$matches->{$tag}->{$pep}->{$mod}->{'proteinid'}}=$results;

#			}
			}
		}
	}


	

	
##################################### Version 12.1.0###################################################
# output the best peptide using the combination of multiple tags
	my $peptide_ref;
	my $peptide_ref_name;
	my $i=0;
	my %PeptideEvaluehash;
	
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{

		foreach my $tag (keys %{$sort_results{$pvalue}})
		{

			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}{$tag}})
			{

				$i++;
				my $j=1;
				foreach my $protein (keys %{$sort_results{$pvalue}{$tag}{$pep}})
				{

					if($j==1)
					{
						my ($M_H, $lSideMass, $rSideMass, $Peptidematchedratio, $TagSeq, $TagEvalue, $TagRankNum, $tagsidemass, $ions ,$Reference, $Peptide) = (split/\s+/,$sort_results{$pvalue}{$tag}{$pep}{$protein})[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11];

						my ($matched,$total)=split(/\//,$ions);
############## speed #########
		

=head
						if(!defined($PeptideEvaluehash{$total}{$matched}))
						{
							$PeptideEvaluehash{$total}{$matched} = sprintf("%.2f",$self->get_peptide_pvalue($matched,$total));
	
							$PeptideEvalue = $PeptideEvaluehash{$total}{$matched};
						}
						else
						{
							$PeptideEvalue = $PeptideEvaluehash{$total}{$matched};						
						}
=cut
						$PeptideEvalue = $pvalue;
						
						
						
						$sort_results{$pvalue}{$tag}{$pep}{$protein} = sprintf("%-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-0.2f    %-2s  %-0.4f  %-10s  %-20s    %-20s", $M_H, $lSideMass, $rSideMass, $PeptideEvalue, $TagSeq, $TagEvalue,  $TagRankNum, $tagsidemass, $ions ,$Reference, $Peptide);

						push @{$peptide_ref->{$pep}->{'info'}}, 
						[$j, $M_H, $lSideMass, $rSideMass, $PeptideEvalue, $TagSeq, $TagEvalue, $TagRankNum];
						$peptide_ref->{$pep}->{'ref'}->{$protein}++;
						
						$peptide_ref_name->{$pep}->{'ref'} = $protein;
						$j++;
					}
					else
					{
						$peptide_ref->{$pep}->{'ref'}->{$protein}++;
						$j++;
					}
				}
			
			}
		}
#		last if ($i>=$parameter->{'number_of_detailed_result'});			
	}
######################### Output the selected candidate peptide ####################

	my $WeightedEvalue = $self->calculate_weighted_score($peptide_ref,$mod_symbol_used);

	
	print OUTPUT "[Selected identified peptide(s)]\n";
    print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq    TagRank  TagNum  Jscore      Reference               Peptide    \n";
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $order = 0;
    for my $Peptide (sort {
            $WeightedEvalue->{$b}->{'Weighted_Evalue'}
            <=>
            $WeightedEvalue->{$a}->{'Weighted_Evalue'} 
			||
            $peptide_ref_name->{$b}->{'ref'}
            cmp
            $peptide_ref_name->{$a}->{'ref'} 			
		} keys %$WeightedEvalue) #TODO: any better sort sub?
    {
		my $Tag_info    = $WeightedEvalue->{$Peptide}->{'longest_Tag_info'};
		
        my $each_info   = join "    ", (@$Tag_info)[1, 2, 3, 4, 5];
        my $Evalue      = $WeightedEvalue->{$Peptide}->{'Weighted_Evalue'};
        my $Num_of_Tags = $WeightedEvalue->{$Peptide}->{'Num_of_Tags'};
		my $tag_rank_num = $WeightedEvalue->{$Peptide}->{'Tag_Rank_Num'};
        my $refs        = $peptide_ref->{$Peptide}->{'ref'};
        my @References  = ();		
		foreach my $protein_id (keys %$refs)
		{
			my @protein_id_array = split(/\|/,$protein_id);
			foreach my $pro_id (@protein_id_array)
			{
				next if($pro_id eq '|');
				my ($protein,$proteinDesc) = $index->getProtein($pro_id,$protein_index_file);
				push (@References,$protein);
			}
		}	

        my $first_ref   = shift @References;
		$order++;

		print OUTPUT $order,"      ";
		print OUTPUT $each_info,"       ";
		printf OUTPUT ("%-2s    %-2s     %-0.2f       %-20s  %-20s\n", $tag_rank_num, $Num_of_Tags,$Evalue,$first_ref,$Peptide);
=head
		if($order==1)
		{
			my $ion_out = $outfile;
			$ion_out=~s/spout/ion/;
			open(IONOUT,">$ion_out") || die "can not open the output file: $output!\n";
			print IONOUT $Peptide,"\n";
			foreach (keys %{$matched_ions{$Peptide}})
			{
				print IONOUT $_,"\t",$matched_ions{$Peptide}{$_},"\n";
			}
			close(IONOUT);
		}
=cut		
#        print OUTPUT "$order  $each_info  $Evalue  $Num_of_Tags  $first_ref   $Peptide\n";
        print OUTPUT "                                                                                            $_\n" for @References;
		last if ($order >= $parameter->{'number_of_selected_result'});		
    }	
	
	print OUTPUT "\n\n";
	print OUTPUT "[All identified peptide(s)]\n";	
	print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  TagEvalue  TagRank  TagSideMass  Ions   Reference                 Peptide   \n";  
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";



######################### Output in detail #####################################	
	my $i=0;
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
		foreach my $tag (keys %{$sort_results{$pvalue}})
		{
			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}{$tag}})
			{
				$i++;
				my $j=1;				
				
				foreach my $protein_id (keys %{$sort_results{$pvalue}{$tag}{$pep}})
				{

		########### 10/31/2014 after reduction of database
	#			my ($proteinName,$proteinDesc) = $index->getProtein($protein_id,$protein_index_file);
					my @protein_id_array = split(/\|/,$protein_id);

					foreach my $pro_id (@protein_id_array)
					{
						next if($pro_id eq '|');
						my ($protein,$proteinDesc) = $index->getProtein($pro_id,$protein_index_file);

					
						if($j==1)
						{
							my @result_array = split(/\s+/,$sort_results{$pvalue}{$tag}{$pep}{$protein_id});
							printf OUTPUT ("%2d    %0.3f     %0.4f    %0.4f      %.2f       %5s      %0.2f      %1s      %.2f  %5s %20s  %20s \n",$i,$result_array[0],$result_array[1],$result_array[2],$result_array[3],$result_array[4],$result_array[5],$result_array[6],$result_array[7],$result_array[8],$protein,$pep);							
							$j++;
						}
						else
						{
							print OUTPUT "                                                                                                  ",$protein,"\n";
							$j++;
						}
					}
					last if ($i>=$parameter->{'number_of_detailed_result'});
				}
			}
		}
	}
	close(OUTPUT);
}

sub WriteResults4NoTags
{
	my ($self)=@_;
	my $database = $self->get_database();
	my $AA_mass = $self->get_AA_mass();	
	my $protein_index_file = $database . ".prdx";	
	my $index = new Spiders::BuildIndex();
	my ($self,$outfile,$matches,$tag_rank_num,$prec_peak_int,$ms2_signal_noise_ratio,$searchtagarray) = @_;
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";
	my $parameter = $self->get_parameter();
	#my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");	
	my %mod_symbol = ('M'=>"@",'S'=>"#",'T'=>"%",'Y'=>"*",'G'=>"^",'K'=>"&",'D'=>'?','A'=>'~','Q'=>'!','P'=>"(",'E'=>")",'E'=>"{",'V'=>"}","V"=>"[","H"=>"]","C"=>":","F"=>",","I"=>';',"L"=>',',"R"=>"<","N"=>">","W"=>"'");	
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";


	print OUTPUT "\nJUMP version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon,$mday,$hour,$min,$sec;	
	my $dtafile = $self->get_dtafile();
	#$dtafile =~ s/tag/dta/;
	print OUTPUT "DTA file = ", $dtafile,"\n";
	
	print OUTPUT "Database = ", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass = ",$self->get_precMass(),"   ", "Percentage of precursor peak intensity = ", sprintf("%0.2f",$prec_peak_int*100),"\%\n";
	print OUTPUT "Precursor matched peptides = ", $self->{'matched_prec_peptide_num'},"\, ","Peptides mass tolerance =", $self->get_mass_tolerance(), "  Peptides mass tolerance units = $parameter->{'peptide_tolerance_units'}\n";
	print OUTPUT "MS2 signal noise ratio = ",sprintf("%0.2f",$ms2_signal_noise_ratio),"\, ","Tag number = ",$self->get_tag_number(),"\n";
	
	print OUTPUT "ion series ABCDVWXYZ: ", $parameter->{'ion_series'},"\n";
	print OUTPUT "ion losses H2O H3PO4 Ammonia: ", $parameter->{'ion_losses_MS2'},"\n";
	my $massutils = new Spiders::MassUtils();
	my $params = $self->get_parameter();
	$massutils->set_parameter($params);	
	
	our ($AA_mass,$dynamic_mod,$fragType,$series,$loss) = $massutils->get_fragment_parameter();	
###### output modification information ###############################
	my $i=0;
	my $mod_symbol_used="";
	foreach (keys %$parameter)
	{
		if($_ =~/add_/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,"=",$parameter->{$_},"  ";
		}

		if($_ =~/dynamic_/ and $parameter->{$_} != 0)
		{
			#print OUTPUT $_,$mod_symbol[$i],"=",$parameter->{$_},"  ";
			if($_ =~/dynamic_(\w+)/ and $parameter->{$_} != 0)
			{
				print OUTPUT $_,$mod_symbol{$1},"=",$parameter->{$_},"  ";
				$i++;
			}
		}
	}

	print OUTPUT "\n\n";	

	
	my $i=0;
	my %sort_results=();
	my %comb_results=();
	my %peptide_pvalue=();
	my %PeptideEvaluehash=();
	my %matched_pep_ions;
	
	foreach my $pep (keys %{$matches})
	{	
		
		my $j=1;
		foreach my $mod (keys %{$matches->{$pep}})
		{	

			my $peptideseq = $matches->{$pep}->{$mod}->{'matchPeptide'};

			if(!defined($peptide_pvalue{$pep}))
			{	

				my ($matched,$match_ion) = $self->get_peptide_matched($peptideseq,$mod,$matches->{$pep}->{$mod}->{'theoretical_MH'}/1000,$AA_mass);
				next if ($matched<3);

				my $pep_pvalue = 0;

				$ions{$pep} = "$self->{'Matched_ion_num'}\/$self->{'Theor_ion_num'}";

						
				if(!defined($matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}}))
				{
					$matched_ions{$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'}} = $match_ion;
					next if ($matched==0);
													
	############### seperately scoring for different ion types	##################					
					my $ion_scoring = 1;
					$ion_scoring = $parameter->{'ion_scoring'};

					if($ion_scoring == 1)
					{
						$pep_pvalue = sprintf("%0.2f",$self->get_peptide_pvalue($self->{'Matched_ion_num'}, $self->{'Theor_ion_num'}));
					}
					elsif($ion_scoring == 2)
					{
						my %merge_ion_losses;
						foreach my $ion_type (keys %$match_ion)
						{
							my ($ion,$loss)=split(/-/,$ion_type);
							my ($matched_num,$theor_num) = split(/\|/,$match_ion->{$ion_type}); 
							$merge_ion_losses{$ion}{'matched_ion_num'} += $matched_num;
							$merge_ion_losses{$ion}{'Theor_ion_num'} += $theor_num;								
						}
						foreach my $merged_ion (keys %merge_ion_losses)
						{
							next if($merge_ion_losses{$merged_ion}{'matched_ion_num'}==0);
						
							$pep_pvalue += sprintf("%0.2f",$self->get_peptide_pvalue($merge_ion_losses{$merged_ion}{'matched_ion_num'}, $merge_ion_losses{$merged_ion}{'Theor_ion_num'}));
						}
					}
					else
					{
						print "please input the right ion scoring methods\n";
					}
					$matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}}=$pep_pvalue;
							
				}
				else
				{
					$pep_pvalue = $matched_pep_ions{$self->{'Theor_ion_num'}}{$self->{'Matched_ion_num'}};
				}

=head				
				if(!defined($PeptideEvaluehash{$total}{$matched}))
				{	
					$pep_pvalue = sprintf("%.2f",$self->get_peptide_pvalue($self->{'Matched_ion_num'} , $self->{'Theor_ion_num'}));

					$PeptideEvaluehash{$total}{$matched} = $pep_pvalue;
				}
				else
				{
					$pep_pvalue=$PeptideEvaluehash{$total}{$matched};
				}
=cut				
##################### 6/12/2013 Add the modification weight for the peptide with 				
				#$pep_pvalue = $self->weight_p_value_notag($pep,$pep_pvalue,$mod_symbol_used);
				next if (!defined($pep_pvalue));
				next if ($pep_pvalue==0);				
				next if($pep_pvalue=~/nan/);
				next if($pep_pvalue !~ /\d+/);
###########################################################################					
				$peptide_pvalue{$pep} = $pep_pvalue;
			
#				print $peptideseq,"\t",$mod,"\t",$matched,"\t",$total,"\t",$self->{'Matched_ion_num'},"\t",$self->{'Theor_ion_num'},"\t",$pep_pvalue,"\t",$peptide_pvalue{$pep},"\n";
			}	
			my $ions = $self->{'Matched_ion_num'} . "/$self->{'Theor_ion_num'}";
			
			my $results = sprintf("%-0.3f    %-0.4f    %-0.4f       %-0.2f         %-6s     %-0.2f   %-2d   %-0.4f       %-6s     %-20s    %-20s", 
			$matches->{$pep}->{$mod}->{'theoretical_MH'}/1000,$matches->{$pep}->{$mod}->{'lthreot_Mass'},$matches->{$pep}->{$mod}->{'rthreot_Mass'},
			$peptide_pvalue{$pep},$matches->{$pep}->{$mod}->{'tagSeq'},$matches->{$pep}->{$mod}->{$parameter->{'tag_select_method'}},0,$matches->{$pep}->{$mod}->{'sideMass'},$ions,$matches->{$pep}->{$mod}->{'proteinid'},
				$matches->{$pep}->{$mod}->{'matchPeptide_orig'});
#			print $peptide_pvalue{$pep},"aa\t",$tag,"bb\t",$matches->{$pep}->{$mod}->{'matchPeptide_orig'},"\t",	$matches->{$pep}->{$mod}->{'proteinid'},"\n";
#			$sort_results{$peptide_pvalue{$pep}}{$tag}{$matches->{$pep}->{$mod}->{'matchPeptide_orig'}}{$matches->{$pep}->{$mod}->{'proteinid'}}=$results;
			$sort_results{$peptide_pvalue{$pep}}{$matches->{$pep}->{$mod}->{'matchPeptide_orig'}}{$matches->{$pep}->{$mod}->{'proteinid'}}=$results;			
		}
	}


	print OUTPUT "[Selected identified peptide(s)]\n";
    print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  TagRank   TagNum   Jscore     Reference                       Peptide    \n";
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $order = 0;
	my $i=0;
	
	LOOP:
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
#		foreach my $tag (keys %{$sort_results{$pvalue}})
#		{
			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}})
			{
				$i++;
				my $j=1;
				foreach my $protein_id (keys %{$sort_results{$pvalue}{$pep}})
				{	
				#	my $protein_id = $sort_results{$pvalue}{$pep}{$protein_key};
				
	########### 10/31/2014 after reduction of database
#			my ($proteinName,$proteinDesc) = $index->getProtein($protein_id,$protein_index_file);
					my @protein_id_array = split(/\|/,$protein_id);
				#	print "hhhh",$protein_id,"ffff\t",@protein_id_array,"dddddd\n";
					foreach my $pro_id (@protein_id_array)
					{
						next if($pro_id eq '|');
						my ($protein,$proteinDesc) = $index->getProtein($pro_id,$protein_index_file);


				
#				foreach my $protein (keys %{$sort_results{$pvalue}{$pep}})
#				{
						if($j==1)
						{
			############## check whether we can find other tags matched to this first peptide

							my @result_array = split(/\s+/,$sort_results{$pvalue}{$pep}{$protein_id});		
	#						if($i==1)
	#						{
								my $flag = 0;
								for(my $l=0;$l<=$#$searchtagarray;$l++)
								{
									my $tags =  $searchtagarray->[$l];
									my $pepseq = $pep;


									my $mod = $self->get_mod($pepseq);
									my $frag_tolerance = $self->get_frag_tolerance();
									my @seq_array = split(/\./,$pepseq);
									my $pepseq_core = $seq_array[1];		
									my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$pepseq_core,$mod,$tags,$frag_tolerance);
									if ($sideMassMatchResult->{'lSide_match'})
									{
										printf OUTPUT ("%2d    %0.3f     %0.4f    %0.4f      %.2f   %5s    %1s    %1s    %.2f   %20s  %20s \n",$i,$result_array[0],$tags->{'sideMass'},$result_array[0]-$tags->{'sideMass'},$result_array[3],$tags->{'tagSeq'},$l,1,$result_array[3],$protein,$pep);							
										$flag = 1;		
										$l=$#$searchtagarray+1;
									}
									elsif($sideMassMatchResult->{'rSide_match'})
									{
										printf OUTPUT ("%2d    %0.3f     %0.4f    %0.4f      %.2f   %5s    %1s    %1s    %.2f   %20s  %20s \n",$i,$result_array[0],$result_array[0]-$tags->{'sideMass'},$tags->{'sideMass'},$result_array[3],$tags->{'tagSeq'},$l,1,$result_array[3],$protein,$pep);							
										$flag = 1;		
										$l=$#$searchtagarray+1;									
									}
								}
								if($flag == 0)
								{
									printf OUTPUT ("%2d    %0.3f     %0.4f    %0.4f      %.2f   %5s    %1s    %1s    %.2f   %20s  %20s \n",$i,$result_array[0],$result_array[1],$result_array[2],$result_array[3],$result_array[4],0,0,$result_array[3],$protein,$pep);							
								}							
								
	#						}

					#		printf OUTPUT ("%2d    %10s\n",$i,$sort_results{$pvalue}{$tag}{$pep}{$protein});
							$j++;
						}
						else
						{
							print OUTPUT "                                                                                  ",$protein,"\n";
							$j++;
						}
					}
					last LOOP if ($i>=$parameter->{'number_of_selected_result'});
				}
			}
#		}
	}

	print OUTPUT "\n\n";

	
	print OUTPUT "[All identified peptide(s)]\n";	
	print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  TagEvalue  TagRank  TagSideMass  Ions   Reference                 Peptide   \n";  
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";



######################### Output in detail #####################################	
	my $i=0;
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
		foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}})
		{
			$i++;
			my $j=1;
			foreach my $protein_id (keys %{$sort_results{$pvalue}{$pep}})
			{
	########### 10/31/2014 after reduction of database
#			my ($proteinName,$proteinDesc) = $index->getProtein($protein_id,$protein_index_file);
				my @protein_id_array = split(/\|/,$protein_id);
				#	print "hhhh",$protein_id,"ffff\t",@protein_id_array,"dddddd\n";
				foreach my $pro_id (@protein_id_array)
				{
					next if($pro_id eq '|');
					my ($protein,$proteinDesc) = $index->getProtein($pro_id,$protein_index_file);

				
					if($j==1)
					{
						my @result_array = split(/\s+/,$sort_results{$pvalue}{$pep}{$protein_id});
						printf OUTPUT ("%2d    %0.3f     %0.4f    %0.4f      %.2f       %5s      %0.2f      %1s      %.2f  %5s %20s  %20s \n",$i,$result_array[0],$result_array[1],$result_array[2],$result_array[3],$result_array[4],$result_array[5],$result_array[6],$result_array[7],$result_array[8],$protein,$pep);							
						$j++;
					}
					else
					{
						print OUTPUT "                                                                                                  ",$protein,"\n";
						$j++;
					}
				}
				last if ($i>=$parameter->{'number_of_detailed_result'});
			}
		}
	}

	close(OUTPUT);
}


sub calculate_weighted_score # private method
{   
	my ($self, $peptide_ref,$mod_symbol_used) = @_;
	my $parameter = $self->get_parameter();	
	my $charge = $self->get_precCharge();
	
	my $coeff = $parameter->{'tag_coeff_evalue'};
	my $method=3;
    my $WE;
    ((caller)[0] eq ref $self) or die "Private method called.\n";
    for my $Peptide (keys %$peptide_ref)
    {   
		my $info        = $peptide_ref->{$Peptide}->{'info'};
        #my $Reference   = $peptide_ref->{'ref'};
        #--------
        $WE->{$Peptide}->{'Num_of_Tags'} = scalar @$info;
        my $PeptideEvalue       = $info->[0]->[4]; # all are equal for one peptide
		
        my @sorted_lines        = sort { length $b->[5] <=> length $a->[5] } @$info;
        my $longest_tag_line    = shift @sorted_lines;
        $WE->{$Peptide}->{'longest_Tag_info'} = $longest_tag_line;
        my ($longest_Tag_seq, $longest_Tag_Evalue) = (@$longest_tag_line)[5, 6];
######## Get the highest score  
		@sorted_lines        = sort { $b->[6] <=> $a->[6] } @$info;
        $highestscore_tag_line    = shift @sorted_lines;
        my ($highest_Tag_seq,$highest_Tag_Evalue,$highest_Tag_Rank) = (@$highestscore_tag_line)[5, 6,7];
#####################################################		
		
#        $WE->{$Peptide}->{'Weighted_Evalue'} = $PeptideEvalue + $highest_Tag_Evalue;
        $WE->{$Peptide}->{'Tag_Evalue'} = $highest_Tag_Evalue;
		$WE->{$Peptide}->{'Tag_Rank_Num'} = $highest_Tag_Rank;
        my @TagSeq_Evalue       = map { [$_->[5], $_->[6]] } @sorted_lines;
	#	print $PeptideEvalue,"\t",$longest_Tag_Evalue,"\t",$highest_Tag_Evalue,"\n";
        for my $each_pair (@TagSeq_Evalue)
        {   
			my ($TagSeq, $TagEvalue) = @$each_pair;
        #    $WE->{$Peptide}->{'Weighted_Evalue'} += ($TagSeq =~ /$longest_Tag_seq/) ?  $TagEvalue : ($coeff  * $TagEvalue);
			if($method==1)
			{
				next if ($TagEvalue == $highest_Tag_Evalue);
				
				
				next if (length($highest_Tag_seq)>length($TagSeq) and  $highest_Tag_seq=~/$TagSeq/);
				my @substr = $self->lc_substr($highest_Tag_seq,$TagSeq);
				if(scalar(@substr)==0)
				{
					$WE->{$Peptide}->{'Tag_Evalue'} += $TagEvalue;
					#print $highest_Tag_seq,"aaa\t",$TagSeq,"bbbbb\n";
				}
				else
				{			
					#$WE->{$Peptide}->{'Tag_Evalue'} += ($TagEvalue == $highest_Tag_Evalue) ?  0 : ($coeff  * $TagEvalue);
				}
			}
			elsif($method==2)
			{
########## compare the other tag sequence with tag with highest score, if the latter sequence already covers the tag sequences, then ignore
####################################################################### if the new tag is longer than the tag with highest score, then calculate a coefficient for the tag score 

				next if ($TagEvalue == $highest_Tag_Evalue);

	
				next if (length($highest_Tag_seq)>length($TagSeq) and  $highest_Tag_seq=~/$TagSeq/);
############# If the highest_Tag_seq is derived from y ion, there are some tags derived from b ions
				my $rev_tag_seq = reverse $TagSeq;
				$rank_tag =~ s/([^\w])([\w])/$2$1/g;
				next if (length($highest_Tag_seq)>length($TagSeq) and  $highest_Tag_seq=~/$rev_tag_seq/);
				

				my @substr = $self->lc_substr($highest_Tag_seq,$TagSeq);
				if(scalar(@substr)==0)
				{
					$WE->{$Peptide}->{'Tag_Evalue'} += 0.1 *  $TagEvalue; 					
				}
				else
				{
					$WE->{$Peptide}->{'Tag_Evalue'} += 0.1 * (length($TagSeq) - scalar(@substr)) / length($highest_Tag_seq) * $TagEvalue;
				}

			}
			elsif($method==3)
			{
			
			}

        }
		$WE->{$Peptide}->{'Weighted_Evalue'} = $PeptideEvalue + $coeff * $WE->{$Peptide}->{'Tag_Evalue'};
		
####### 6/11/2013   Add weighted score of H3PO4-neutral loss 
			
    }
    return $WE;
}

sub weight_p_value_notag
{
	my ($self, $pep,$pep_pvalue,$mod_symbol_used) = @_;
	my $parameter = $self->get_parameter();	
	
	if($parameter->{'pho_neutral_loss'})
	{
		if($self->get_pho_neutral_loss())
		{
			if($pep=~/$mod_symbol_used/)
			{
				$pep_pvalue = $pep_pvalue + $parameter->{'pho_neutral_loss'} * $pep_pvalue;
			}
		}
	}
	return $pep_pvalue;
}

sub lc_substr {
  my ($self, $str1, $str2) = @_; 
  my $l_length = 0; # length of longest common substring
  my $len1 = length $str1; 
  my $len2 = length $str2; 
  my @char1 = (undef, split(//, $str1)); # $str1 as array of chars, indexed from 1
  my @char2 = (undef, split(//, $str2)); # $str2 as array of chars, indexed from 1
  my @lc_suffix; # "longest common suffix" table
  my @substrings; # list of common substrings of length $l_length
 
  for my $n1 ( 1 .. $len1 ) { 
    for my $n2 ( 1 .. $len2 ) { 
      if ($char1[$n1] eq $char2[$n2]) {
        # We have found a matching character. Is this the first matching character, or a
        # continuation of previous matching characters? If the former, then the length of
        # the previous matching portion is undefined; set to zero.
        $lc_suffix[$n1-1][$n2-1] ||= 0;
        # In either case, declare the match to be one character longer than the match of
        # characters preceding this character.
        $lc_suffix[$n1][$n2] = $lc_suffix[$n1-1][$n2-1] + 1;
        # If the resulting substring is longer than our previously recorded max length ...
        if ($lc_suffix[$n1][$n2] > $l_length) {
          # ... we record its length as our new max length ...
          $l_length = $lc_suffix[$n1][$n2];
          # ... and clear our result list of shorter substrings.
          @substrings = ();
        }
        # If this substring is equal to our longest ...
        if ($lc_suffix[$n1][$n2] == $l_length) {
          # ... add it to our list of solutions.
          push @substrings, substr($str1, ($n1-$l_length), $l_length);
        }
      }
    }
  }   
 
  return @substrings;
}



sub exportMatches
{
	my ($self, $f, $matches) = @_;
#### changed by xusheng ########

	open OUT, '>'.$f or die;
	foreach my $tid (sort keys %$matches)
	{
		foreach my $pep (keys %{$matches->{$tid}})
		{

			next if ($matches->{$tid}->{$pep}->{'tagSeqMatch'}==0);
			print OUT $matches->{$tid}->{$pep}->{'scanNum'},"\t",$matches->{$tid}->{$pep}->{'precMass'},"\t",$matches->{$tid}->{$pep}->{'threotical_MH'},"\t",$matches->{$tid}->{$pep}->{'tagSeq'},"\t",$matches->{$tid}->{$pep}->{'matchPeptide'}->{'seq'},"\t",$pep,"\t"; #$matches->{$tid}->{$pep}->{'matchPeptide'}->{'desc'},"\t";
			if($matches->{$tid}->{$pep}->{'tagSeqMatch'}==1)
			{
				print OUT "Tag sequence matched\t";
			}
			elsif($matches->{$tid}->{$pep}->{'tagSeqMatch'}>1)
			{
				print OUT "Partial tag matched\t";
			}
			else
			{
				print OUT "Tag sequence unmatched\t";
			}
			print OUT $matches->{$tid}->{$pep}->{'sideMass'},"\t";
			if($matches->{$tid}->{$pep}->{'lthreot_Mass'})
			{
				print OUT $matches->{$tid}->{$pep}->{'lthreot_Mass'},"\t";
			}
			else
			{
				print OUT "\t";
			}
				
			if($matches->{$tid}->{$pep}->{'lsideMassMatch'}>0 or $matches->{$tid}->{$pep}->{'rsideMassMatch'}>0)
			{
				print OUT "lSide mass matched\t";	
			}
			else
			{
				print OUT "lSide mass unmatched\t";
			}
			
			if($matches->{$tid}->{$pep}->{'rthreot_Mass'})
			{
				print OUT $matches->{$tid}->{$pep}->{'rthreot_Mass'},"\t";
			}
			else
			{
				print OUT "\t";
			}
			if($matches->{$tid}->{$pep}->{'rsideMassMatch'}>0)
			{
				print OUT "rSide mass matched\t";
			}
			else
			{
			
				print OUT "rSide mass unmatched\t";
			}
			
			print OUT $matches->{$tid}->{$pep}->{'rankP'},"\t",$matches->{$tid}->{$pep}->{'hyperP'},"\n";
		}
	}
	close OUT;
}	


sub get_isotopic_distribution
{
	my ($self, $mass) = @_;
	my $loop = 0;
	if($mass<1500)
	{
		$loop = 1;
	}
	elsif($mass<3000 && $loop<=2)
	{
		$loop = 2;
	}
	elsif($mass<4500 && $loop<=4)
	{
		$loop = 4;
	}
	elsif($mass<6000 && $loop<=6)
	{
		$loop = 6;
	}	
	return $loop;	
	
}

sub get_mod
{
	my ($self,$seq) = @_;
	my $mod;

	my %mod_pos=();
					
	my $nomod_pep = $seq;

	my @seq_array = split(/\./,$nomod_pep);
	$nomod_pep = $seq_array[1];

	while($nomod_pep =~ /[^a-zA-Z]+/g)
	{
		$mod_pos{$-[0]}=1;
		$nomod_pep =~ s/[^a-zA-Z]+//;
	}				
	if((scalar keys %mod_pos)==0)
	{
	############### 8/16/2013 fixed a length issue 				
		my $length=length($seq);
		$length=length($nomod_pep);
						
		$mod = ":" x ($length+1);
	}
	else
	{
	############## if the peptide has modification, get the modification info ##################	
	# Fix a bug: K.QIVWKYCGR.M@			
	#					my $nomod_pep = $peptidehash->{$PeptideID}->{'seq'};
	#					my @seq_array = split(/\./,$nomod_pep);
	#					$nomod_pep = $seq_array[1];			
		my @nomodseq = split(//,$nomod_pep);
					
		my $length=length($nomod_pep);

		for(my $i=0;$i<=$length;$i++)
		{
			if($mod_pos{$i+1})
			{
				$mod .= ":";
				$mod .= $nomodseq[$i];
			}
			else
			{
				$mod .= ":";
			}
		}
	}
	return $mod;
	
}


1;

__END__
