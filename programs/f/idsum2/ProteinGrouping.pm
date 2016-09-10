#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: idsum2::ProteinGrouping

package idsum2::ProteinGrouping;

sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);

  return $self;
}

# function to get all the proteins with shared peptides into a group, or shared protein hash 
# %sharedpros{$shared_protein}=1;
sub get_shared_proteins{

	# read references of shared protein hash, protein hash, peptide hash and target protein
	my ($sharedpros, $protein_hash, $peptide_hash, $protein) = @_;
	
	# find peptides of the target protein from the protein hash
	# for each peptide
	for my $peptide (keys %{$$protein_hash{$protein}{'peptides'}})
	{
		# if the peptide is unique, then go to the next peptide without further going down
		next if ($$peptide_hash{$peptide}{'unique'} == 1);

		# find all the proteins for the peptide from the peptide hash
		# for each protein
  		for my $shared (keys %{$$peptide_hash{$peptide}{'proteins'}})
		{
		
			next if ($shared=~/Decoy/);
			
			# if the target protein from the function input is not defined in the protein hash, exit the program
			# this part is for debug purpose only, not related to the function
			if (!defined($$protein_hash{$protein}))
			{
				print Dumper($$peptide_hash{$peptide});
				exit;
				delete $$peptide_hash{$peptide}{'proteins'}{$protein};
				next;
			}
			
			# if the protein is already defined in the shared protein hash, go to the next protein
    		next if (defined($$sharedpros{$shared}));
			# add the protein in the shared protein hash
   			$$sharedpros{$shared} = 1;
			# recursively call the current function to search the downstream shared proteins for the current protein or current shared protein
   	 		get_shared_proteins($sharedpros, $protein_hash, $peptide_hash, $shared);
  		}
	}
}

# function to group proteins
sub group_proteins{
	shift @_;
	# take the input variables, protein hash reference, peptide hash reference, subgroup percent, parameter hash reference
	my ($protein_hash, $peptide_hash, $param_hash) = @_;
	# take bypass_grouping parameter value
	my $bypass = $$param_hash{'bypass_grouping'};
	# take bypass filtering parameter value
	my $bypassfil = $$param_hash{'bypass_filtering'};
	# hash for proteins that have been used for grouping
	my %used;
	# define group number variable, starting from 1
	my $group_num = 1;
	# define subgroup number variable, starting from 1
	my $subgroup_num = 1;
	# if no protein defined in protein hash, report error
		
	my $order = 1;
	
	my $grouphash;
		
	if (scalar(keys %$protein_hash) == 0)
	{
		print "No proteins were identified for this group !!!\n";
		return (0,0);
	}
	# sort the protein in protein hash by spectral count, total peptides, and unique peptides
	# for each protein after sorting

	#print "\nGrouping proteins ...\n";
	for my $protein (sort {	$$protein_hash{$b}{'occurrence'} <=> $$protein_hash{$a}{'occurrence'} || 
							$$protein_hash{$b}{'total_nomod'} <=> $$protein_hash{$a}{'total_nomod'} ||
							$$protein_hash{$b}{'unique_nomod'} <=> $$protein_hash{$a}{'unique_nomod'} ||
							$$protein_hash{$b}{'max_xcorr'} <=> $$protein_hash{$a}{'max_xcorr'}
						   } keys %$protein_hash)
	{

		next if ($protein=~/Decoy/); # include decoy proteins for protein grouping
		# if the protein is already defined in the used protein hash, then go to the next protein
		next if (defined($used{$protein}));
#		print $protein,"\t",$$protein_hash{$protein}{'occurrence'},"\t",$$protein_hash{$protein}{'total_nomod'},"\t",$$protein_hash{$protein}{'unique_nomod'},"\n\n";		
		# if the protein doesn't have any peptide, then delete the protein, go to the next protein
		if (scalar (keys %{$$protein_hash{$protein}{'peptides'}}) == 0)
		{
			delete ($$protein_hash{$protein});
			next;
		}
	#	print "\rGrouping Proteins $group_num: $protein                     ";
		# if the protein doesn't have any shared peptide
		if ($$protein_hash{$protein}{'shared'} == 0)
		{
			$$protein_hash{$protein}{'group'} = $group_num;
			$$protein_hash{$protein}{'subgroup'} = 1;
			$used{$protein} = $group_num;			
			$grouphash->{$protein}->{'group'} = $group_num;			
			$grouphash->{$protein}->{'subgroup'} = 1;			
			$grouphash->{$protein}->{'order'} = 1;
			$subgroup_num += 1;	
		} 
		else 
		{
			my %shared; 
			$shared{$protein} = 1;
			get_shared_proteins(\%shared, $protein_hash, $peptide_hash, $protein);

			($subgrouphash) = subgroup_proteins(\%shared, $protein_hash);
		
################# if subgroup is 1 and then classify the group by name with "uncharacterize"			
			$subgrouphash = classify_protein($subgrouphash,\%shared,$protein_hash);
			
		#		print Dumper($$protein_hash{$protein}) if ($debug);
			if ($bypass == 1)
			{
				for my $pro (keys %shared)
				{
					$$protein_hash{$pro}{'group'} = $group_num;
					$used{$pro} = $group_num;
				}

			}
			
			foreach $protein (keys %$subgrouphash)
			{


				$grouphash->{$protein}->{'group'} = $group_num;
				$grouphash->{$protein}->{'subgroup'} = $subgrouphash->{$protein}->{'subgroup'};			
				$grouphash->{$protein}->{'order'} = $subgrouphash->{$protein}->{'order'};
				
			}
			
		}
		
		# add 1 on protein group number
		$group_num++;
	}
	
	#print "\n";
	return ($grouphash,$group_num,$subgroup_num);
}

# build %subgrouphash: assign subgroup and order numbers
# $subgrouphash->{$protein}->{'subgroup'} = $subgroup;
# $subgrouphash->{$protein}->{'order'} = order the proteins (which can not be distinguished from each other) within one subgroup
sub subgroup_proteins{
	my ($group_hash, $protein_hash) = @_;
	my $subgroup = 0;
	my $order = 1;
	my $subgrouphash;

	
	
	for my $protein (sort {$$protein_hash{$b}{'occurrence'} <=> $$protein_hash{$a}{'occurrence'} ||  
					$$protein_hash{$b}{'total_nomod'} <=> $$protein_hash{$a}{'total_nomod'} ||
                    $$protein_hash{$b}{'unique_nomod'} <=> $$protein_hash{$a}{'unique_nomod'} ||
					$$protein_hash{$b}{'max_xcorr'} <=> $$protein_hash{$a}{'max_xcorr'}
                    } keys %$group_hash)
	{
		if($subgroup==0)
		{
			$subgroup = 1;
			$subgrouphash->{$protein}->{'subgroup'} = $subgroup;
			$subgrouphash->{$protein}->{'order'} = 1;			
		}
		elsif($$protein_hash{$protein}{'unique_nomod'}>0)
		{
			$subgroup++;
			$subgrouphash->{$protein}->{'subgroup'} = $subgroup;
			$subgrouphash->{$protein}->{'order'} = 1;
		}
		else
		{
			$order++;
			$subgrouphash->{$protein}->{'subgroup'} = 1;
			$subgrouphash->{$protein}->{'order'} = $order;			
		}
	}
	return ($subgrouphash);
}

# ??
sub classify_protein
{
	my ($subgrouphash,$sharedhash,$protein_hash) = @_;
	my $order = 0;

	for my $protein (sort {$protein_hash->{$a}->{'Unchar'} <=> $protein_hash->{$b}->{'Unchar'} ||
					$protein_hash->{$b}->{'occurrence'} <=> $protein_hash->{$a}->{'occurrence'} ||  
					$protein_hash->{$b}->{'total_nomod'} <=> $protein_hash->{$a}->{'total_nomod'} ||
                    $protein_hash->{$b}->{'unique_nomod'} <=> $protein_hash->{$a}->{'unique_nomod'} ||
					$protein_hash->{$b}->{'max_xcorr'} <=> $protein_hash->{$a}->{'max_xcorr'}
                    } keys %$subgrouphash)				
	{
		if($subgrouphash->{$protein}->{'subgroup'}==1)
		{
			$order++;
			$subgrouphash->{$protein}->{'order'}=$order;
		}
	}
	return $subgrouphash;
}


1;

