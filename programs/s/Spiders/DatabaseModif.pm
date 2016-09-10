#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::DatabaseModif

######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Coverter MzXML to dta files	          ****    #
#       ****					                      ****    #
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #
#       ****all rights reserved.		              ****    #
#       ****xusheng.wang@stjude.org		              ****    #
#       ****					                      ****    #
#       ****					                      ****    #
#       **************************************************    #
###############################################################

###### The package of task1 has been renamed because there is same fuction in Algorithm::Combinatorics. 
###### USAGE: my %result = changeValues(\%Motif, $comb_num, \%massHash);
#!/usr/bin/perl -I /home/xwang4/scripts
########################### Example #################################
#use strict;
#use warnings;
#use Spiders::DatabaseModif;

#my %modif = ( 'M' => 57, 'T' => 78, 'S' => 79 );
#my %new_modif = ('M'=>'M@','T'=>'T#','S'=>'S$');
#my $comb_num = 6;

#my %pephash = (
#    '1050.25124' => ['PITWGSDVAR'],
#);

#my $cmdatabase = new Spiders::DatabaseModif();
#$cmdatabase->GenerateModifSeq(\%modif,$comb_num,\%pephash,\%new_modif);

#foreach my $key ( keys(%pephash) ) {
#    print "$key=>";
#    foreach ( @{ $pephash{$key} } ) {
#        print "$_,";
#    }
#    print "\n";
#}
#######################################################################

package Spiders::DatabaseModif;
      
use strict;
use warnings;
use Exporter;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();


sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}

### This subroutine is used to change the mass value for all possible modifications for database########
### 
sub GenerateModifSeq
{
	my ($self,$modif,$params,$pephash,$new_modif)=@_;
####### generate combination of modified amino acid 
	my $combNUM = $params->{'max_modif_num'}-1;

    my %mid_hash =$self->modify_combinations(\%$modif, $combNUM+1);
############ For each peptide in the peptide hash, change the mass value #######		
    foreach my $pepmass (keys(%$pephash))
    {
        foreach my $peptide (@{$$pephash{$pepmass}})
        {
			print $peptide,"\n";

			my $flag        = 0;
		
			my @modif_array = ();
########## The following codes might not efficient enough for large database 			
			my @peptide_aa = split(//, $peptide);
			my $i=0;
			my %modif_pos_aa;
			foreach my $aa (@peptide_aa)
			{
				$i++;
				if (exists($$modif{$aa}))
				{
					$flag++;
############# save the position, rather than amino acid itself ################					
					$modif_pos_aa{$i}=$aa;			
				}
			
			}
			if($flag>10)
			{
				$flag = 2;
			}
			if($flag>($combNUM+1))
			{
				$flag = $combNUM + 1;
			}
######### @modif_array is used to save amino acid 			
			@modif_array = map {$_} sort {$a<=>$b} keys %modif_pos_aa;	
####### If there is more than one modification in the peptide ############
			if ($flag != 0)
			{
				for (my $i = 1 ; $i <= $flag ; $i++)
				{
	#				print $i,"bbbbbbbbbbbb\n";
################## do not pass the @array as a parameter (??????????) ###########################					
					my $iter = $self->comb($i,\@modif_array,$combNUM);

					foreach my $c (@$iter)
					{

						
						my $changed_peptide = $peptide;
						
						my @pos_aa = map {$modif_pos_aa{$_}} @$c;
	#					print @pos_aa,"\n";
						my $join = join("", sort(@pos_aa));	
						my $modif_pepmass = 0;

						if (exists($mid_hash{$join}))
						{
				############### mass = modif_mass * 1000		
							$modif_pepmass = $pepmass + $mid_hash{$join} * 1000;
#							print $modif_pepmass,"\t",$params->{'max_peptide_mass'},"\n";
############## Version 12.1.0 6/3/2013########################### 						
# excluding those peptides with modification > max mass
# otherwise, Eat a lot of memory
							next if ($modif_pepmass > ($params->{'max_peptide_mass'} * 1000));	
						}	
						
						my $j=0;
						foreach my $pos (@$c)
						{
########replace the amino acid with modified amino acid ##############################
								substr($changed_peptide,$pos-1+$j,1,$new_modif->{$modif_pos_aa{$pos}});
								$j++;
						}
#						print $changed_peptide,"\t",$modif_pepmass,"\t",$join,"\n";					
############ @pos_aa is used to save amino acid, rather than position ##################	

						if (exists($mid_hash{$join}))
						{					
					#		print $modif_pepmass,"\t",$changed_peptide,"\n";
							
							push (@{$pephash->{$modif_pepmass}}, $changed_peptide);
						}							
					}
					undef @$iter;

				}
			}
		}
    }
}

sub modify_combinations
{
	my ($self,$hash,$comb_num) = @_;
    my %new_hash;

    for (my $i = 0 ; $i < $comb_num ; $i++)
    {
        if ($i == 0)
        {
            foreach my $key (keys(%$hash))
            {
                $new_hash{$key} = $$hash{$key};
            }
        }
        else
        {
            foreach my $hash_key (keys(%new_hash))
            {
                foreach my $key (keys(%$hash))
                {
                    my $key_join = "$hash_key" . "$key";
                    my @alpha    = sort(split(//, $key_join));
                    my $key_comb = join('', @alpha);

                    $new_hash{$key_comb} = $$hash{$key}  + $new_hash{$hash_key};

                }
            }
        }
    }
    return (%new_hash);
}

sub comb{
	my ($self, $n, $dataref,$combNUM)=@_;
	my @data=@$dataref;
	my @result;
	undef @result;
#	foreach (@data)
#	{
#		print $_," ";
#	}
#	print $n,"\n";
	return [map {[$_]} @data] if $n==1;
	while(1)
	{
		last if @data<$n;
		my $item=shift @data;
		my $ret=$self->comb($n-1,\@data,$combNUM);
		for(@$ret)
		{
	#	 unshift @$_,$item;
####### added on 6/3/2013 ##########
## to avoid eating a lot of memory for those peptide with many aa that can be modified		 
		 next if (scalar(@$_)>$combNUM);
		 unshift @$_,$item;
##########################################		 
		 push @result,$_;
		}

		undef $ret;
	}
#	print $n,"\t",$#result,"bbbbbb\n";
	return \@result;
}

1;

