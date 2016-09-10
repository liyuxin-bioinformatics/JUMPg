#!/usr/bin/perl


######### Deisotope ##########################################
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
# this value is used to remove the precursor ion

############################################################### 
# getFragmentMasses subroutine is obtained from  InSilicoSpectro 
# which is developed by Geneva Bioinformatics www.genebio.com
##############################################################

package Spiders::MassUtils;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);
use Spiders::MathUtils;
use Spiders::Params;

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();


sub new{
	my ($class,%arg)=@_;
    my $self = {
		_C_value =>undef,
		_H_value=>undef,
		_amino_acid_mass=>undef,
    };
    bless $self, $class;
	return $self;
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
		$self->{_C_value}=1.00335;
	}
	return $self->{_C_value};
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


sub set_amino_acid_mass
{
	my ($self,$AA_mass);
	return $self->{'_amino_acid_mass'} = $AA_mass;
}

############ static modification was added in the mass of AA#############
sub get_AA_mass
{
	my ($self) = shift;

	my $AA_mass;
	my $parameter = $self->get_parameter();
	if(!defined($self->{'_amino_acid_mass'}))
	{	
		$AA_mass->{'G'} = 57.02146372057;	
		$AA_mass->{'D'} = 115.02694302383;
		$AA_mass->{'A'} = 71.03711378471;
		$AA_mass->{'Q'} = 128.05857750528;
		$AA_mass->{'S'} = 87.03202840427;
		$AA_mass->{'K'} = 128.094963014;
		$AA_mass->{'P'} = 97.05276384885;
		$AA_mass->{'E'} = 129.04259308797;
		$AA_mass->{'V'} = 99.06841391299;
		$AA_mass->{'M'} = 131.04048491299;
		$AA_mass->{'T'} = 101.04767846841;
		$AA_mass->{'H'} = 137.05891185845;
		$AA_mass->{'C'} = 103.00918478471;
		$AA_mass->{'F'} = 147.06841391299;
		$AA_mass->{'I'} = 113.08406397713;		
		$AA_mass->{'L'} = 113.08406397713;
		$AA_mass->{'J'} = 113.08406397713;
		$AA_mass->{'R'} = 156.1011110236;
		$AA_mass->{'N'} = 114.04292744114;	
		$AA_mass->{'Y'} = 163.06332853255;
		$AA_mass->{'W'} = 186.07931294986;
	}
### add the static modification ################## 		
	foreach my $param (keys %$parameter)
	{
		if($param =~/add\_(\w+)\_/ and $parameter->{$param}>0)
		{
			$AA_mass->{$1} += $parameter->{$param};
		}
	}

######### Get modification 	
	return $AA_mass;
}


sub get_AA_mass_with_dyn_mod
{
	my ($self) = shift;
	my %mod_symbol = ('M'=>"@",'S'=>"#",'T'=>"%",'Y'=>"*",'G'=>"^",'K'=>"&",'D'=>'?','A'=>'~','Q'=>'!','P'=>"(",'E'=>")",'E'=>"{",'V'=>"}","V"=>"[","H"=>"]","C"=>":","F"=>",","I"=>';',"L"=>',',"R"=>"<","N"=>">","W"=>"'");
						
#	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");	
	my $parameter = $self->get_parameter();
	my $p = Spiders::Params->new();
	my ($dynamic_modif, $largest_modif) = $p->get_dynamic_modifications($parameter);
	my $AA_mass = $self->get_AA_mass();
	my $i=0;
	foreach my $modif_aa (keys %$dynamic_modif)
	{
#		$AA_mass->{$modif_aa . $mod_symbol[$i]} =  $AA_mass->{$modif_aa} + $dynamic_modif->{$modif_aa};
		$AA_mass->{$modif_aa . $mod_symbol{$modif_aa}} =  $AA_mass->{$modif_aa} + $dynamic_modif->{$modif_aa};

		$i++;
	}
	return $AA_mass;
}

sub get_Mod_mass
{
	my ($self) = shift;
	my $Mod_mass;
	
	if(!defined($self->{'_mod_mass'}))
	{	
		$Mod_mass->{'M'} = 15.99492;	
		$Mod_mass->{'Y'} = 79.96633;
		$Mod_mass->{'S'} = 79.96633;
		$Mod_mass->{'T'} = 79.96633;
	}
	return $Mod_mass;
}

sub get_immo_mass
{
	my $self=shift;
	my %immo_mass = (
	#				'F'=>1,
	#				'H'=>1,
	#				'P'=>1,
	#				'I'=>1,
					'Y'=>1,
	#				'L'=>1,
	#				'V'=>1,
					);
	return \%immo_mass;
}

sub get_mol_mass
{
	my ($self) = shift;
	my $Mol_mass;
	if(!defined($self->{'_mol_mass'}))
	{	
		$Mol_mass->{'NH3'} = 15.99492;	
		$Mol_mass->{'H'} = 1.007276466812;

	}
	return $Mol_mass;

}

=head
sub get_amino_acid_mass
{
	my ($self) = shift;
	my $AA_mass;
	my $modification = $self->get_parameter();
	
	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");
	if(!defined($self->{'_amino_acid_mass'}))
	{
		$AA_mass->{'G'} = 57.021464;	
		$AA_mass->{'D'} = 115.02694;
		$AA_mass->{'A'} = 71.037114;
		$AA_mass->{'Q'} = 128.05858;
		$AA_mass->{'S'} = 87.032029;
		$AA_mass->{'K'} = 128.09496;
		$AA_mass->{'P'} = 97.052764;
		$AA_mass->{'E'} = 129.04259;
		$AA_mass->{'V'} = 99.068414;
		$AA_mass->{'M'} = 131.04048;
		$AA_mass->{'T'} = 101.04768;
		$AA_mass->{'H'} = 137.05891;
		$AA_mass->{'C'} = 103.00919;
		$AA_mass->{'F'} = 147.06841;
		$AA_mass->{'I'} = 113.08406;
		$AA_mass->{'L'} = 113.08406;
		$AA_mass->{'J'} = 113.08406;		
		$AA_mass->{'R'} = 156.10111;
		$AA_mass->{'N'} = 114.04293;	
		$AA_mass->{'Y'} = 163.06333;
		$AA_mass->{'W'} = 186.07931;
		
		
		if(defined $modification->{'STY'} and (($modification->{'STY'} + 0.00) != 0.00))
		{
			my $S_sym = shift @mod_symbol;
			$S_sym = "S" . $S_sym;
			$AA_mass->{$S_sym} = 87.032029 + $modification->{'STY'};
			my $T_sym = shift @mod_symbol;
			$T_sym = "T" . $T_sym;			
			$AA_mass->{$T_sym} = 101.04768 + $modification->{'STY'};
			
			my $Y_sym = shift @mod_symbol;
			$Y_sym = "Y" . $Y_sym;							
			$AA_mass->{$Y_sym} = 163.06333 + $modification->{'STY'};			
		}
		if(defined $modification->{'M'} and (($modification->{'M'}+0.00) != 0.00))
		{
			my $M_sym = shift @mod_symbol;
			$M_sym = "M" . $M_sym;							
			
			$AA_mass->{$M_sym} = 131.04048+$modification->{'M'};			
		}
		if(defined $modification->{'K'})
		{
			my @num_mod = split(/\s+/,$modification->{'K'});
			
			for(my $i=0; $i<@num_mod; $i++)
			{
				next if($num_mod[$i] == 0);
				my $K_sym = shift @mod_symbol;
				$K_sym = "K" . $K_sym;							
		
				$AA_mass->{$K_sym} = 128.09496+$num_mod[$i];			
			}
		}
		if(defined $modification->{'C'} and (($modification->{'C'}+0.00) != 0.00))
		{
			my $C_sym = shift @mod_symbol;
			$C_sym = "C" . $C_sym;							
		
			$AA_mass->{$C_sym} = 103.00919+$modification->{'C'};			
		}			
	}
	
	$self->{'_amino_acid_mass'} = $AA_mass;
	return $self->{'_amino_acid_mass'};
}
=cut

sub generate_mass_combination
{
	my ($self,$mass) = @_;
	my %mass2;

	foreach my $aa (keys %$mass)
	{
		foreach my $aa2 (keys %$mass)
		{
			my $key = $aa . $aa2;
			my $rev_key = $aa2 . $aa;
			next if($mass2{$rev_key});
			$mass2{$key} = $mass->{$aa} + $mass->{$aa2}; 
		}		
	}
	return \%mass2;
}

sub get_peptide_mass
{
	my ($self,$peptSeq,$modif,$terminus,$AA_mass,$dynamic_mod)=@_;
	my $param = $self->get_parameter();

	my $H = 1.007276466812;;
	my $m = 0;
  # Computes peptide mass

  if(defined $peptSeq)
  {
		my $len = length($peptSeq);
		my @modif;
		if (defined($modif)){
			if ((my $ref = ref($modif)) eq 'ARRAY'){
			# The modifs are given as a vector directly
				@modif = @$modif;
				croak("Vector @$modif too long[".join(',',@modif)."]") if (@modif > $len+2);
			}
			elsif (!$ref){
			# Scalar, assume string
				@modif = split(/:/, $modif);
			}
			else{
				croak("Unknown format for specifying modifications [$modif]");
			}
		}

		if($terminus eq 'N')
		{
			$m = $H;

			if(defined($param->{"add_Nterm_peptide"}) and $param->{"add_Nterm_peptide"}>0)
			{

				$m += $param->{"add_Nterm_peptide"};
			}
			
		}
		elsif($terminus eq 'C')
		{
			$m=19.017806;
			if(defined($param->{"add_Cterm_peptide"}) and $param->{"add_Cterm_peptide"}>0)
			{
				$m += $param->{"add_Cterm_peptide"};
			}			
			
		}
		foreach (@modif){
			next if(!$_);
			if($dynamic_mod->{$_})
			{
				$m += $dynamic_mod->{$_};
			}
			else
			{	
				my $parameter = new Spiders::Params();
				my ($dynamic_mod,$largest_mass) = $parameter->get_dynamic_modifications($param);
				if($dynamic_mod->{$_})
				{
					$m += $dynamic_mod->{$_};
				}
				else
				{
						print %$dynamic_mod,"aa\taa",$peptSeq,"\t",$_,"\t",$terminus,"\n";
						print "Have you defined right dynamic modification?\n";
						exit;
				}
			}
		}
		if(defined $peptSeq)
		{
			foreach (split(//, $peptSeq)){
				if(!($AA_mass->{$_}))
				{
				#	print "\nfff","not existing amino acid: $peptSeq","aa\tbb",$_,"cc\n";
					last;
				}			
				$m += $AA_mass->{$_};
			}
		}
	}
    return $m; 
}


sub getFragmentMasses {
    my ($self, %h) = @_;
    my ( $pept, $modif, $frags, $spectrum ) = ( $h{pept}, $h{modif}, $h{fragTypes}, $h{spectrum} );

	my $param = $self->get_parameter();
	my $parameter = new Spiders::Params();
	my ($dynamic_mod,$largest_mass) = $parameter->get_dynamic_modifications($param);
	  
#	my $modif_mass = $self->get_Mod_mass(); 
	my $AA_mass = $self->get_AA_mass();
	my $mol_mass = $self->get_mol_mass();
	my $immo_mass = $self->get_immo_mass();

	my %fragType = %{$self->get_fragType()};
	my %series = %{$self->get_series()};
	my %loss = %{$self->get_loss()};
	
	my $immoDelta =  -26.9871;
	my $nTerm = 1.007825;

	# if mono = 0; if average =1;
	my $massType = 0;
    # Cleans the spectrum hash just in case
    undef(%$spectrum);

    # Gets peptide sequence
    my $peptSeq = $pept;
 
    # Computes peptide mass
    my $len = length($peptSeq);
    my @modif = split( /:/, $modif );
	
    my $peptideMass =$self->get_peptide_mass($peptSeq, $modif,"C",$AA_mass,$dynamic_mod)-$mol_mass->{'H'};
	
    $spectrum->{peptideMass} = $peptideMass;
    $spectrum->{peptide}     = $pept;
    $spectrum->{modif}       = [@modif];

    # Computes the sums of amino acid masses
    my @pept = split( //, $peptSeq );
    my $mass = 0;
    $mass += $dynamic_mod->{$modif[0]} if ( $modif[0] );    # N-term modif
    my @base;
    push( @base, 0 );    # for complete C-Term ions
    for ( my $i = 0 ; $i < @pept ; $i++ ) {

		
		if($AA_mass->{$pept[$i]})
		{
			$mass += $AA_mass->{$pept[$i]};
			#print $AA_mass->{$pept[$i]},"bb\n";
			$mass += $dynamic_mod->{$modif[$i+1]} if ( $modif[ $i + 1 ] );    # internal modif
			#print $dynamic_mod->{$modif[$i+1]},"aa\n" if ( $modif[ $i + 1 ] );

		}
		else
		{
			print "Not existing amino acids: $pept[$i] \n";
		}
        push( @base, $mass );
    }
    $base[-1] += $dynamic_mod->{$modif[$len+1]} if ( $modif[ $len + 1 ] );      # C-term modif

    # Computes the fragments of each requested fragment type
    foreach my $frag (@$frags) {
        if ( $frag eq 'immo' ) 
		{
			
	
            # Immonium ions
            my %already;
            for ( my $i = 1 ; $i < @pept - 1 ; $i++ ) {
                if ( defined( $immo_mass->{ $pept[$i] } ) ) {
					 if(!defined($modif[$i+1]))
					 {
						$modif[$i+1] = '';
					}
                    my $actualAA = "$pept[$i]|$modif[$i+1]";

                    next if ( defined( $already{$actualAA} ) );

#                    if (!$modif[ $i + 1 ] || (   ( $pept[$i] eq 'C' )  || ( $pept[$i] eq 'M' ) || ( $pept[$i] eq 'H' )))
#                   {
                        my $mass     = $AA_mass->{$pept[$i]} + $immoDelta;
						print $mass,"eeee\n";
                        my $immoName = $pept[$i];
                        if ( $modif[ $i + 1 ] ) {
                            $immoName .= "+$modif[$i+1]";
                            $mass += $dynamic_mod->{"$modif[$i+1]"};
                        }
						print $mass,"fff\n",;
                        $spectrum->{mass}{intern}{$frag}{$immoName} = $mass;
   #                     if ( $pept[$i] eq 'K' ) {

                            # Consider a possible extra mass with ammonia loss
     #                       $mass -= $mol_mass->{'NH3'};
      #                      $spectrum->{mass}{intern}{$frag}{"$immoName-NH3"} = $mass;
      #                  }
                        $already{$actualAA} = 1;
                        $spectrum->{ionType}{$frag}[0] = 'immo';
           #         }
                }
            }			
			
        }
        else {

            # Regular fragment types

            my $series    = $fragType{$frag}{series};
            my $charge    = $fragType{$frag}{charge};
            my $loss      = $fragType{$frag}{loss};
            my $firstFrag = $series{$series}{firstFrag};
            my $lastFrag  = $series{$series}{lastFrag};
            my $delta     = $series{$series}{delta}[$massType];

            if ( !defined($loss) ) {

                # no loss, straightforward computation
                $delta += ( $charge - 1 ) * $mol_mass->{'H'};

                if ( $series{$series}{terminus} eq 'N' ) {

                    # N-term ions
                    $delta += $nTerm;
					if($param->{"add_Nterm_peptide"}>0)
					{
						$delta += $param->{"add_Nterm_peptide"};
					}					
                    for (my $i = $firstFrag; $i <= $len - $lastFrag +1; $i++)
                    {
                        $spectrum->{mass}{term}{$frag}[ $i - 1 ] =( $base[$i] + $delta) / $charge;

                    }
                    $spectrum->{ionType}{$frag}[0] = $frag;
                }
                else 
				{

                    # C-term ions, reverse and complement masses
                    for (my $i = $firstFrag - 1 ; $i <= $len - $lastFrag +1; $i++)
                    {
				
                        $spectrum->{mass}{term}{$frag}[$i] = ( $peptideMass - $base[ $len - $i - 1 ] + $delta ) / $charge;
                    }
                    $spectrum->{ionType}{$frag}[0] = $frag;
                }
            }
            else {

                # Losses, possibly multiple and combined

                # Locates available positions for loss for each loss type (and checks all residues are different)
                my @loss = @$loss;
				my $nTotLoss=0;
                my ( @avail, @distinctPos );
                my ( $startLoop, $endLoop ) = ( $series{$series}{terminus} eq 'N' ) ? ( 0, $len - $lastFrag ) : ( $lastFrag - 1, $len - 1 );
				my $H2O_num = 0;						
				my $NH3_num = 0;
				my $HPO3_num = 0;
				my $H3PO4_num = 0;				
                for ( my $i = $startLoop ; $i <= $endLoop ; $i++ ) {
                    for ( my $j = 0 ; $j < @loss ; $j++ ) {
						next if(!defined($loss{ $loss[$j] }{residues}{ $pept[$i]}));

						if($loss[$j] eq 'H2O')
						{
							$H2O_num++;
						}
						if($loss[$j] eq 'NH3')
						{
							$NH3_num++;
						}
						if($loss[$j] eq 'HPO3')
						{
							$HPO3_num++;
						}
						if($loss[$j] eq 'H3PO4')
						{
							$H3PO4_num++;
						}
						
                        if ( $loss{ $loss[$j] }{residues}{ $pept[$i] } ) {
							next if (!defined($modif[ $i + 1 ]));
                            push( @{ $avail[$j] }, $i );
                            $nTotLoss++;
                            $distinctPos[$i]++;
                        }
                    }
                }

				
				
                next if ( $nTotLoss == 0 );    # this ion type is not possible for this peptide

                for ( my $i = 0 ; $i < @distinctPos ; $i++ ) {
					next if (!defined($distinctPos[$i]));
                    if ( $distinctPos[$i] > 1 ) {
                        croak("Multiple loss at one single amino acid [$pept[$i]] in [$frag]");
                    }
                }

                # Computes maximum number of losses for each loss type
                my @maxLoss;
                for ( my $j = 0 ; $j < @loss ; $j++ ) {
                    my $repeat = $fragType{$frag}{repeat}[$j];
                    $repeat = $len if ( $repeat == -1 );
                    $maxLoss[$j] = defined( $avail[$j] ) ? (( @{ $avail[$j] } < $repeat ) ? scalar( @{ $avail[$j] } ) : $repeat ) : 0;
                }

                # Reverses the loss positions for C-term ions
                if ( $series{$series}{terminus} eq 'C' ) {
                    for ( my $j = 0 ; $j < @loss ; $j++ ) {
                        @{ $avail[$j] } = reverse( @{ $avail[$j] } );
                    }
                }

                # Generates every combination of number of possible losses
                my $comb = 0;
                my @nLoss = split( //, '0' x @maxLoss );
                $nLoss[0] = 1;
                $delta += $nTerm if ( $series{$series}{terminus} eq 'N' );

						
				
                while (1) {

                    # Computes the fragments of the current combination

                    # Adapt delta to the number of losses
                    my $d = $delta;
					
                    for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
                        $d += $nLoss[$i] * $loss{ $loss[$i] }{delta}[$massType];
                    }

                    if ( $series{$series}{terminus} eq 'N' ) 
					{

                        # N-term ions

                        # First position that includes as many as required possible losses
                        my $rightMost = 0;
						


						if($param->{"add_Nterm_peptide"}>0)
						{
							$delta += $param->{"add_Nterm_peptide"};
						}			
						
					##### xusheng on 11/25/2014 because loss from N terminal
						if(($H2O_num > 1) and ($frag =~  /H2O/))
						{
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}
						if(($NH3_num > 0) and ($frag =~  /NH3/))
						{
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}						
						if(($HPO3_num > 1) and ($frag =~  /HPO3/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] > $rightMost ) )
								{
									$rightMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}
						if(($H3PO4_num > 1) and ($frag =~  /H3PO4/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] > $rightMost ) )
								{
									$rightMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}
					
                    }
                    else {

                        # C-term ions

                        # Last position that includes as many as required possible losses (from the right)
                       # my $leftMost = $len - 1;
						 my $leftMost = $len-1;
					 
					##### xusheng on 11/25/2014 because loss from N terminal
						#print $H2O_num,"\t",$frag,"dfdfddfd\n";
						if(($H2O_num>0) and ($frag =~ /H2O/))
						{
                        # Computes the fragments and check they are visible (firstFrag/lastFrag)
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}						
						}
						if(($NH3_num > 1) and ($frag =~ /NH3/))
						{
							# Computes the fragments and check they are visible (firstFrag/lastFrag)
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}						
						}
						if(($HPO3_num > 1) and ($frag =~ /HPO3/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] < $leftMost ))
								{
									$leftMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}								
						}
						if(($H3PO4_num > 1) and ($frag =~ /H3PO4/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] < $leftMost ))
								{
									$leftMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}							
						}						
					}

                    # Computes the exact ion type and saves its name
					my $ionType = $series;
					for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
						if ( $nLoss[$i] > 1 ) {
							$ionType .= "-$nLoss[$i]($loss[$i])";
						}
						elsif ( $nLoss[$i] == 1 ) {
							$ionType .= "-$loss[$i]";
						}
					}
					$spectrum->{ionType}{$frag}[$comb] = $ionType;
					$comb++;

					# Gets the next combination
					my $i;
					for ($i = 0 ;( $i < @maxLoss ) && ( $nLoss[$i] == $maxLoss[$i] ) ; $i++)
					{
					}
					last if ( $i == @maxLoss );
					$nLoss[$i]++;
					for ( my $j = 0 ; $j < $i ; $j++ ) 
					{
						$nLoss[$j] = 0;
					}
					
                }
            }
        }
    }

}    # getFragmentMasses

=head

sub getFragmentMasses {
    my ($self, %h) = @_;
    my ( $pept, $modif, $frags, $spectrum,$peptideMass,$AA_mass) = ( $h{pept}, $h{modif}, $h{fragTypes}, $h{spectrum},$h{peptidemass},$h{AA_mass} );

	my $param = $self->get_parameter();
	my $parameter = new Spiders::Params();
	my ($dynamic_mod,$largest_mass) = $parameter->get_dynamic_modifications($param);
	  
#	my $modif_mass = $self->get_Mod_mass(); 
#	my $AA_mass = $self->get_AA_mass();
	my $mol_mass = $self->get_mol_mass();
	my $immo_mass = $self->get_immo_mass();

	my %fragType = %{$self->get_fragType()};
	my %series = %{$self->get_series()};
	my %loss = %{$self->get_loss()};
	
	my $immoDelta =  -26.9871;
	my $nTerm = 1.007825;

	# if mono = 0; if average =1;
	my $massType = 0;
    # Cleans the spectrum hash just in case
    undef(%$spectrum);

    # Gets peptide sequence
    my $peptSeq = $pept;
 
    # Computes peptide mass
    my $len = length($peptSeq);
    my @modif = split( /:/, $modif );
## version 11.0.2	
 #   my $peptideMass =$self->get_peptide_mass($peptSeq, $modif,"C",$AA_mass,$dynamic_mod)-$mol_mass->{'H'};
	if($param->{"add_Nterm_peptide"}>0)
	{
		$peptideMass -= $param->{"add_Nterm_peptide"};
	}				
	
    $spectrum->{peptideMass} = $peptideMass;

    $spectrum->{peptide}     = $pept;
    $spectrum->{modif}       = [@modif];

    # Computes the sums of amino acid masses
    my @pept = split( //, $peptSeq );
    my $mass = 0;
    $mass += $dynamic_mod->{$modif[0]} if ( $modif[0] );    # N-term modif
    my @base;
    push( @base, 0 );    # for complete C-Term ions
    for ( my $i = 0 ; $i < @pept ; $i++ ) {

		
		if($AA_mass->{$pept[$i]})
		{
			$mass += $AA_mass->{$pept[$i]};
			#print $AA_mass->{$pept[$i]},"bb\n";
			$mass += $dynamic_mod->{$modif[$i+1]} if ( $modif[ $i + 1 ] );    # internal modif
			#print $dynamic_mod->{$modif[$i+1]},"aa\n" if ( $modif[ $i + 1 ] );

		}
		else
		{
			print "Not existing amino acids: $pept[$i] \n";
		}
        push( @base, $mass );
    }
    $base[-1] += $dynamic_mod->{$modif[$len+1]} if ( $modif[ $len + 1 ] );      # C-term modif

    # Computes the fragments of each requested fragment type
    foreach my $frag (@$frags) {
        if ( $frag eq 'immo' ) 
		{
			
	
            # Immonium ions
            my %already;
            for ( my $i = 1 ; $i < @pept - 1 ; $i++ ) {
                if ( defined( $immo_mass->{ $pept[$i] } ) ) {
					 if(!defined($modif[$i+1]))
					 {
						$modif[$i+1] = '';
					}
                    my $actualAA = "$pept[$i]|$modif[$i+1]";

                    next if ( defined( $already{$actualAA} ) );

#                    if (!$modif[ $i + 1 ] || (   ( $pept[$i] eq 'C' )  || ( $pept[$i] eq 'M' ) || ( $pept[$i] eq 'H' )))
#                   {
                        my $mass     = $AA_mass->{$pept[$i]} + $immoDelta;
					
                        my $immoName = $pept[$i];
                        if ( $modif[ $i + 1 ] ) {
                            $immoName .= "+$modif[$i+1]";
                            $mass += $dynamic_mod->{"$modif[$i+1]"};
                        }
						
                        $spectrum->{mass}{intern}{$frag}{$immoName} = $mass;
   #                     if ( $pept[$i] eq 'K' ) {

                            # Consider a possible extra mass with ammonia loss
     #                       $mass -= $mol_mass->{'NH3'};
      #                      $spectrum->{mass}{intern}{$frag}{"$immoName-NH3"} = $mass;
      #                  }
                        $already{$actualAA} = 1;
                        $spectrum->{ionType}{$frag}[0] = 'immo';
           #         }
                }
            }			
			
        }
        else {

            # Regular fragment types

            my $series    = $fragType{$frag}{series};
            my $charge    = $fragType{$frag}{charge};
            my $loss      = $fragType{$frag}{loss};
            my $firstFrag = $series{$series}{firstFrag};
            my $lastFrag  = $series{$series}{lastFrag};
            my $delta     = $series{$series}{delta}[$massType];

            if ( !defined($loss) ) {

                # no loss, straightforward computation
                $delta += ( $charge - 1 ) * $mol_mass->{'H'};
				if($charge == 2)
				{
					$delta += $mol_mass->{'H'}/2; 
				}	
                if ( $series{$series}{terminus} eq 'N' ) {

                    # N-term ions
                    $delta += $nTerm;
					if($param->{"add_Nterm_peptide"}>0)
					{
						$delta += $param->{"add_Nterm_peptide"};
					}					
                    for (my $i = $firstFrag; $i <= $len - $lastFrag +1; $i++)
                    {
                        $spectrum->{mass}{term}{$frag}[ $i - 1 ] =( $base[$i] + $delta) / $charge;

                    }
                    $spectrum->{ionType}{$frag}[0] = $frag;
                }
                else 
				{

                    # C-term ions, reverse and complement masses
                    for (my $i = $firstFrag - 1 ; $i <= $len - $lastFrag +1; $i++)
                    {
				
                        $spectrum->{mass}{term}{$frag}[$i] = ( $peptideMass - $base[ $len - $i - 1 ] + $delta ) / $charge;
                    }
                    $spectrum->{ionType}{$frag}[0] = $frag;
                }
            }
            else {

                # Losses, possibly multiple and combined

                # Locates available positions for loss for each loss type (and checks all residues are different)
                my @loss = @$loss;
				my $nTotLoss=0;
                my ( @avail, @distinctPos );
                my ( $startLoop, $endLoop ) = ( $series{$series}{terminus} eq 'N' ) ? ( 0, $len - $lastFrag ) : ( $lastFrag - 1, $len - 1 );
				my $H2O_num = 0;						
				my $NH3_num = 0;
				my $HPO3_num = 0;
				my $H3PO4_num = 0;				
                for ( my $i = $startLoop ; $i <= $endLoop ; $i++ ) {
                    for ( my $j = 0 ; $j < @loss ; $j++ ) {
						next if(!defined($loss{ $loss[$j] }{residues}{ $pept[$i]}));

						if($loss[$j] eq 'H2O')
						{
							$H2O_num++;
						}
						if($loss[$j] eq 'NH3')
						{
							$NH3_num++;
						}
						if($loss[$j] eq 'HPO3')
						{
							$HPO3_num++;
						}
						if($loss[$j] eq 'H3PO4')
						{
							$H3PO4_num++;
						}
						
                        if ( $loss{ $loss[$j] }{residues}{ $pept[$i] } ) {
							next if (!defined($modif[ $i + 1 ]));
                            push( @{ $avail[$j] }, $i );
                            $nTotLoss++;
                            $distinctPos[$i]++;
                        }
                    }
                }

				
				
                next if ( $nTotLoss == 0 );    # this ion type is not possible for this peptide

                for ( my $i = 0 ; $i < @distinctPos ; $i++ ) {
					next if (!defined($distinctPos[$i]));
                    if ( $distinctPos[$i] > 1 ) {
                        croak("Multiple loss at one single amino acid [$pept[$i]] in [$frag]");
                    }
                }

                # Computes maximum number of losses for each loss type
                my @maxLoss;
                for ( my $j = 0 ; $j < @loss ; $j++ ) {
                    my $repeat = $fragType{$frag}{repeat}[$j];
                    $repeat = $len if ( $repeat == -1 );
                    $maxLoss[$j] = defined( $avail[$j] ) ? (( @{ $avail[$j] } < $repeat ) ? scalar( @{ $avail[$j] } ) : $repeat ) : 0;
                }

                # Reverses the loss positions for C-term ions
                if ( $series{$series}{terminus} eq 'C' ) {
                    for ( my $j = 0 ; $j < @loss ; $j++ ) {
                        @{ $avail[$j] } = reverse( @{ $avail[$j] } );
                    }
                }

                # Generates every combination of number of possible losses
                my $comb = 0;
                my @nLoss = split( //, '0' x @maxLoss );
                $nLoss[0] = 1;
                $delta += $nTerm if ( $series{$series}{terminus} eq 'N' );

						
				
                while (1) {

                    # Computes the fragments of the current combination

                    # Adapt delta to the number of losses
                    my $d = $delta;
					
                    for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
                        $d += $nLoss[$i] * $loss{ $loss[$i] }{delta}[$massType];
                    }

                    if ( $series{$series}{terminus} eq 'N' ) 
					{

                        # N-term ions

                        # First position that includes as many as required possible losses
                        my $rightMost = 0;
						


						if($param->{"add_Nterm_peptide"}>0)
						{
							$delta += $param->{"add_Nterm_peptide"};
						}			
						
					##### xusheng on 11/25/2014 because loss from N terminal
						if(($H2O_num > 1) and ($frag =~  /H2O/))
						{
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}
						if(($NH3_num > 0) and ($frag =~  /NH3/))
						{
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}						
						if(($HPO3_num > 1) and ($frag =~  /HPO3/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] > $rightMost ) )
								{
									$rightMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}
						if(($H3PO4_num > 1) and ($frag =~  /H3PO4/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] > $rightMost ) )
								{
									$rightMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for ( my $i = $rightMost + 1 ; $i <= $len - $lastFrag + 1 ; $i++ )
							{
								if ( $i >= $firstFrag ) {
									$spectrum->{mass}{term}{$frag} [ ( $comb * $len ) + $i - 1 ] = ( $base[$i] + $d ) / $charge;
								}
							}							
						}
					
                    }
                    else {

                        # C-term ions

                        # Last position that includes as many as required possible losses (from the right)
                       # my $leftMost = $len - 1;
						 my $leftMost = $len-1;
					 
					##### xusheng on 11/25/2014 because loss from N terminal
						#print $H2O_num,"\t",$frag,"dfdfddfd\n";
						if(($H2O_num>0) and ($frag =~ /H2O/))
						{
                        # Computes the fragments and check they are visible (firstFrag/lastFrag)
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}						
						}
						if(($NH3_num > 1) and ($frag =~ /NH3/))
						{
							# Computes the fragments and check they are visible (firstFrag/lastFrag)
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}						
						}
						if(($HPO3_num > 1) and ($frag =~ /HPO3/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] < $leftMost ))
								{
									$leftMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}								
						}
						if(($H3PO4_num > 1) and ($frag =~ /H3PO4/))
						{
							for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
								if (   ( $nLoss[$i] > 0 ) && ( $avail[$i][ $nLoss[$i] - 1 ] < $leftMost ))
								{
									$leftMost = $avail[$i][ $nLoss[$i] - 1 ];
								}
							}
							for (my $i = $len - $leftMost - 1 ; $i < $len - $lastFrag +1; $i++)
							{
								if ( $i >= $firstFrag - 1 ) {
									$spectrum->{mass}{term}{$frag}[ ( $comb * $len ) + $i ] = ( $peptideMass - $base[ $len - $i - 1 ] + $d ) / $charge;
								}
							}							
						}						
					}

                    # Computes the exact ion type and saves its name
					my $ionType = $series;
					for ( my $i = 0 ; $i < @maxLoss ; $i++ ) {
						if ( $nLoss[$i] > 1 ) {
							$ionType .= "-$nLoss[$i]($loss[$i])";
						}
						elsif ( $nLoss[$i] == 1 ) {
							$ionType .= "-$loss[$i]";
						}
					}
					$spectrum->{ionType}{$frag}[$comb] = $ionType;
					$comb++;

					# Gets the next combination
					my $i;
					for ($i = 0 ;( $i < @maxLoss ) && ( $nLoss[$i] == $maxLoss[$i] ) ; $i++)
					{
					}
					last if ( $i == @maxLoss );
					$nLoss[$i]++;
					for ( my $j = 0 ; $j < $i ; $j++ ) 
					{
						$nLoss[$j] = 0;
					}
					
                }
            }
        }
    }

}    # getFragmentMasses
=cut


sub get_AA_combination
{
	my $self = shift;
	my %combination = (
		  'W'  => ['EG','GE','AD','DA'],
		  'R'  => ['GV','VG'],
		  'Q'  => ['AG','GA'],
		  'N'  => ['GG'],
        );
	return \%combination;	
}


sub output_mass_combination
{
	my ($self,$mass2,$mass)=@_;

	foreach my $aa ( sort {$mass2->{$a}<=>$mass2->{$b}} keys %$mass2)
	{
#		print $aa,"\t",$mass2->{$aa};
		foreach my $single_aa (keys %$mass)
		{
			if($mass->{$single_aa} < ($mass2->{$aa}+0.02) and $mass->{$single_aa} > ($mass2->{$aa}-0.02))
			{
				print "\t",$single_aa,"\t",$mass->{$single_aa};
			}
		}
		print "\n";
	}

}

sub get_mass_range
{
	my ($self,$mass) = @_;

	my $mathutil = new Spiders::MathUtils();
	
	my $max_mass = $mathutil->max($mass);
	my $min_mass = $mathutil->min($mass);
	my $range = $max_mass - $min_mass;
	return $range;
}

sub get_max_mass
{
	my ($self,$mass) = @_;

	my $mathutil = new Spiders::MathUtils();
	
	my $max_mass = $mathutil->max($mass);
	return $max_mass;
}

sub get_min_mass
{
	my ($self,$mass) = @_;
	my $mathutil = new Spiders::MathUtils();
	
	my $min_mass = $mathutil->min($mass);
	return $min_mass;
}


sub get_amino_acid_number
{
	my ($self,$mass) = @_;
	return scalar keys %$mass;
}

sub get_mis_cleavage_number
{
	my ($self,$peptide) = @_;
	
	my $mis = ($peptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;
	return $mis;
} 

sub get_mod_site_number
{
	my ($self,$peptide) = @_;
	my $mods = ($peptide =~ s/([\@\#\*\^\~\$]+)/$1/g) || 0;
	return $mods;
}
sub istryptic
{  #returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
        my ($self,$Sequence) = @_;

        my $Nterm = 0;
        my $Cterm = 0;
        if ($Sequence =~ m/^[RK]\.[A-OQ-Z].*/ || $Sequence =~ m/^[\-]\..*/) 
		{
			$Nterm = 1
		}
		if (($Sequence =~ m/.\.[A-Z\#\@\*\^\~\$]*[KR][\#\@\*\^\~\$]*\.[A-OQ-Z]\Z/) || ($Sequence =~ m/.\.[A-Z\#\@\*\^\~\$]*\.\-\Z/))
		{
			$Cterm = 1
		}
        if ($Nterm && $Cterm) 
		{
                return 3;
        }
		elsif ($Nterm || $Cterm){
                return 1 if $Nterm;
                return 2 if $Cterm;
        } else {
                return 0;
        }
}




sub get_fragType
{
	my ($self) = shift;
	
my %fragType = ('z++-NH3' => {
                       'charge' => '2',
                       'loss' => ['NH3'],
                       'repeat' => ['1'],
                       'series' => 'z'
                     },
				'y-NH3' => {
                     'charge' => '1',
                     'loss' => ['NH3'],
                     'repeat' => ['1'],
                     'series' => 'y'
                   },
				'y++' => {
                   'charge' => '2',
                   'series' => 'y'
						},
				'c-H2O*' => {
                      'charge' => '1',
                      'loss' => ['H2O'],
                      'repeat' => ['-1'],
                      'series' => 'c'
                    },
          'a' => {
                 'charge' => '1',
                 'series' => 'a'
               },
          'b-H2O*-NH3*' => {
                           'charge' => '1',
                           'loss' => ['H2O','NH3'],
                           'repeat' => ['-1','-1'],
                           'series' => 'b'
							},
          'y++-H2O' => {
                       'charge' => '2',
                       'loss' => ['H2O'],
                       'repeat' => ['1'],
                       'series' => 'y'
                     },
          'y++-H3PO4' => {
                       'charge' => '2',
                       'loss' => ['H3PO4'],
                       'repeat' => ['1'],
                       'series' => 'y'
                     },	
          'y++-HPO3' => {
                       'charge' => '2',
                       'loss' => ['HPO3'],
                       'repeat' => ['1'],
                       'series' => 'y'
                     },						 
          'a++-H2O' => {
                       'charge' => '2',
						'loss' => ['H2O'],
                       'repeat' => ['1'],
                       'series' => 'a'
                     },
          'a++-H3PO4' => {
                       'charge' => '2',
						'loss' => ['H3PO4'],
                       'repeat' => ['1'],
                       'series' => 'a'
                     },	
          'a++-HPO3' => {
                       'charge' => '2',
						'loss' => ['HPO3'],
                       'repeat' => ['1'],
                       'series' => 'a'
                     },						 
          'x++-H2O' => {
                       'charge' => '2',
                       'loss' => ['H2O'],
                       'repeat' => ['1'],
                       'series' => 'x'
                     },
          'x++-H3PO4' => {
                       'charge' => '2',
                       'loss' => ['H3PO4'],
                       'repeat' => ['1'],
                       'series' => 'x'
                     },
          'x++-HPO3' => {
                       'charge' => '2',
                       'loss' => ['HPO3'],
                       'repeat' => ['1'],
                       'series' => 'x'
                     },					 
          'y' => {
                 'charge' => '1',
                 'series' => 'y'
               },
          'b++-NH3*' => {
                        'charge' => '2',
                        'loss' => ['NH3'],
                        'repeat' => ['-1'],
                        'series' => 'b'
                      },
          'a+++' => {
                    'charge' => '3',
                    'series' => 'a'
                  },
          'z++-H2O' => {
                       'charge' => '2',
                       'loss' => ['H2O'],
                       'repeat' => ['1'],
                       'series' => 'z'
                     },
          'z++-H3PO4' => {
                       'charge' => '2',
                       'loss' => ['H3PO4'],
                       'repeat' => ['1'],
                       'series' => 'z'
                     },
          'z++-HPO3' => {
                       'charge' => '2',
                       'loss' => ['HPO3'],
                       'repeat' => ['1'],
                       'series' => 'z'
                     },						 
          'z-H2O' => {
                     'charge' => '1',
                     'loss' => ['H2O'],
                     'repeat' => ['1'],
                     'series' => 'z'
                   },
          'z-H3PO4' => {
                     'charge' => '1',
                     'loss' => ['H3PO4'],
                     'repeat' => ['1'],
                     'series' => 'z'
                   },
          'z-HPO3' => {
                     'charge' => '1',
                     'loss' => ['HPO3'],
                     'repeat' => ['1'],
                     'series' => 'z'
                   },				   
          'x++-NH3*' => {
                        'charge' => '2',
                        'loss' => ['NH3'],
                        'repeat' => ['-1'],
                        'series' => 'x'
                      },
          'x++' => {
                   'charge' => '2',
                   'series' => 'x'
                 },
          'c++' => {
                   'charge' => '2',
                   'series' => 'c'
                 },
          'b-H2O' => {
                     'charge' => '1',
                     'loss' => ['H2O'],
                     'repeat' => ['1'],
                     'series' => 'b'
                   },
          'b-H3PO4' => {
                     'charge' => '1',
                     'loss' => ['H3PO4'],
                     'repeat' => ['1'],
                     'series' => 'b'
                   },
          'b-HPO3' => {
                     'charge' => '1',
                     'loss' => ['HPO3'],
                     'repeat' => ['1'],
                     'series' => 'b'
                   },				   
          'x-H2O' => {
                     'charge' => '1',
                     'loss' => ['H2O'],
                     'repeat' => ['1'],
                     'series' => 'x'
                   },
          'x-H3PO4' => {
                     'charge' => '1',
                     'loss' => ['H3PO4'],
                     'repeat' => ['1'],
                     'series' => 'x'
                   },				   
          'x-HPO3' => {
                     'charge' => '1',
                     'loss' => ['HPO3'],
                     'repeat' => ['1'],
                     'series' => 'x'
                   },				   
          'a++-H2O*' => {
                        'charge' => '2',
                        'loss' => ['H2O'],
                        'repeat' => ['-1'],
                        'series' => 'a'
                      },
          'y++-NH3' => {
                       'charge' => '2',
                       'loss' => ['NH3'],
                       'repeat' => ['1'],
                       'series' => 'y'
                     },
          'x++-H2O*' => {
                        'charge' => '2',
                        'loss' => ['H2O'],
                        'repeat' => ['-1'],
                        'series' => 'x'
                      },
          'b++-NH3' => {
                       'charge' => '2',
                       'loss' => ['NH3'],
                       'repeat' => ['1'],
                       'series' => 'b'
                     },
          'y+++' => {
                    'charge' => '3',
                    'series' => 'y'
                  },
          'x++-NH3' => {
                       'charge' => '2',
                       'loss' => ['NH3'],
                       'repeat' => ['1'],
                       'series' => 'x'
                     },
          'b-H2O*' => {
                      'charge' => '1',
                      'loss' => ['H2O'],
                      'repeat' => ['-1'],
                      'series' => 'b'
                    },
          'b-NH3' => {
                     'charge' => '1',
                     'loss' => ['NH3'],
                     'repeat' => ['1'],
                     'series' => 'b'
                   },
          'x-H2O*' => {
                      'charge' => '1',
                      'loss' => ['H2O'],
                      'repeat' => ['-1'],
                      'series' => 'x'
                    },
          'c' => {
                 'charge' => '1',
                 'series' => 'c'
               },
          'b-NH3*' => {
                      'charge' => '1',
                      'loss' => ['NH3'],
                      'repeat' => ['-1'],
                      'series' => 'b'
                    },
          'y++-H2O*' => {
                        'charge' => '2',
                        'loss' => ['H2O'],
                        'repeat' => ['-1'],
                        'series' => 'y'
                      },
          'b' => {
                 'charge' => '1',
                 'series' => 'b'
               },
          'a-H2O' => {
                     'charge' => '1',
                     'loss' => ['H2O'],
                     'repeat' => ['1'],
                     'series' => 'a'
                   },
          'a-H3PO4' => {
                     'charge' => '1',
                     'loss' => ['H3PO4'],
                     'repeat' => ['1'],
                     'series' => 'a'
                   },
          'a-HPO3' => {
                     'charge' => '1',
                     'loss' => ['HPO3'],
                     'repeat' => ['1'],
                     'series' => 'a'
                   },				   
          'z' => {
                 'charge' => '1',
                 'series' => 'z'
               },
          'a++-NH3' => {
                       'charge' => '2',
                       'loss' => ['NH3'],
                       'repeat' => ['1'],
                       'series' => 'a'
                     },
          'c-NH3*' => {
                      'charge' => '1',
                      'loss' => ['NH3'],
                      'repeat' => ['-1'],
                      'series' => 'c'
                    },
          'b++-H2O' => {
                       'charge' => '2',
                       'loss' => ['H2O'],
                       'repeat' => ['1'],
                       'series' => 'b'
                     },
          'b++-H3PO4' => {
                       'charge' => '2',
                       'loss' => ['H3PO4'],
                       'repeat' => ['1'],
                       'series' => 'b'
                     },	
          'b++-HPO3' => {
                       'charge' => '2',
                       'loss' => ['HPO3'],
                       'repeat' => ['1'],
                       'series' => 'b'
                     },						 
          'z-H2O*' => {
                      'charge' => '1',
                      'loss' => ['H2O'],
                      'repeat' => ['-1'],
                      'series' => 'z'
                    },
          'b++-H2O*' => {
                        'charge' => '2',
                        'loss' => ['H2O'],
                        'repeat' => ['-1'],
                        'series' => 'b'
                      },
          'z++' => {
                   'charge' => '2',
                   'series' => 'z'
                 },
          'z-NH3' => {
                     'charge' => '1',
                     'loss' => ['NH3'],
                     'repeat' => ['1'],
                     'series' => 'z'
                   },
          'b+++' => {
                    'charge' => '3',
                    'series' => 'b'
                  },
          'x' => {
                 'charge' => '1',
                 'series' => 'x'
               },
          'b++' => {
                   'charge' => '2',
                   'series' => 'b'
                 },
          'c-H2O' => {
                     'charge' => '1',
                     'loss' => ['H2O'],
                     'repeat' => ['1'],
                     'series' => 'c'
                   },
          'c-H3PO4' => {
                     'charge' => '1',
                     'loss' => ['H3PO4'],
                     'repeat' => ['1'],
                     'series' => 'c'
                   },
          'c-HPO3' => {
                     'charge' => '1',
                     'loss' => ['HPO3'],
                     'repeat' => ['1'],
                     'series' => 'c'
                   },				   
          'c++-NH3' => {
                       'charge' => '2',
                       'loss' => ['NH3'],
                       'repeat' => ['1'],
                       'series' => 'c'
                     },
          'z++-NH3*' => {
                        'charge' => '2',
                        'loss' => ['NH3'],
                        'repeat' => ['-1'],
                        'series' => 'z'
                      },
          'y-H2O*-NH3*' => {
                           'charge' => '1',
                           'loss' => ['H2O','NH3'],
                           'repeat' => ['-1','-1'],
                           'series' => 'y'
                         },
          'z-NH3*' => {
                      'charge' => '1',
                      'loss' => ['NH3'],
                      'repeat' => ['-1'],
                      'series' => 'z'
                    },
          'a-NH3' => {
                     'charge' => '1',
                     'loss' => ['NH3'],
                     'repeat' => ['1'],
                     'series' => 'a'
                   },
          'z++-H2O*' => {
                        'charge' => '2',
                        'loss' => ['H2O'],
                        'repeat' => ['-1'],
                        'series' => 'z'
                      },
          'y-NH3*' => {
                      'charge' => '1',
                      'loss' => ['NH3'],
                      'repeat' => ['-1'],
                      'series' => 'y'
                    },
          'y-H2O' => {
                     'charge' => '1',
                     'loss' => ['H2O'],
                     'repeat' => ['1'],
                     'series' => 'y'
                   },
          'y-H3PO4' => {
                     'charge' => '1',
                     'loss' => ['H3PO4'],
                     'repeat' => ['1'],
                     'series' => 'y'
                   },	
          'y-HPO3' => {
                     'charge' => '1',
                     'loss' => ['HPO3'],
                     'repeat' => ['1'],
                     'series' => 'y'
                   },					   
          'y++-NH3*' => {
                        'charge' => '2',
                        'loss' => ['NH3'],
                        'repeat' => ['-1'],
                        'series' => 'y'
                      },
          'a++' => {
                   'charge' => '2',
                   'series' => 'a'
                 },
          'b-H2O-NH3' => {
                         'charge' => '1',
                         'loss' => ['H2O','NH3'],
                         'repeat' => ['1',
   '1'],
                         'series' => 'b'
                       },
          'a++-NH3*' => {
                        'charge' => '2',
                        'loss' => ['NH3'],
                        'repeat' => ['-1'],
                        'series' => 'a'
                      },
          'x-NH3' => {
                     'charge' => '1',
                     'loss' => ['NH3'],
                     'repeat' => ['1'],
                     'series' => 'x'
                   },
          'a-NH3*' => {
                      'charge' => '1',
                      'loss' => ['NH3'],
                      'repeat' => ['-1'],
                      'series' => 'a'
                    },
          'a-H2O*' => {
                      'charge' => '1',
                      'loss' => ['H2O'],
                      'repeat' => ['-1'],
                      'series' => 'a'
                    },
          'y-H2O*' => {
                      'charge' => '1',
                      'loss' => ['H2O'],
                      'repeat' => ['-1'],
                      'series' => 'y'
                    },
          'x-NH3*' => {
                      'charge' => '1',
                      'loss' => ['NH3'],
                      'repeat' => ['-1'],
                      'series' => 'x'
                    },
          'c++-NH3*' => {
                        'charge' => '2',
                        'loss' => ['NH3'],
                        'repeat' => ['-1'],
                        'series' => 'c'
                      },
          'c-NH3' => {
                     'charge' => '1',
                     'loss' => ['NH3'],
                     'repeat' => ['1'],
                     'series' => 'c'
                   },
          'c++-H2O*' => {
                        'charge' => '2',
                        'loss' => ['H2O'],
                        'repeat' => ['-1'],
                        'series' => 'c'
                      },
          'c++-H2O' => {
                       'charge' => '2',
                       'loss' => ['H2O'],
                       'repeat' => ['1'],
                       'series' => 'c'
                     },
          'c++-H3PO4' => {
                       'charge' => '2',
                       'loss' => ['H3PO4'],
                       'repeat' => ['1'],
                       'series' => 'c'
                     },
          'c++-HPO3' => {
                       'charge' => '2',
                       'loss' => ['HPO3'],
                       'repeat' => ['1'],
                       'series' => 'c'
                     }						 
        );
		return \%fragType;
}

sub get_series
{
	my %series = (
          'y' => {'lastFrag' => '2',
                  'delta' => ['1.007825','1.00797594155'],
                 'terminus' => 'C',
                 'firstFrag' => '1'
               },
          'c' => {
                 'lastFrag' => '2',
                 'delta' => ['17.026549','17.03069085415'],
                 'terminus' => 'N',
                 'firstFrag' => '2'
               },
          'a' => {
                 'lastFrag' => '2',
                 'delta' => ['-27.994915','-28.01034199508'],
                 'terminus' => 'N',
                 'firstFrag' => '2'
               },
          'b' => {
                 'lastFrag' => '2',
                 'delta' => ['0','0'],
                 'terminus' => 'N',
                 'firstFrag' => '2'
               },
          'x' => {
                 'lastFrag' => '2',
                 'delta' => ['26.98709','27.00236605353'],
                 'terminus' => 'C',
                 'firstFrag' => '1'
               },
          'z' => {
                 'lastFrag' => '2',
                 'delta' => ['-15.010899','-15.01473897105'],
                 'terminus' => 'C',
                 'firstFrag' => '1'
               }
        );
		return \%series;
}



sub get_loss
{
	my $self = shift;
	
	my %loss = (
          'NH3' => {
                   'residues' => {
                                   'Q' => 1,
                                   'N' => 1,
                                   'R' => 1
                                 },
                   'delta' => ['-17.026549','-17.03069085415']
                 },
          'H2O' => {
                   'residues' => {
                                   'S' => 1,
                                   'T' => 1,
                                   'D' => 1,
                                   'E' => 1,								   
                                 },
                   'delta' => ['-18.010565','-18.01525697318']
                 },
          'HPO3' => {
                   'residues' => {
                                   'S' => 1,
                                   'T' => 1
                                 },
                   'delta' => ['-79.96633052075','-79.96633052075']
                 },				 				 
          'H3PO4' => {
                   'residues' => {
                                   'S' => 1,
                                   'T' => 1
                                 },
                   'delta' => ['-97.984171671262','-97.984171671262']
                 }				 
        );
	return \%loss;
}

1;
