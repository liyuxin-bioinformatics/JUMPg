#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Deisotope

######### Deisotope ##########################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2012 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
# this value is used to remove the precursor ion

package Spiders::Deisotope;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_C_value =>undef,
		_H_value=>undef,
		_mass_error=>undef,
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

sub set_pho_neutral_loss
{
	my ($self,$c_value)=@_;
	$self->{_pho_value}=$c_value;	
}

sub get_pho_neutral_loss
{
	my $self=shift;
	if(!defined($self->{_pho_value}))
	{
		$self->{_pho_value}=97.984171671262;
	}
	return $self->{_pho_value};
}

sub set_dta
{
	my ($self,$dta)=@_;
	$self->{'_dta_file'} = $dta;
	return $self->{'_dta_file'};
}


sub get_dta
{
	my $self=shift;
	return $self->{'_dta_file'};
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

sub MS2_deisotope
{
	my $self=shift;
	my $H=$self->get_H_value();

	my $dta = $self->get_dta();
	my $parameter=$self->get_parameter();

# define the mz hash for storing the mz and intensity
    my %mz_hash;
# open each dta file to load the data 

    open(DTAFILE,$dta) || die "can not open the dta file: $dta";
# get the precusor mz and charge 
    my $prec_mz_charge = <DTAFILE>;
    my ($prec_mass,$prec_charge) = split(/\s+/,$prec_mz_charge);

    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
		my @data =split(/\s+/,$_);

		$mz_hash{$data[0]} = $data[1];
	}
	
### check whether there are phosphorylation losses
	my $pho_loss_num = 0;
	if(defined($parameter->{'ion_losses_MS1'}))
	{
		$pho_loss_num = $self->check_pho_loss($prec_mass,$prec_charge,\%mz_hash,10);
	}
	
	
	
# sort the mz according to the intensity

    my @mz_array;
	my %prec_mz_hash;
# remove un-fragmented ion mass (3-unit), save the precursor value in %prec_mz_hash
	my $prec_mz = (($prec_mass-$H)/$prec_charge)+$H;
	$self->remove_prec(\%mz_hash,$prec_mz,\%prec_mz_hash);

## sort mz according to the intensity


#	$self->sort_hash(\%mz_hash,\@mz_array);

	my %charge_hash;
#decharge and deiotope for each peaks 
#	foreach my $mz (@mz_array)
	foreach my $mz (reverse sort {$mz_hash{$a}<=>$mz_hash{$b}} keys %mz_hash)
	{
		if(abs($mz-216.042020900552)<(216.042020900552*$parameter->{'ppm'})/1000000)
		{

			open(OUT,">$dta.y-immo");
			print OUT "dtafile\tmz\tintensity\n";			
			print OUT $dta,"\t",$mz,"\t",$mz_hash{$mz},"\n";
			close(OUT);
		}		
		next if (!defined($mz_hash{$mz}));
# find charge and isotopes
		my $charge = $self->define_charge(\%mz_hash,$parameter->{'ppm'},$mz,$prec_charge);
## define the isotopic peaks using the charge, tolerance and orignal peaks
		next if ($charge == 0);
		if(defined($charge))
		{
			$self->deisotope(\%mz_hash,$charge,$mz,$parameter->{'ppm'});
			$charge_hash{$mz}=$charge;
		}
	}

### change the mass if the charge is larger than 1
#	print "debuging........................\n";
	foreach my $mz (keys %mz_hash)
	{
#		next if (!defined $charge_hash{$mz});
		if(defined $charge_hash{$mz} && $charge_hash{$mz} >1)
		{

				$self->changeMH(\%mz_hash,$charge_hash{$mz},$mz);

		}

	}

	my $ms2_signal_noise_ratio = $self->calculate_signal_noise_ratio($prec_mass,$prec_charge,\%mz_hash);
	
	open(OUT,">$dta");
	print OUT $prec_mz_charge;

	foreach my $mz (sort {$a<=>$b} keys %mz_hash)
	{
		if(defined($mz_hash{$mz}))
		{
			print OUT $mz," ",$mz_hash{$mz},"\n";
		}
	}
	close(OUT);
	return ($ms2_signal_noise_ratio, $pho_loss_num);
	
}

sub MS1_deisotope
{
	my $self=shift;
	my $H=$self->get_H_value();

	my $dta = $self->get_dta();
	my $parameter=$self->get_parameter();

# define the mz hash for storing the mz and intensity
    my %mz_hash;
# open each dta file to load the data 

    open(DTAFILE,$dta) || die "can not open the dta file: $dta";
# get the precusor mz and charge 
    my $prec_mz_charge = <DTAFILE>;


    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
		my @data =split(/\s+/,$_);
		$mz_hash{$data[0]} = $data[1];
	}
# sort the mz according to the intensity

    my @mz_array;
	my %prec_mz_hash;



## sort mz according to the intensity


#	$self->sort_hash(\%mz_hash,\@mz_array);

	my %charge_hash;
#decharge and deiotope for each peaks 
#	foreach my $mz (@mz_array)
	foreach my $mz (reverse sort {$mz_hash{$a}<=>$mz_hash{$b}} keys %mz_hash)
	{

		next if (!defined($mz_hash{$mz}));
# find charge and isotopes
		my $charge = $self->define_charge(\%mz_hash,$parameter->{'ppm'},$mz,3);
## define the isotopic peaks using the charge, tolerance and orignal peaks
		next if ($charge == 0);
		if(defined($charge))
		{
			#$self->deisotope(\%mz_hash,$charge,$mz,$parameter->{'ppm'});
			$charge_hash{$mz}=$charge;
		}
	}

### change the mass if the charge is larger than 1
#	print "debuging........................\n";
	foreach my $mz (keys %mz_hash)
	{
#		next if (!defined $charge_hash{$mz});
		if(defined $charge_hash{$mz} && $charge_hash{$mz} >1)
		{
				$self->changeMH(\%mz_hash,$charge_hash{$mz},$mz);
		}

	}

	#my $ms2_signal_noise_ratio = $self->calculate_signal_noise_ratio($prec_mass,3,\%mz_hash);
	
	open(OUT,">$dta");
	print OUT $prec_mz_charge;

	foreach my $mz (sort {$a<=>$b} keys %mz_hash)
	{
		if(defined($mz_hash{$mz}))
		{
			print OUT $mz," ",$mz_hash{$mz},"\n";
		}
	}
	close(OUT);
#	return $ms2_signal_noise_ratio;
}


sub calculate_signal_noise_ratio
{
	my ($self,$prec_mass,$charge,$mz_hash) = @_;
	my $charge_value = ($charge >= 2) ? 2 : 1; 
	my $number_ions_select = int(($prec_mass / 118)*0.4) * $charge_value;
	
	my @signal_intensity_array=();
	my @noise_intensity_array=();
	my $i=0;
	if(!(defined(%$mz_hash)))
	{
		return "N/A";
	}
	elsif(scalar(keys %$mz_hash)<3)
	{
		return 0;
	}
	foreach my $mz (keys %$mz_hash)
	{
		if(!defined($mz_hash->{$mz}))
		{
			delete $mz_hash->{$mz};			
		}	
		if($mz_hash->{$mz} eq '')
		{
			delete $mz_hash->{$mz};
		}
	}
	foreach my $mz (reverse sort {$mz_hash->{$a}<=>$mz_hash->{$b}} keys %$mz_hash)
	{
		next if ($mz<150);
		if($i<$number_ions_select)
		{
			push(@signal_intensity_array,$mz_hash->{$mz});
		}
		else
		{
			push(@noise_intensity_array,$mz_hash->{$mz});			
		}
		$i++;
	}
	if(scalar (@signal_intensity_array) < 2 || scalar (@noise_intensity_array)<2)
	{
		return "N/A";
	}
	
	my $signal_median = $self->median(\@signal_intensity_array);
	my $noise_median = $self->median(\@noise_intensity_array);
	if($noise_median ==0)
	{
		return "N/A";
	}
	my $signal_noise_ratio = $signal_median / $noise_median;
	return $signal_noise_ratio;
}

sub median
{
	my ($self,$array)=@_;
    my @vals = sort {$a <=> $b} @$array;
    my $len = @vals;
	return 0 if($len == 0);
	
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
		if(!defined($vals[int($len/2)]))
		{
			return $vals[int($len/2)-1];
		}
		elsif(!defined($vals[int($len/2)-1]))
		{
			return 0;
		}
		else
		{
			return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
		}
    }
}

######### Remove un-fragmented ions #################################
# remove the precursor mass according the window defined by user
# normally, the window size is 3
# save the precursor mass in %prec_mz_hash
# remove the keys and values from %mz_hash

sub remove_prec
{
	my ($self, $mz_hash,$prec_mz,$prec_mz_hash) = @_;
	my $C = $self->get_C_value();
	my $parameter=$self->get_parameter();

	my ($low, $high) = ($prec_mz- $parameter->{'prec_window'}/2*$C, $prec_mz + $parameter->{'prec_window'}/2*$C);
	while (my ($mz, $intensity) = each  %$mz_hash)
	{
		if($mz>$low && $mz<$high)
		{
			$$prec_mz_hash{$mz} = $intensity;
			delete $$mz_hash{$mz};
		}
	}
}

sub sort_hash
{
	my $self=shift;
	my ($mz_hash,$mz_array) = @_;

	my @mz_array=();
	foreach  my $mz (reverse sort {$$mz_hash{$a}<=>$$mz_hash{$b}} keys %$mz_hash)
	{
		push (@$mz_array,$mz);
	}
}


# M=1/M 

########## within one unit ################################################
# the largest charge can be obtained from precursor ion charge
# if the charge between two peaks > precusor ion charge, then go to next peak
#############################################################################



sub define_charge
{
### input #######
# $mz_hash is the reference of %mz_hash
# $mz is the specific mass that is used to define the window for it
# $ppm is the mass tolerence used for defining window
#################
	my $self=shift;
	my ($mz_hash,$ppm,$orign_mz,$prec_charge) = @_;
	my $parameter=$self->get_parameter();
	my $C = $self->get_C_value();

## The max charge can not larger than that of precursor ion
	my $maxcharge = $prec_charge;
	my $charge;

	my ($low, $high) = ($orign_mz-$C-($parameter->{'ppm'}/1000000)*$orign_mz, $orign_mz+$C+($parameter->{'ppm'}/1000000)*$orign_mz);
	my @lowarray = split('\.', $low);
	my @higharray = split('\.', $high);
	my ($lowint, $highint) = ($lowarray[0], $higharray[0]);

	my %found;
	while (my ($mz, $intensity) = each  %$mz_hash)
	{
		next if($mz<$low || $mz>$high); 
	    $found{$mz} = $intensity;
	}
  # Get mz with the highest intensity
	my $strongest_mz = $orign_mz;
	
	for my $mz (sort {$found{$b}<=>$found{$a}} keys %found)
	{
		next if ($mz==$orign_mz);
		my $diff = 1/abs($mz-$strongest_mz);
		my $round_diff = sprintf("%.0f", $diff);
		next if ($round_diff > $maxcharge);
		my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));
		next if ($var > ($parameter->{'ppm'}/1000000)*$strongest_mz);
		my $diffsum += ($var*1000000)/$mz;
		$charge = $round_diff;
		last;
	}
	if(!defined($charge))
	{
		$charge=0;
	}
	return $charge;

}

sub deisotope
{
	my $self=shift;
	my $parameter=$self->get_parameter();
	my $C = $self->get_C_value();

	my ($mz_hash,$charge,$select_mz,$ppm) = @_;

	my $isotopic_peaks_mass_error=$self->{'_mass_error'};
	my $TMT_data_start_mz = 0;
	if(defined($parameter->{'TMT_data'}))
	{
		$TMT_data_start_mz = $parameter->{'TMT_data'};		
	}	
# search the following peaks 

	my $search_loop = 1;
	my $flag=1;
	my $previous_int = $mz_hash->{$select_mz};
	my $selected_int = $mz_hash->{$select_mz};
	while($search_loop && $flag)
	{
		
		my ($peak_mz_low,$peak_mz_high) =  ($select_mz + ($C/$charge)*$search_loop-($parameter->{'ppm'}/1000000)*$select_mz, $select_mz + ($C/$charge)*$search_loop +($parameter->{'ppm'}/1000000)*$select_mz);
#		my ($peak_mz_low,$peak_mz_high) =  ($select_mz + ($C/$charge)*$search_loop - 0.02, $select_mz+ ($C/$charge)*$search_loop + 0.02);  
		
		$flag=0;
		foreach my $mz (keys %$mz_hash)
		{

			next if($mz<100);
			next if($mz<$TMT_data_start_mz);
			my $previous_mz;
			if($mz>$peak_mz_low && $mz<$peak_mz_high)
			{
				$flag=1;
				
#				if(($mz_hash->{$mz} / $previous_int) > $parameter->{'Mn_Mn1'})
				if((($mz_hash->{$mz} / $previous_int) > 0) and (($mz_hash->{$mz} / $previous_int) < 1.2))
				{
			
					$previous_int = $mz_hash->{$mz};
					$previous_mz = $mz;
					
					$mz_hash->{$select_mz} += $mz_hash->{$mz};
					$isotopic_peaks_mass_error->{$select_mz}->{'error'} = abs($select_mz - $mz) - ($C/$charge)*$search_loop;

					$isotopic_peaks_mass_error->{$select_mz}->{'intensity'} = $mz_hash->{$select_mz} + $mz_hash->{$mz};

					
					$search_loop++;

					delete $mz_hash->{$mz};
					
				}				
				else
				{
					$search_loop=0;
					$previous_mz=0;
				}
			}
		
		}
		
	}


# question is -- do we need to use charge to correct the tolerance
#	my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge)-($parameter->{'ppm'}/1000000)*$select_mz, $select_mz - ($C/$charge) +($parameter->{'ppm'}/1000000)*$select_mz);  
	my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge) - ($parameter->{'ppm'}/1000000)*$select_mz, $select_mz - ($C/$charge) + ($parameter->{'ppm'}/1000000)*$select_mz);  
	foreach my $mz (keys %$mz_hash)
	{

		next if($mz<$TMT_data_start_mz);
		next if ($select_mz eq $mz);
# search the previous peak ( only one peak) 
		if($mz>$previous_peak_mz_low && $mz<$previous_peak_mz_high)
		{
#			if(($mz_hash->{$mz} < $mz_hash->{$select_mz}))
#			{
#				delete $mz_hash->{$mz};
#			}
			if(defined ($mz_hash->{$mz}) && defined ($mz_hash->{$select_mz}) )
			{



# use the orignal intensity of selected peaks instead of peak in the hash
			#				if(($mz_hash->{$mz} / $mz_hash->{$select_mz}) > $parameter->{'M_M1'})
				if(($mz_hash->{$mz} / $mz_hash->{$select_mz}) > $parameter->{'M_M1'})
				{
			
					$mz_hash->{$mz} += $mz_hash->{$select_mz};
					$isotopic_peaks_mass_error->{$mz}->{'error'} = abs($select_mz - $mz) - $C/$charge;
					$isotopic_peaks_mass_error->{$mz}->{'intensity'} = $mz_hash->{$select_mz} + $mz_hash->{$mz};
		
					delete $mz_hash->{$select_mz};
				}
			}
		}
	}

	$self->{'_mass_error'} = $isotopic_peaks_mass_error;
}       

sub get_isotopic_peaks_mass_error
{
	my $self = shift;
	return $self->{'_mass_error'};
}

sub print_mass_error
{
	my ($self,$filename) = @_;
	my $mass_error = $self->get_isotopic_peaks_mass_error();
	open(FILE,">$filename") || die "can not open the file";
	foreach my $mz (keys %$mass_error)
	{
		
		print FILE $mz,"\t",$mass_error->{$mz}->{'intensity'},"\t",$mass_error->{$mz}->{'error'},"\n";
	}
}

sub changeMH
{
	my $self=shift;
	my $H=$self->get_H_value();
	my ($mz_hash,$charge,$mz)=@_;

	my $new_mz = ($mz-$H)*$charge+$H;
	$new_mz = sprintf("%.6f", $new_mz);
	$mz_hash->{$new_mz}=$mz_hash->{$mz};
	delete $mz_hash->{$mz};
## Merge peaks with deisotope and without deisotope		
	$self->deisotope_charge_1($mz_hash,$new_mz);
}


sub deisotope_charge_1
{
	my ($self,$mz_hash,$selected_mz) = @_;
	my $parameter=$self->get_parameter();
	my $ppm = $parameter->{'charge12_ppm'};
	
	my ($peak_mz_low,$peak_mz_high) =  ($selected_mz - ($ppm/1000000)*$selected_mz, $selected_mz +($ppm/1000000)*$selected_mz);  
#	print $peak_mz_low,"\t",$peak_mz_high,"aaa\n";
	foreach my $mz (sort {$a<=>$b} keys %$mz_hash)
	{
### check the next
		
		if($mz>$peak_mz_low && $mz<$peak_mz_high)
		{
			next if ($mz == $selected_mz);

########## need discussion about it		
#			if(($mz_hash->{$mz} / $mz_hash->{$selected_mz}) > $parameter->{'M_M1'})
#			{
				if(defined($mz_hash->{$mz}))
				{
					$mz_hash->{$selected_mz} += $mz_hash->{$mz};
### remove the un-deisotope peak				
					delete $mz_hash->{$mz};
				}
#			}
		}
	}

}	

######### add on 6/10/2013 for checking the phosophorylation neutral loss 

sub check_pho_loss
{
	my ($self,$prec_mass,$charge,$dtahash_orig,$loss_ratio) = @_;
	my $param = $self->get_parameter();
	my $pho_mass = $self->get_pho_neutral_loss();
	my $H = $self->get_H_value();


	my $error_unit = $param->{'frag_mass_tolerance_unit'};
	
	my $pho_loss_time = 0;
	my $sum_int = 0;
	my $average_int = 0;
	my $number = 0;
	foreach my $mz (keys %$dtahash_orig)
	{
		$sum_int += $dtahash_orig->{$mz};
		$number++;
	}
	$average_int = $sum_int / $number;
	my $flag=0;	
	my $flag_M1=0;
	for(my $i=1;$i<4;$i++)
	{
		#my $pho_peak = $prec_mz - $pho_mass*$i;
		my $prec_mz = ($prec_mass - $H)	/ $charge + $H;
		my $pho_peak = $prec_mz - $i*($pho_mass / $charge);


		my $prec_mz_M1 = ($prec_mass) / $charge + $H;
		my $pho_peak_M1 = $prec_mz_M1 - $i*($pho_mass / $charge);

		my $error = $param->{'frag_mass_tolerance'};		
		if($error_unit==2)
		{
			$error = $error*$prec_mz/1000000;
		}		
		foreach my $mz (keys %{$dtahash_orig})
		{
			next if (($mz-1)>$pho_peak);
		
			if(abs($pho_peak - $mz) < $error)
			{
				next if($dtahash_orig->{$mz} < ($average_int*5));				
				$flag=1;
			}
			if(abs($pho_peak_M1 - $mz) < $error)
			{

				next if($dtahash_orig->{$mz} < ($average_int*5));
				$flag_M1=1;
			}			
		}

		if($flag==1 and $flag_M1==1)
		{
			$pho_loss_time++;
			$flag=0;	
			$flag_M1=0;			
		}
	}
	return $pho_loss_time;
}



1;
