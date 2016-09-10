#!/usr/local/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::PIP

####################### PIP ###################################
#                                                             #
#       **************************************************    #  
#       **** JUMP program           		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::PIP;
        
use strict;
use warnings;
use File::Basename;
use Storable;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 2.01;
@ISA	 = qw(Exporter);
@EXPORT      = ();

#### Version 12.1.0 ############
# add a parameter to control the PPI 

########### Example ################
# use Spiders::PIP;
# my $pip = new Spiders::PIP;
# $pip->set_parameter($param);
# $pip->set_origmz_array();
# $pip->set_origmsms_hash();
# $pip->set_isolation_window();
# my $PIPref = $pip->Calculate_PIP();
########################################

sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}

sub set_parameter
{
	my ($self,$parameter) = @_;
	$self->{'_parameter'} = $parameter;
}

sub get_parameter
{
	my ($self) = shift;
	return $self->{'_parameter'};
}

sub set_origmz_array
{
	my ($self,$origmz_array) = @_;
	$self->{'_origmz_array'} = $origmz_array;
}

sub get_origmz_array
{
	my ($self) = shift;
	return $self->{'_origmz_array'};
}

sub set_origmsms_hash
{
	my ($self,$origmsms_hash) = @_;
	$self->{'_origmsms_hash'} = $origmsms_hash;
}

sub get_origmsms_hash
{
	my ($self) = shift;
	return $self->{'_origmsms_hash'};
}

sub set_isolation_window
{
	my ($self,$isolation_window) = @_;
	$self->{'_isolation_window'} = $isolation_window;
}

sub get_isolation_window
{
	my ($self) = shift;
	return $self->{'_isolation_window'};
}

sub set_isolation_offset
{
	my ($self,$isolation_window_offset) = @_;
	$self->{'_isolation_window_offset'} = $isolation_window_offset;
}

sub get_isolation_offset
{
	my ($self) = shift;
	return $self->{'_isolation_window_offset'};
}

sub set_isolation_variation
{
	my ($self,$isolation_window_variation) = @_;
	$self->{'_isolation_window_variation'} = $isolation_window_variation;
}

sub get_isolation_variation
{
	my ($self) = shift;
	return $self->{'_isolation_window_variation'};
}


sub set_dta_path
{
	my ($self,$dir)=@_;
	$self->{'_dta_path'} = $dir;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{'_dta_path'};
}


sub Calculate_PIP
{

#my $hash_dir = $ARGV[0];
#my $mzarray = retrieve($hash_dir."origmz_array");
#my $msmshash = retrieve($hash_dir."./origmsms_hash");

#	my $window = 1;
	my $self = shift;
	my %PIP;
	my $window = $self->get_isolation_window()/2;
	my $isolation_offset = $self->get_isolation_offset();
	my $isolation_variation = $self->get_isolation_variation();

	
	my $msmshash = $self->get_origmsms_hash();
	my $mzarray = $self->get_origmz_array();
	my $parameter = $self->get_parameter();
#	print "MS2_scan\tMS1_scan\tPrec_Mass\tStrongest_Mass\tCharge\tStrongest_Total_Intensity\tTotal_Intensity\tPIP\n";
	foreach my $scan (sort {$a<=>$b} keys %$msmshash)
	{
    
		my $hash = $msmshash->{$scan};

		next if (!defined($$hash{'prec_mz'}));

		my $survey = $$hash{'survey'}; 
		my $origmz = $$hash{'prec_mz'};
		my ($lowmz, $highmz) = ($origmz-$window+$isolation_offset-$isolation_variation, $origmz+$window+$isolation_offset+$isolation_variation);

		my @lowarray = split('\.', $lowmz); my @higharray = split('\.', $highmz);
		my ($lowint, $highint) = ($lowarray[0], $higharray[0]);
	
# Create a hash with all mz values between lowint and highint
		my %mzhash;
		my @mz_array;

		for (my $i = $lowint; $i<=$highint; $i++)
		{
			next if (!defined(%{$$mzarray[$survey][$i]}));
			while (my ($key, $value) = each %{$$mzarray[$survey][$i]})
			{
				if($key>=$lowmz and $key<=$highmz)
				{
					$mzhash{$key} = $value;
				}
			}
		}

		my %temp_mzhash = %mzhash;
		foreach my $mz (reverse sort {$mzhash{$a}<=>$mzhash{$b}} keys %mzhash)
		{
			push(@mz_array,$mz);
		}


###### decharge #########################

#	$parameter->{'ppm'} = 10;

	#decharge and deiotope for each peaks 
		my $intensity_sum = 0;
		my %mz_hash_saved;
		my %charge_hash;
		my %mono_hash;
		my %mono_array_hash;
		
		foreach my $mz (@mz_array)
		{
			next if (!defined($mzhash{$mz}));
	# find charge and isotopes

			my ($charge) = $self->define_charge(\%mzhash,$mzarray,$survey,$mz);

	## define the isotopic peaks using the charge, tolerance and orignal peaks
	#	next if ($charge == 0);
	#		print $charge,"\t",$mz,"hhhhh\n";	

			if(defined($charge) and $charge != 0)
			{
	#
	#		print $mzhash{$strongest_mz},"uuuuuuu\n";
				my ($mono_mz,$mono_mz_array) = $self->deisotope(\%mzhash,$charge,$mz,\%mz_hash_saved);
	#			print $scan,"\t","mono_mz_1:$mono_mz\n";
	#		print $mzhash{$strongest_mz},"vvvvvvvv\n";
				if($temp_mzhash{$mz})
				{	
					$intensity_sum += $mz_hash_saved{$mz};
					$charge_hash{$mz}=$charge;
					
					$mono_hash{$mz}=$mono_mz;
					$mono_array_hash{$mz} = $mono_mz_array;
				}
			}
			else
			{
				if($temp_mzhash{$mz})
				{	
					$charge_hash{$mz}=0;		
					$mz_hash_saved{$mz} = $mzhash{$mz};

					$mono_hash{$mz} = $mz;

					$intensity_sum += $mz_hash_saved{$mz};
					delete $mzhash{$mz};			
				}			
			}
		}

	######## sort by intensity ############
		my $i=0;

		foreach my $mz (reverse sort {$mz_hash_saved{$a}<=>$mz_hash_saved{$b}} keys %mz_hash_saved)
		{
			if($temp_mzhash{$mz})
			{
				my $duplicate = 0;
############# if same charge and same mono have been saved, ignore it ############			
				foreach my $order (keys %{$PIP{$scan}})
				{
					if($PIP{$scan}{$order}{'mono_mz'}==$mono_hash{$mz} and $PIP{$scan}{$order}{'charge'}==$charge_hash{$mz})
					{
						$duplicate = 1;
					}
				}
				if($duplicate == 0)
				{
					$i++;
					$PIP{$scan}{$i}{'mz'} = $mz;
					$PIP{$scan}{$i}{'mono_mz'} = $mono_hash{$mz};
					$PIP{$scan}{$i}{'mono_mz_array'} = $mono_array_hash{$mz};					
					$PIP{$scan}{$i}{'charge'} = $charge_hash{$mz};
								
					$PIP{$scan}{$i}{'pip'} = $mz_hash_saved{$mz} / $intensity_sum;
			#		print $scan,"\t","mono_mz_2:$PIP{$scan}{$i}{'mono_mz'}\n";
				}
				

			}
		}

		undef @mz_array;
	}
	return \%PIP;
}
#######################################################################	
=head	
	# Create a hash of only the peaks between low and high
	my %found;
	my $sum_int = 0;
	while (my ($mz, $intensity) = each  %mzhash){
	  next if($mz<$lowmz || $mz>$highmz);
	  $found{$mz} = $intensity;
	  $sum_int += $intensity;
	}
	my $largest_int = 0;
	for my $mz (sort {$found{$b}<=>$found{$a}} keys %found)
	{
		$largest_int = $found{$mz};
		last;
	}
	my $pip = $largest_int/$sum_int;
}
=cut

######################### Deisotoping #############################

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



sub define_charge
{
### input #######
# $mz_hash is the reference of %mz_hash
# $mz is the specific mass that is used to define the window for it
# $ppm is the mass tolerence used for defining window
#################
	my ($self,$mshash, $mzarray, $scan, $origmz) = @_;

	my $charge=0;

	my $parameter =$self->get_parameter();
	#my $intrappm = 10; 
	my $intrappm = $parameter->{'intrascanppm'};

	my $C = $self->get_C_value();
	my $H = $self->get_H_value();

	my $maxcharge = 5;


	my $strongest_mz = $origmz;
	my $strongest_scan = $scan;
####### Version 12.1.0: define mass region: First step to find whether there is any charge between 2 to 5 using $C/1.8
#################### 
#############
	my ($low, $high) = ($strongest_mz-$C-($intrappm/1000000)*$strongest_mz, $strongest_mz+$C+($intrappm/1000000)*$strongest_mz);
  
	my @lowarray = split('\.', $low);
	my @higharray = split('\.', $high);
	my ($lowint, $highint) = ($lowarray[0], $higharray[0]);

  # Create a hash with all mz values between lowint and highint
	my %mzhash;

   #    print "$strongest_mz low = $lowint high = $highint\n";
	for (my $i = $lowint; $i<=$highint; $i++)
	{
		next if (!defined(%{$$mzarray[$strongest_scan][$i]}));
		while (my ($key, $value) = each %{$$mzarray[$strongest_scan][$i]}){

###### selected the peaks within the exact window, not the expanded window	
			if($key>=$low and $key<=$high)
			{
#				print $key,"\t",$value,"\t",$low,"\t",$high,"\n";
				$mzhash{$key} = $value;	
			}
	
		}
	}

####### Fix a bug: the strongest peak ####################	
#	foreach my $mz (keys %mzhash)
#	{
#		if($mzhash{$mz}>$mzhash{$strongest_mz})
#		{
#			$strongest_mz = $mz;
#		}
#	}
######### get the strongest, then selected again
	
	($low, $high) = ($strongest_mz-$C-($intrappm/1000000)*$strongest_mz, $strongest_mz+$C+($intrappm/1000000)*$strongest_mz);
	@lowarray = split('\.', $low);
	@higharray = split('\.', $high);	
	($lowint, $highint) = ($lowarray[0], $higharray[0]);

	undef %mzhash;
	for (my $i = $lowint; $i<=$highint; $i++){
		next if (!defined(%{$$mzarray[$strongest_scan][$i]}));
		while (my ($key, $value) = each %{$$mzarray[$strongest_scan][$i]})
		{
				$mzhash{$key} = $value;	
		}
	}	
#################################################			
########## save mzhash for MS1 deisotope ####################
	$self->{'_deisotope_mz_hash'}=\%mzhash;	
##################################################################
	
  # Create a hash of only the peaks between low and high
	my %found;
	while (my ($mz, $intensity) = each  %mzhash){
############ check if it has been removed from the peaks within isolation window		
#		next if (!defined($mshash->{$mz}));		
		next if($mz<$low || $mz>$high);
		$found{$strongest_scan}{$mz} = $intensity;
	}

	my $diffsum=0;

  # Get mz with the highest intensity
	for my $mz (sort {$found{$strongest_scan}{$b}<=>$found{$strongest_scan}{$a}} keys %{$found{$strongest_scan}}){
		next if ($mz <= $strongest_mz);
		my $diff = 1/abs($mz-$strongest_mz);
		my $round_diff = sprintf("%.0f", $diff);
		next if ($round_diff==0);
#		print $lowint,"\t",$highint,"\t",$strongest_scan,"\t",$mz,"\t",$strongest_mz,"\t",$round_diff,"ffffffff\n";
		next if ($round_diff > $maxcharge);
		my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));
		next if ($var > ($intrappm/1000000)*$strongest_mz);
                
		$diffsum += ($var*1000000)/$mz;
		$charge = $round_diff;
		
#		$$realcharge{'dechargedscan'} = $strongest_scan;
		last;
	}
	if($charge==0 )
	{
		for my $mz (sort {$found{$strongest_scan}{$b}<=>$found{$strongest_scan}{$a}} keys %{$found{$strongest_scan}})
		{
			next if ($mz >= $strongest_mz);
			my $diff = 1/abs($mz-$strongest_mz);
			my $round_diff = sprintf("%.0f", $diff);
			next if ($round_diff==0);
#			print $lowint,"\t",$highint,"\t",$strongest_scan,"\t",$mz,"\t",$strongest_mz,"\t",$round_diff,"ggggg\n";
			next if ($round_diff > $maxcharge);
			my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));
			next if ($var > ($intrappm/1000000)*$strongest_mz);
					
			$diffsum += ($var*1000000)/$mz;
			$charge = $round_diff;
	#		$$realcharge{'dechargedscan'} = $strongest_scan;
			last;
		}	
	}
	if($charge==1)
	{
		for my $mz (sort {$found{$strongest_scan}{$b}<=>$found{$strongest_scan}{$a}} keys %{$found{$strongest_scan}})
		{
			next if ($mz >= $strongest_mz);
			my $diff = 1/abs($mz-$strongest_mz);
			my $round_diff = sprintf("%.0f", $diff);
			next if ($round_diff==0);
#			print $lowint,"\t",$highint,"\t",$strongest_scan,"\t",$mz,"\t",$strongest_mz,"\t",$round_diff,"ggggg\n";
			next if ($round_diff > $maxcharge);
			my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));
			next if ($var > ($intrappm/1000000)*$strongest_mz);
					
			$diffsum += ($var*1000000)/$mz;
			$charge = $round_diff;
## search other charges 			
			next if $charge == 1;
	#		$$realcharge{'dechargedscan'} = $strongest_scan;
			last;
		}	
	}
	return ($charge) if(defined($charge));

}


sub deisotope
{
	my ($self,$mz_hash,$charge,$select_mz,$mz_hash_saved) = @_;

	my $mono_mz=$select_mz;
	my @mono_mz_array=();
	
	
	my $parameter=$self->get_parameter();
	my $ppm = $parameter->{'intrascanppm'};
#	$parameter->{'Mn_Mn1'} = 0.5;
	my $C = $self->get_C_value();

#	my $isotopic_peaks_mass_error=$self->{'_mass_error'};
	
# search the following peaks 
	my $search_loop = 1;
	my $flag=1;
	my $previous_int = $mz_hash->{$select_mz};
	my $selected_int = $mz_hash->{$select_mz};
	$mz_hash_saved->{$select_mz} = $mz_hash->{$select_mz};
	while($search_loop && $flag)
	{
		my ($peak_mz_low,$peak_mz_high) =  ($select_mz + ($C/$charge)*$search_loop-($ppm/1000000)*$select_mz, $select_mz + ($C/$charge)*$search_loop +($ppm/1000000)*$select_mz);
#		my ($peak_mz_low,$peak_mz_high) =  ($select_mz + ($C/$charge)*$search_loop - 0.02, $select_mz+ ($C/$charge)*$search_loop + 0.02);  
		
		$flag=0;
		foreach my $mz (keys %$mz_hash)
		{
			next if($mz<100);
			my $previous_mz;
			if($mz>$peak_mz_low && $mz<$peak_mz_high)
			{
				$flag=1;
#				if(($mz_hash->{$mz} / $previous_int) < $parameter->{'Mn_Mn1'})
#				{
					$previous_int = $mz_hash->{$mz};
					$previous_mz = $mz;
					
					$mz_hash_saved->{$select_mz} += $mz_hash->{$mz};
				
#					$isotopic_peaks_mass_error->{$select_mz}->{'error'} = abs($select_mz - $mz) - ($C/$charge)*$search_loop;
#					print $isotopic_peaks_mass_error->{$select_mz}->{'error'},"\t",$select_mz,"\t",$mz,"\t",$charge,"\t",$search_loop,"\t1\n";
#					$isotopic_peaks_mass_error->{$select_mz}->{'intensity'} = $mz_hash->{$select_mz} + $mz_hash->{$mz};

#					print $select_mz,"\t",$mz,"\t",$mz_hash_saved->{$select_mz},"?????????????\n";
					$search_loop++;
					delete $mz_hash->{$mz};
					delete $mz_hash->{$select_mz};

#				}
#				else
#				{
#					$search_loop=0;
#					$previous_mz=0;
#				}
			}
		}
		
	}

# question is -- do we need to use charge to correct the tolerance
#	my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge)-($parameter->{'ppm'}/1000000)*$select_mz, $select_mz - ($C/$charge) +($parameter->{'ppm'}/1000000)*$select_mz);  
	my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge) - ($ppm/1000000)*$select_mz, $select_mz - ($C/$charge) + ($ppm/1000000)*$select_mz);  

	foreach my $mz (keys %$mz_hash)
	{
		next if ($select_mz eq $mz);
# search the previous peak ( only one peak) 

		if($mz>$previous_peak_mz_low && $mz<$previous_peak_mz_high)
		{

			if(defined ($mz_hash->{$mz}) )
			{
# use the orignal intensity of selected peaks instead of peak in the hash
			#				if(($mz_hash->{$mz} / $mz_hash->{$select_mz}) > $parameter->{'M_M1'})
#				if(($mz_hash->{$mz} / $selected_int) > $parameter->{'M_M1'})
#				{
					#$mz_hash->{$mz} += $mz_hash->{$select_mz};
					$mz_hash_saved->{$select_mz} += $mz_hash->{$mz};					
#					$isotopic_peaks_mass_error->{$mz}->{'error'} = abs($select_mz - $mz) - $C/$charge;
#					$isotopic_peaks_mass_error->{$mz}->{'intensity'} = $mz_hash->{$select_mz} + $mz_hash->{$mz};
#					print $isotopic_peaks_mass_error->{$select_mz}->{'error'},"\t",$select_mz,"\t",$mz,"\t",$charge,"\t2\n";
################### Here it differs from the deisotoping function which needs to keep the monoisotopic peak				
				#	delete $mz_hash->{$select_mz};
					delete $mz_hash->{$mz};
					delete $mz_hash->{$select_mz};
#				}
			}
		}
	}
########### creating mz_hash for PIP calculation ###########

 	
####### to find the mono isotopic peak ###############	
	my $isotopic_mz_hash = $self->{'_deisotope_mz_hash'};	

	my $max_loop = $self->get_isotopic_distribution($select_mz*$charge);

	for(my $search_loop=0; $search_loop <= $max_loop; $search_loop++)
	{

		my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge)*$search_loop - ($ppm/1000000)*$select_mz, $select_mz - ($C/$charge)*$search_loop + ($ppm/1000000)*$select_mz); 
#		print $previous_peak_mz_low,"\t",$previous_peak_mz_high,"\n";
		foreach my $mz (keys %{$isotopic_mz_hash})
		{
	# search the previous peak ( only one peak) 

			if($mz>$previous_peak_mz_low && $mz<$previous_peak_mz_high)
			{

#				if(defined ($mz_hash->{$mz}) && defined ($mz_hash->{$select_mz}) )
#				{
################################## 6/28/2013 ##################
					my $intensity_ratio = $self->get_intensity_ratio($select_mz*$charge,$search_loop);
				#	my $intensity_ratio = 0;
					if(($isotopic_mz_hash->{$mz} / $isotopic_mz_hash->{$select_mz}) > $intensity_ratio)
					{
						$mono_mz = $mz;
						push (@mono_mz_array, $mono_mz);						
					}
#				}
			}
		}
	}	
	return ($mono_mz,\@mono_mz_array);
}       

sub changeMH_folder {
	######################################################
	## Entire routine modified by JCho on 10/10/2014	##
	######################################################
	my ($self,$PIP) = @_;
	my $dir = $self->get_dta_path();
	my $H = $self->get_H_value();
	my $parameter = $self->get_parameter();
	
	my @dta = glob("$dir/*.dta");
	#my $newdir = "${dir}.10000";
	#system(qq(mkdir $newdir));
	my $newdir = $dir;

	my %charge_number;
	$charge_number{'0'}=0; $charge_number{'1'}=0; $charge_number{'2'}=0; $charge_number{'3'}=0; $charge_number{'4'}=0; $charge_number{'5'}=0;	
	my %ppi_number;
	my $dta_number = scalar @dta;

	my $msmshash = $self->get_origmsms_hash();
	my $dtafile_basename = basename($dir);
	$dtafile_basename =~s/(.*)\.(\d+)/$1/;
	my $dtafile;

	foreach my $scan (keys %{$PIP}) {
		foreach my $order (keys %{$PIP->{$scan}}) {
			my $prec_mz =  $PIP->{$scan}->{$order}->{'mono_mz'};
			my $prec_mz_array =  $PIP->{$scan}->{$order}->{'mono_mz_array'};
			my $charge = $PIP->{$scan}->{$order}->{'charge'};

######### 6/12/2013 ##### control the number of ppi selected
			if($parameter->{'max_num_ppi'} > 0) {
				if($order <= $parameter->{'max_num_ppi'}) {
					$charge_number{$charge}++;
					if($charge == 0) {
						if($parameter->{'ppi_charge_0'}) {
							for (my $k = 2; $k <= 3; $k++) {
								my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
								my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
								while ($charge < 4) {
									my $temp_MH = sprintf("%.5f",(($prec_mz - $H) * $charge + $H));
									my ($larger,$smaller) = (0,0);
									for (my $i = 0; $i <= $#msms_mz_array; $i++) {
										if($temp_MH < $msms_mz_array[$i]) {
											$larger++;
										} else {
											$smaller++;
										}
									}
									next if (($larger + $smaller) == 0);
									my $fraction = $larger / ($larger + $smaller);
									if($fraction > 0.05) {
										$charge++; 
									} else {
										last;
									}
								}
								$dtafile = $dtafile_basename.".".$scan.".".$order.".".$k.".dta";
								$ppi_number{$order}++;
								open (OUT, ">$newdir/$dtafile");
								my $MH = sprintf("%.5f",(($prec_mz - $H) * $k + $H));
								print OUT "$MH $k\n";
								for (my $i = 0; $i <= $#msms_mz_array; $i++) {
									print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
								}
								close(OUT);
							}
						}
					} elsif($charge == 1) { ######## when charge ==1, we got more decoy than target for peptides with high scores 
						if($parameter->{'ppi_charge_1'}) {
							my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
							my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
							while(1) {
								my $temp_MH = sprintf("%.5f",(($prec_mz - $H) * $charge + $H));
								my ($larger, $smaller) = (0,0);
								for (my $i = 0;$i <= $#msms_mz_array; $i++) {
									if ($temp_MH < $msms_mz_array[$i]) {
										$larger++;
									} else {
										$smaller++;
									}
								}
								my $fraction = $larger / ($larger + $smaller);
								if($fraction > 0.05) {
									$charge++; 
								} else {
									last;
								}
							}
							$charge_number{$charge}++;
							$dtafile = $dtafile_basename.".".$scan.".".$order.".".$charge.".dta";
							$ppi_number{$order}++;
							open (OUT, ">$newdir/$dtafile");
							my $MH = sprintf("%.5f",(($prec_mz - $H) * $charge + $H));
							print OUT "$MH $charge\n";
							for (my $i = 0; $i <= $#msms_mz_array; $i++) {
								print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
							}
							close(OUT);
						}
					} else {
						my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
						my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
						$dtafile = $dtafile_basename.".".$scan.".".$order.".".$charge.".dta";
						$charge_number{$charge}++;
						$ppi_number{$order}++;
						open (OUT, ">$newdir/$dtafile");
						my $MH = sprintf("%.5f",(($prec_mz - $H) * $charge + $H));
						print OUT "$MH $charge\n";
						for (my $i = 0; $i <= $#msms_mz_array; $i++) {
							print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
						}
						close(OUT);
=head
						for(my $j=0;$j<=$#$prec_mz_array;$j++)
						{
							next if($prec_mz_array->[$j] == $prec_mz);
							my $order_mono = 100 * $order + 1 + $j;
							$dtafile = $2 . "." . $specscan . "." . $order_mono . "." . $charge . ".dta";
							$charge_number{$charge}++;
							$ppi_number{$order}++;
							open (OUT, ">$newdir/$dtafile");
							shift @inarray;
							my $MH = sprintf("%.5f",(($prec_mz_array->[$j]-$H)*$charge+$H));
							print OUT "$MH $charge\n";
							print OUT @inarray;
							close(OUT);							
						}
=cut						
					}
				}				
#####################################################			
			} elsif(($PIP->{$scan}->{$order}->{'pip'}*100)>$parameter->{'percentage_ppi'} and $parameter->{'max_num_ppi'} == 0) {
				if($charge == 0) {
					$charge_number{$charge}++;
					if($parameter->{'ppi_charge_0'}) {
						my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
						my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
						for (my $k = 2; $k <= 3; $k++) {
							$dtafile = $dtafile_basename.".".$scan.".".$order.".".$k.".dta";
							$ppi_number{$order}++;
							open (OUT, ">$newdir/$dtafile");
							my $MH = sprintf("%.5f",(($prec_mz - $H) * $k + $H));
							print OUT "$MH $k\n";
							for (my $i = 0; $i <= $#msms_mz_array; $i++) {
								print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
							}
							close(OUT);
						}
					}
				} elsif($charge == 1 and $parameter->{'ppi_charge_1'}) { ######## when charge ==1, we got more decoy than target for peptides with high scores
					my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
					my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
					$charge_number{$charge}++;
					$dtafile = $dtafile_basename.".".$scan.".".$order.".".$charge.".dta";
					$ppi_number{$order}++;
					open (OUT, ">$newdir/$dtafile");
					my $MH = sprintf("%.5f",(($prec_mz - $H) * $charge + $H));
					print OUT "$MH $charge\n";
					for (my $i = 0; $i <= $#msms_mz_array; $i++) {
						print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
					}
					close(OUT);
				}
				else 
				{
					my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
					my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
					$dtafile = $dtafile_basename.".".$scan.".".$order.".".$charge.".dta";
					$charge_number{$charge}++;
					$ppi_number{$order}++;
					open (OUT, ">$newdir/$dtafile");
					my $MH = sprintf("%.5f",(($prec_mz-$H)*$charge+$H));
					print OUT "$MH $charge\n";
					for (my $i = 0; $i <= $#msms_mz_array; $i++) {
						print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
					}
					close(OUT);
=head
					for(my $j=0;$j<=$#$prec_mz_array;$j++)
					{
						next if($prec_mz_array->[$j] == $prec_mz);
						my $order_mono = 100 * $order + 1 + $j;
						$dtafile = $2 . "." . $specscan . "." . $order_mono . "." . $charge . ".dta";
						$charge_number{$charge}++;
						$ppi_number{$order}++;
						open (OUT, ">$newdir/$dtafile");
						shift @inarray;
						my $MH = sprintf("%.5f",(($prec_mz_array->[$j]-$H)*$charge+$H));
						print OUT "$MH $charge\n";
						print OUT @inarray;
						close(OUT);
					}
=cut
				}
			} elsif($parameter->{'max_num_ppi'} < 0) {
				next if ($order <= 1);
				if($order <= abs($parameter->{'max_num_ppi'})) {
					$charge_number{$charge}++;
					if($charge == 0) {
						if($parameter->{'ppi_charge_0'}) {
							my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
							my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
							for (my $k = 2; $k <= 3; $k++) {
								$dtafile = $dtafile_basename.".".$scan.".".$order.".".$k.".dta";
								$ppi_number{$order}++;
								open (OUT, ">$newdir/$dtafile");
								my $MH = sprintf("%.5f",(($prec_mz - $H) * $k + $H));
								print OUT "$MH $k\n";
								for (my $i = 0; $i <= $#msms_mz_array; $i++) {
									print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
								}
								close(OUT);
							}
						}
					} elsif($charge == 1) {######## when charge ==1, we got more decoy than target for peptides with high scores
						$charge_number{$charge}++;
						if($parameter->{'ppi_charge_1'}) {
							my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
							my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
							$dtafile = $dtafile_basename.".".$scan.".".$order.".".$charge.".dta";
							$ppi_number{$order}++;
							open (OUT, ">$newdir/$dtafile");
							my $MH = sprintf("%.5f",(($prec_mz-$H)*$charge+$H));
							print OUT "$MH $charge\n";
							for (my $i = 0; $i <= $#msms_mz_array; $i++) {
								print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
							}
							close(OUT);
						}
					} else {
						my @msms_mz_array = @{$msmshash -> {$scan} -> {'msms_mz'}};
						my @msms_int_array = @{$msmshash -> {$scan} -> {'msms_int'}};
						$dtafile = $dtafile_basename.".".$scan.".".$order.".".$charge.".dta";
						$charge_number{$charge}++;
						$ppi_number{$order}++;
						open (OUT, ">$newdir/$dtafile");
						my $MH = sprintf("%.5f",(($prec_mz-$H)*$charge+$H));
						print OUT "$MH $charge\n";
						for (my $i = 0; $i <= $#msms_mz_array; $i++) {
							print OUT "$msms_mz_array[$i] $msms_int_array[$i]\n";
						}
						close(OUT);
=head
						for(my $j=0;$j<=$#$prec_mz_array;$j++)
						{
							next if($prec_mz_array->[$j] == $prec_mz);
							my $order_mono = 100 * $order + 1 + $j;
							$dtafile = $2 . "." . $specscan . "." . $order_mono . "." . $charge . ".dta";
							$charge_number{$charge}++;
							$ppi_number{$order}++;
							open (OUT, ">$newdir/$dtafile");
							shift @inarray;
							my $MH = sprintf("%.5f",(($prec_mz_array->[$j]-$H)*$charge+$H));
							print OUT "$MH $charge\n";
							print OUT @inarray;
							close(OUT);							
						}
=cut						
					}
				}
			}
		}
	}
	my $charge_dist = "";
	my $ppi_dist = "";
	my $sum_charge = 0;

	foreach my $key (sort keys %charge_number) {
		$sum_charge += $charge_number{$key};
		$charge_dist .= "$key = $charge_number{$key} ";
	}
	print "  Decharged $sum_charge precursor ions: 0=$charge_number{'0'} 1=$charge_number{'1'} 2=$charge_number{'2'} 3=$charge_number{'3'} 4=$charge_number{'4'} 5=$charge_number{'5'}\n";
	print "  PPI distribution: ";
	foreach my $ppi_key (sort {$a<=>$b} keys %ppi_number) {
		print "$ppi_key = $ppi_number{$ppi_key} ";
		$ppi_dist .= "$ppi_key = $ppi_number{$ppi_key} ";
	}
	print "\n";
#	system(qq(rm -rf $dir/*));
#	system(qq(cp $dir/sequest.params $newdir >/dev/null 2>&1));
#	system(qq(rm -R $dir/ >/dev/null 2>&1));
#	system(qq(rm -rf $dir >/dev/null 2>&1));	
#	system(qq(mv $newdir $dir >/dev/null 2>&1));
	return ($charge_dist,$ppi_dist);
}



sub get_isotopic_distribution
{
	my ($self, $mz) = @_;
	my $loop = 0;
	if($mz<1500)
	{
		$loop = 1;
	}
	elsif($mz<3000 && $loop<=2)
	{
		$loop = 2;
	}
	elsif($mz<4500 && $loop<=4)
	{
		$loop = 4;
	}
	elsif($mz<6000 && $loop<=6)
	{
		$loop = 6;
	}	
	return $loop;	
	
}

sub get_intensity_ratio
{
	my ($self,$mz,$loop) = @_;
	my $ratio=2;
### if the mass <1500, there is no isotopic peaks preceeding the monoisotope	
	if($mz<1500)
	{
		$ratio = 0.3;
	}
	elsif($mz<3000 && $loop<=2)
	{
		if($loop == 1)
		{
			$ratio = 0.2;
		}
		elsif($loop == 2)
		{
			$ratio = 0.1;
		}
	}
	elsif($mz<4500 && $loop<=4)
	{
		if($loop == 1)
		{
			$ratio = 0.2;
		}
		elsif($loop == 2)
		{
			$ratio = 0.1;
		}
		elsif($loop == 3)
		{
			$ratio = 0.1;
		}
		elsif($loop == 4)
		{
			$ratio = 0.1;
		}
	}
	elsif($mz<6000 && $loop<=6)
	{
		if($loop == 1)
		{
			$ratio = 0.2;
		}
		elsif($loop == 2)
		{
			$ratio = 0.1;
		}
		elsif($loop == 3)
		{
			$ratio = 0.1;
		}
		elsif($loop == 4)
		{
			$ratio = 0.1;
		}
		elsif($loop == 5)
		{
			$ratio = 0.1;
		}
		elsif($loop == 6)
		{
			$ratio = 0.1;
		}
	}	
	return $ratio;
}



1;
	
	


	



	

				
