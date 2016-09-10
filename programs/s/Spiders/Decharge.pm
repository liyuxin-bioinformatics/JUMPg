#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0/
## Module name: Spiders::Decharge
############################################################# 
#                                                           #
#       **************************************************  #    
#       **** Decharge program for MS1 and MS2         ****  #    
#       ****                                          ****  #    
#       ****Copyright (C) 2012 - Xusheng Wang        ****  #    
#       ****all rights reserved.                      ****  #    
#       ****xusheng.wang@stjude.org                   ****  #    
#       ****                                          ****  #    
#       ****                                          ****  #    
#       **************************************************  #    
#############################################################

package Spiders::Decharge;

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



sub decharge{
	my ($self, $dtafile, $dta_hash, $msmshash, $mshash, $mzarray, $realcharge) = @_;

	# Decharge
	my $source = "prec";
	my $sourcemz = "prec_mz";
	my $sourceint = "prec_int";
	
	$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
#	print $dtafile,"\n";
	my ($scan, $specscan) = ($3, $3);
	
	$specscan =~ s/^0*//;

	my ($survey, $mz) = ($$msmshash{$specscan}{'survey'}, $$msmshash{$specscan}{$sourcemz});

	$$realcharge{'origscan'} = $scan;
	$$realcharge{'origmz'} = $mz;

  # Use original scan to begin deisotoping
	my $strongest_mz = $mz;
	my $strongest_scan = $scan;
	$$realcharge{'dechargedmz'} = $strongest_mz;

	my $charge = $self->find_charge($mshash, $mzarray, $survey, $mz);
=head
    if (defined($charge)){
		$$dta_hash{$scan}{'decharged'} = $charge;
		if($charge>1)
		{
			$$dta_hash{$scan}{'decharged'} = $charge;
			$dtafile = $self->changeMH($source, $dta_hash, $msmshash, $scan, $charge);
		}
		else
		{
			$$dta_hash{$scan}{'decharged'} = 2;
			$dtafile = $self->changeMH($source, $dta_hash, $msmshash, $scan, 2);		
		}
	}
=cut
#	return $dtafile;
	return $charge;
}

sub MS1_deisotope
{
	my ($self,$select_mz,$charge,$mz_hash)=@_;

	my $C = $self->get_C_value();

	my $parameter = $self->get_parameter();
### keep consistent with parameter file	
	my $peptide_tolerance = $parameter->{'frag_mass_tolerance'};

	my $prec_mz = $select_mz;
	my $flag = 1;
	my $max_loop = $self->get_isotopic_distribution($select_mz*$charge);
	
	for(my $search_loop=0; $search_loop < $max_loop; $search_loop++)
	{
		my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge)*$search_loop - $peptide_tolerance, $select_mz - ($C/$charge)*$search_loop + $peptide_tolerance); 
#		print $previous_peak_mz_low,"\t",$previous_peak_mz_high,"\n";
		foreach my $mz (keys %{$mz_hash})
		{
	# search the previous peak ( only one peak) 

			if($mz>$previous_peak_mz_low && $mz<$previous_peak_mz_high)
			{

				if(defined ($mz_hash->{$mz}) && defined ($mz_hash->{$select_mz}) )
				{
					my $intensity_ratio = $self->get_intensity_ratio($select_mz*$charge,$search_loop);

					if(($mz_hash->{$mz} / $mz_hash->{$select_mz}) > $intensity_ratio)
					{
						$prec_mz = $mz;
					}
				}
			}
		}
	}	
	return $prec_mz;
}

sub get_isotopic_distribution
{
	my ($self, $mz) = @_;
	my $loop = 0;
	if($mz<1500)
	{
		$loop = 0;
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
		$ratio = 0.8;
	}
	elsif($mz<3000 && $loop<=2)
	{
		if($loop == 1)
		{
			$ratio = 0.4;
		}
		elsif($loop == 2)
		{
			$ratio = 0.3;
		}
	}
	elsif($mz<4500 && $loop<=4)
	{
		if($loop == 1)
		{
			$ratio = 0.6;
		}
		elsif($loop == 2)
		{
			$ratio = 0.2;
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
			$ratio = 0.6;
		}
		elsif($loop == 2)
		{
			$ratio = 0.3;
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



sub find_charge{
	my ($self, $mshash, $mzarray, $scan, $origmz) = @_;

	my $charge;

	my $parameter = $self->get_parameter();
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

 #      print "$strongest_mz low = $lowint high = $highint  aaaaaaaaaaaa\n";
	for (my $i = $lowint; $i<=$highint; $i++){
		next if (!defined(%{$$mzarray[$strongest_scan][$i]}));
		while (my ($key, $value) = each %{$$mzarray[$strongest_scan][$i]}){
###### selected the peaks within the exact window, not the expanded window	
			if($key>=$low and $key<=$high)
			{
	#			print $key,"\t",$value,"\t",$low,"\t",$high,"\n";
				$mzhash{$key} = $value;	
			}
	
		}
	}

####### Fix a bug: the strongest peak ####################	
	foreach my $mz (keys %mzhash)
	{
		if($mzhash{$mz}>$mzhash{$strongest_mz})
		{
			$strongest_mz = $mz;
		}
	}
######### get the strongest, then selected again
	
	($low, $high) = ($strongest_mz-$C-($intrappm/1000000)*$strongest_mz, $strongest_mz+$C+($intrappm/1000000)*$strongest_mz);
	@lowarray = split('\.', $low);
	@higharray = split('\.', $high);	
	($lowint, $highint) = ($lowarray[0], $higharray[0]);
  #     print "$strongest_mz low = $lowint high = $highint bbbbbbbbbbbbbbbbbb\n";
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
		next if($mz<$low || $mz>$high);
		$found{$strongest_scan}{$mz} = $intensity;
	}

	my $diffsum=0;

  # Get mz with the highest intensity
	for my $mz (sort {$found{$strongest_scan}{$b}<=>$found{$strongest_scan}{$a}} keys %{$found{$strongest_scan}}){
		next if ($mz==$strongest_mz);
		my $diff = 1/abs($mz-$strongest_mz);
		my $round_diff = sprintf("%.0f", $diff);
		next if ($round_diff==0);
	#	print $lowint,"\t",$highint,"\t",$strongest_scan,"\t",$mz,"\t",$strongest_mz,"\t",$round_diff,"ffffffff\n";
		next if ($round_diff > $maxcharge);
		my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));
		next if ($var > ($intrappm/1000000)*$strongest_mz);
                
		$diffsum += ($var*1000000)/$mz;
		$charge = $round_diff;
		
#		$$realcharge{'dechargedscan'} = $strongest_scan;
		last;
	}
	
	
	return $charge if(defined($charge));
	
  # If not decharged, examine forward and backward scanrange survey scans
	if (!defined($charge)){
	
	
		$charge = $self->find_charge_by_more_scans($mshash,$mzarray,$strongest_mz,$scan);
		return $charge;
	}
}

sub find_charge_by_more_scans
{
	my ($self,$mshash,$mzarray,$strongest_mz,$scan) = @_;

	my $charge;
	my $strongest_scan = $scan;
	my $maxcharge = 4;

	my $parameter = $self->get_parameter();
	my $interppm = $parameter->{interscanppm};
	my $intrappm = $parameter->{intrascanppm};

	my $scanrange = $parameter->{preproc_msscanrange};
	$scanrange = 4 if (!defined($scanrange));

	my $C = $self->get_C_value();

	my ($low, $high) = ($strongest_mz-$C-($intrappm/1000000)*$strongest_mz, $strongest_mz+$C+($intrappm/1000000)*$strongest_mz);
	my @lowarray = split('\.', $low);
	my @higharray = split('\.', $high);
	my ($lowint, $highint) = ($lowarray[0], $higharray[0]);

	my $strongest_order = $$mshash{'surveyhash'}{$strongest_scan}{'scan_order'};

  	for (my $cycle=1;$cycle<=2;$cycle++){
		for (my $i=1;$i<=$scanrange;$i++){
			my $nextorder = $strongest_order+$i;
			$nextorder = $strongest_order-$i if ($cycle==2);
			last if (!defined($$mshash{'orderhash'}{$nextorder}));
			my $nextscan = $$mshash{'orderhash'}{$nextorder};
			my %hash;
			for (my $j = $lowint; $j<=$highint; $j++){
				next if (!defined(%{$$mzarray[$nextscan][$j]}));
				while (my ($key, $value) = each %{$$mzarray[$nextscan][$j]}){
					$hash{$key} = $value;
				}
			}
			
			my %possible; 
			my $precursor = 0;
			while (my ($mz, $intensity) = each  %hash){
				next if($mz<$low || $mz>$high);
				$possible{$nextscan}{$mz} = $intensity;
				$precursor = $mz if ($mz>$strongest_mz-($interppm/1000000)*$strongest_mz && $mz<$strongest_mz+($interppm/1000000)*$strongest_mz);
			}
			next if ($precursor == 0);
			if (defined($possible{$nextscan})){
				my $prev_charge = 0;
				for my $mz (sort keys %{$possible{$nextscan}}){
					next if ($mz == $precursor);
					my $diff = 1/abs($mz-$precursor);
					my $round_diff = sprintf("%.0f", $diff);
					next if ($round_diff > $maxcharge);
					my $var = abs(abs($mz-$precursor)-($C/$round_diff));
					next if ($var > ($intrappm/1000000)*$precursor);
					$prev_charge = $round_diff if ($prev_charge == 0);
					next if ($prev_charge == 1);
					$charge = $round_diff;
########## save mzhash for MS1 deisotope ####################
					$self->{'_deisotope_mz_hash'}=\%hash;	
##################################################################					
					#	$$realcharge{'dechargedscan'} = $nextscan;
					last;
				}
				if ($prev_charge == 1 && !defined($charge)){
					$charge = $prev_charge;
					#$$realcharge{'dechargedscan'} = $nextscan;
				}
			}
			return $charge if (defined($charge));
		}
	}
}


sub changeMH{
	my ($self, $source,  $dta_hash, $msmshash, $scan, $charge) = @_;
	my $dir = $self->get_dta_path();
	my $H = $self->get_H_value();

    my $sourcemz = $source."_mz";
        # Change the MH inside all dta file, keep four decimal points
    my ($num, $total) = (0, scalar(keys %$dta_hash));
	my $specscan = $scan;
	$specscan =~ s/^0*//;
	my $mz = $$msmshash{$specscan}{$sourcemz};
	
######### deisotoping #################

	$mz = $self->MS1_deisotope($mz,$charge);

#######################################
	
	my $previous_charge = 1;
	my $dtafile = $$dta_hash{$scan}{'charges'}{$previous_charge}{'file'};
	return 0 if(!defined($dtafile));
	open (IN, "<$dtafile");
	my @inarray = <IN>;
	close IN;
	system(qq(rm "$dtafile"));

	$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
	$dtafile = $2 . "." . $3 . "." . $3 . "." . $charge . ".dta";

	open (OUT, ">$dir/$dtafile");
	shift @inarray;
	my $MH = sprintf("%.5f",(($mz-$H)*$charge+$H));
	print OUT "$MH $charge\n";
	print OUT @inarray;
	close(OUT);
	$$dta_hash{$scan}{'charges'}{$charge}{'file'} = "$dir/$dtafile";

	return $$dta_hash{$scan}{'charges'}{$charge}{'file'};
}


sub changeMH_folder
{
	my ($self, $msmshash, $realcharge,$mzhash) = @_;
	my $dir = $self->get_dta_path();
	my $H = $self->get_H_value();
	my @dta = glob("$dir/*.dta");
	
	my $newdir = $dir . ".1";
	system(qq(mkdir $newdir));
	
	foreach my $dtafile (@dta)
	{	
#		print $dtafile,"\n";
		$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
		
		my $specscan = $3;
		my $orig_charge = $5;
		
		my $mz = $$msmshash{$specscan}{prec_mz};
######### deisotoping #################
#		print $dtafile,"\t",$mz,"\t",$realcharge->{$specscan},"\t";
		my $updated_mz = $self->MS1_deisotope($mz,$realcharge->{$specscan},$mzhash->{$specscan});
#		print $updated_mz,"\n";		
#######################################
#		next if ($updated_mz == $mz and $realcharge->{$specscan} == $orig_charge);
		open (IN, "<$dtafile");
		my @inarray = <IN>;
		close IN;
#		system(qq(rm "$dtafile"));
		if(defined($realcharge->{$specscan}) and $realcharge->{$specscan} ne "")
		{
			$dtafile = $2 . "." . $specscan . "." . $specscan . "." . $realcharge->{$specscan} . ".dta";

			open (OUT, ">$newdir/$dtafile");
			shift @inarray;
			my $MH = sprintf("%.5f",(($updated_mz-$H)*$realcharge->{$specscan}+$H));
			print OUT "$MH $realcharge->{$specscan}\n";
			print OUT @inarray;
			close(OUT);
		}
		else
		{
=head		
			$dtafile = $2 . "." . $specscan . "." . $specscan . "." . $orig_charge . ".dta";

			open (OUT, ">$newdir/$dtafile");
			print OUT @inarray;
			close(OUT);	
=cut
			$dtafile = $2 . "." . $specscan . "." . $specscan . "." . "2" . ".dta";

			open (OUT, ">$newdir/$dtafile");
			shift @inarray;
			my $MH = sprintf("%.5f",(($mz-$H)*2+$H));
			print OUT "$MH 2\n";
			print OUT @inarray;
			close(OUT);			
		}
	}
	system(qq(rm -rf $dir));
	system(qq(mv $newdir $dir));	
}


sub create_dtahash
{
	my ($self,$dtafile,$msmshash) = @_;

=head
	my $dir = $self->get_dta_path();

	opendir(DIR, $dir) || die "can not open the directory: $dir\n"; 
	my @dta = grep {/\.dta/} readdir(DIR);
	closedir(DIR);


	for my $dtafile (@dta)
	{
=cut	
		my %dta_hash;


		$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
		my ($scan, $specscan, $specscan2, $charge) = ($3, $3, $4, $5); 
		$specscan =~ s/^[0]+//;
		if ($specscan != $specscan2){   
			print "$dtafile\n";     
		}
		if (!defined($scan) || !defined($specscan) || !defined($charge)){       
			print "$dtafile not defined\n"; 
			exit;   
		}
	
		$dta_hash{$scan}{'charges'}{$charge}{'prec_mz'} = $$msmshash{$specscan}{'prec_mz'};
		$dta_hash{$scan}{'charges'}{$charge}{'intensity'} = $$msmshash{$specscan}{'prec_int'};
		$dta_hash{$scan}{'charges'}{$charge}{'peak_num'} = $$msmshash{$specscan}{'peak_num'};
		$dta_hash{$scan}{'charges'}{$charge}{'lpeak_num'} = $$msmshash{$specscan}{'lpeak_num'};
		$dta_hash{$scan}{'charges'}{$charge}{'file'} = $dtafile;
		
#	}
	return \%dta_hash;
}




1;
