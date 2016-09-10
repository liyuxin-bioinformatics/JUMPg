#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Consolidation

######### Simulation ##########################################
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
# added the normalization function (v1.0.1) on 5/11/2012

package Spiders::Consolidation;
        
use strict;
use warnings;
use Spiders::Dta;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.01;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
        _dta_file => $arg{'-dta_file'},
		_keepnum =>$arg{'-keepnum'},
    };
    bless $self, $class;
	return $self;
}

sub get_dta
{
	my $self=shift;
	return $self->{_dta_file};
}

sub get_keepnum
{
	my $self = shift;
	return $self->{_keepnum};
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
		$self->{_pho_value}=97.99532;
	}
	return $self->{_pho_value};
}


sub get_strong_peak_window
{
	my $self=shift;
	return $self->{'_strong_peak_window'};
}

sub get_fragment_tolerance
{
	my ($self,$mass)=@_;
	my $param = $self->get_parameter();
	my $frag_tolerance = $param->{'frag_mass_tolerance'};
	
	my $tolerance_unit = 1;
	if(defined($param->{'frag_mass_tolerance_unit'}))
	{
		$tolerance_unit = $param->{'frag_mass_tolerance_unit'};
	}
	if($tolerance_unit == 2)
	{
		$frag_tolerance = $frag_tolerance * $mass / 1000000;
	}
	return $frag_tolerance;
}

sub Consolidation
{
	my $self = shift;

	my $dta = $self->get_dta();
	my $keep = $self->get_keepnum();
##### change on 6/10/2013
 	
	my $dtahash_orig = $self->read_dta_file($dta);
############## 1/7/2014 #################################################	
########### missing out files when searching for a large data set
########### The reason might be too many dta files in the folder, if so, skip the consolidation for this particular scan	
	return 0 if (!defined($dtahash_orig->{'ms2'}));
#########################################################################	
	my %dtahash= %{$dtahash_orig->{'ms2'}};


	my $param = $self->get_parameter();
	my $TMT_data_start_mz = 0;
	my $startmz = (sort {$a <=> $b} keys %dtahash)[0];
	my $prec_mz = $dtahash_orig->{'prec_mz'};
	my $charge = $dtahash_orig->{'prec_charge'};
	
=head
### due to jumpq uses the mzXML file rather dta files so we do not need to keep TMT peaks in dta file	
	if(defined($param->{'TMT_data'}))
	{
		$TMT_data_start_mz = $param->{'TMT_data'};
		
	}
	foreach my $mz (sort {$a <=> $b} keys %dtahash)
	{
		if($mz>$TMT_data_start_mz)
		{
			$startmz = $mz;
			last;
		}
	}
=cut

	if(defined($param->{'TMT_data'}))
	{
###### remove TMT reporter peaks
		my $TMT_peaks = $self->get_TMT_peaks();
		my $TMTc_peaks = $self->get_TMTc_peaks($prec_mz);
		my $tolerance = $self->get_fragment_tolerance($TMT_peaks->[0]);
		my $TMTc_tolerance = $self->get_fragment_tolerance(($TMTc_peaks->[0]));
		foreach my $mz (sort {$a<=>$b} keys %dtahash)
		{
			foreach (@$TMT_peaks)
			{			
				if($mz <($_+$tolerance) and $mz >($_-$tolerance))
				{
					delete $dtahash{$mz};
				}			
			}
			foreach (@$TMTc_peaks)
			{			
				if($mz <($_+$TMTc_tolerance) and $mz >($_-$TMTc_tolerance))
				{
					delete $dtahash{$mz};
				}
			}				
		}		

	}
	


		
=head
	open (DTA, "$dta") || die "can not open the dta file: $dta\n";
	my %dtahash;
	my $line0 = <DTA>;
	my $line1 = <DTA>;
	chomp $line1;
	my ($startmz, $startint) = split(/ /, $line1);
	$dtahash{$startmz}=$startint;
	
	while(<DTA>)
	{
		chomp $_;
#		$_=~s/[[:cntrl:]]+//g;
		my ($mz,$int)=split(/\s+/,$_);
		$dtahash{$mz}=$int;
	}
	close(DTA);
=cut	
	system(qq(rm $dta));
	open (DTA, ">$dta");
	print DTA $prec_mz," ",$charge,"\n";
#	print DTA "$line0";

	$startmz = sprintf("%f", $startmz);
	my $nextmz = $startmz + 100; 
	my %temp;
	my $num = keys(%dtahash);
	my $i = 0;
	my $j=0;
	my $win_num=0;
# %strong_peak_window used for saving the mass value of the strongest peak within each window
	my $top_peak = $self->get_top_peak(\%dtahash);


	
	my %strong_peak_window;
	foreach my $mz (sort {$a<=>$b} keys %dtahash)
	{
###### $j is used to control the last data point	
		$j++;

		if($mz<$startmz)
		{
			print DTA "$mz $dtahash{$mz}\n";
		}
		
		my $int = $dtahash{$mz};
		if ($mz == 0){ print "Oops....error in converting"; }
	## Get the mass value and save it to the temp hash, but the last number will not put the temp hash	
		if($mz<=$nextmz)
		{
			$temp{$mz}=$dtahash{$mz};
		}
	## goes to next window .. 	

		if($mz>$nextmz)
		{
			if((scalar keys %temp)>1)
			{
				my $k=0;
				my %keepers;
				for my $mz (reverse sort {$temp{$a}<=>$temp{$b}}  keys %temp)
				{
					if($k==0)
					{
						$strong_peak_window{$win_num}=$temp{$mz};
					}
					$keepers{$mz} = $temp{$mz} * $top_peak/$strong_peak_window{$win_num}*(1-0.01*$win_num);
					$k++; 
					last if ($k==$keep);
				}
				for my $mz (sort {$a<=>$b} keys %keepers)
				{
					print DTA "$mz $keepers{$mz}\n";
				}
				$win_num++;
			}
			elsif((scalar keys %temp)==1)
			{
				for my $mz (keys %temp)
				{	
				#	print $win_num,"\t",$mz,"\t",$temp{$mz},"\n";				
					$strong_peak_window{$win_num}=$temp{$mz};
					$temp{$mz} =$temp{$mz}* $top_peak/$strong_peak_window{$win_num}*(1-0.01*$win_num);					
					print DTA "$mz $temp{$mz}\n";
				}
				$win_num++;
			}			
			%temp=();
			$temp{$mz}=$dtahash{$mz};			
			$nextmz += 100;				
		}
		if(scalar( keys %dtahash) == $j)
		{
			my $k=0;
			for my $mz (reverse sort {$temp{$a}<=>$temp{$b}}  keys %temp)
			{
				if($k==0)
				{
					$strong_peak_window{$win_num}=$temp{$mz};
				}
				
				$k++;
			}	
			for my $mz (sort {$a<=>$b} keys %temp)
			{
				$temp{$mz} =$temp{$mz}* $top_peak/$strong_peak_window{$win_num}*(1-0.01*$win_num);
				print DTA "$mz $temp{$mz}\n";
			}					
#			for my $mz (keys %temp)
#			{	
#				print $win_num,"\t",$mz,"\t",$temp{$mz},"\n";			
#				$strong_peak_window{$win_num}=$temp{$mz};	
#				print DTA "$mz $temp{$mz}\n";

#			}			
		}		
	}

	close(DTA);


}

sub Consolidation2
{
	my $self = shift;

	my $dta = $self->get_dta();
	my $keep = $self->get_keepnum();
	
}

sub Normalize_Consolidation
{
	my ($self)=@_;
	my $strong_peak_window = $self->get_strong_peak_window();
	my $dtafile = $self->get_dta();
	my $dtahash = $self->read_dta_file($dtafile);

	my $keep = $self->get_keepnum();
######## check if there are any peaks with pho neutral loss
	my $param = $self->get_parameter();
##############################################################		
	my $top_peak = $self->get_top_peak($strong_peak_window);




	my $rank_strong_peak;
	my $j=0;
	foreach (reverse sort {$strong_peak_window->{$a}<=>$strong_peak_window->{$b}} keys %$strong_peak_window)
	{
		$rank_strong_peak->{$strong_peak_window->{$_}}=$j;

		$j++;
	}

	my $i=0;
	foreach my $mz (sort {$a<=>$b} keys %{$dtahash->{'ms2'}})
	{

		my $window_num = int($i/$keep);
		$i++;
		my $strong_peak_within_window = $strong_peak_window->{$window_num};
		next if (!defined($strong_peak_within_window));

		$dtahash->{'ms2'}->{$mz} = $dtahash->{'ms2'}->{$mz}*($top_peak/$strong_peak_within_window)*(1-0.01*$rank_strong_peak->{$strong_peak_within_window});

	
	}

	$self->write_dta_file($dtafile,$dtahash);

}

sub get_top_peak
{
	my ($self,$hash)=@_;
	my $max = 0;
	foreach my $value (keys %$hash ) {
		$max = $hash->{$value} if ($hash->{$value} > $max);
	}
	return $max;
}

sub remove_neutral_losses
{

}

sub read_dta_file
{
	my ($self,$dtafile) = @_;
	
	open (DTA, "$dtafile") || die "can not open the dta file: $dtafile\n";
	my %dtahash;
	my $line0 = <DTA>;
	my @data=split(/\s+/,$line0);

	$dtahash{'prec_mz'} = $data[0];
	$dtahash{'prec_charge'} = $data[1];

	while(<DTA>)
	{
		chomp $_;
		my ($mz,$int)=split(/\s+/,$_);
		$dtahash{'ms2'}{$mz}=$int;
	}
	close(DTA);
	return \%dtahash;
}

sub write_dta_file
{
	my ($self,$dtafile,$dtahash) = @_;

	if(-e $dtafile)
	{
		system(qq(rm $dtafile));
	}
	open (DTA, ">$dtafile") || die "can not open the dta file: $dtafile\n";
	print DTA $dtahash->{'prec_mz'}," ",$dtahash->{'prec_charge'},"\n";

	foreach my $mz (sort {$a<=>$b} keys %{$dtahash->{'ms2'}})
	{
		my $intensity = sprintf("%.6f",$dtahash->{'ms2'}->{$mz});
		$mz = sprintf("%.6f",$mz);
		print DTA $mz," ",$intensity,"\n";
	}
	close(DTA);
}

sub remove_TMTc_peaks
{
	my ($self,$dta_file,$minPeak,$dis,$peakTol) = @_;
	my $Dta=new Spiders::Dta();
	$Dta -> {'_dta_file'} = $dta_file;
	$Dta->JumpExe($minPeak,$dis,$peakTol);
	my $int =  $Dta->{'_int_hash'};
	$dis =  $Dta->{'_dis_hash'};
	my $array = $Dta->{'_index_array'};
	my $re = $Dta->{'_rm_value'};
	my %rm_value = %{$re};
	my %hash_int=%{$int};
	my %hash_dis=%{$dis};
	open OUT,">$dta_file" || die "can not open the dta file: $dta_file";;
	print "======Removed=====\n";
	foreach my $in(@{$array})
	{
		if(!exists $rm_value{$in})
		{
			print OUT "$in\t$hash_int{$in}\n";	
		}
		else
		{
			print "$in\t$hash_int{$in}\n";
		}
	}
	close(OUT);
}

sub get_TMT_peaks
{
	my ($self) = @_;
	my @TMT_peaks = (126.127726,127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145,131.138180);
	return \@TMT_peaks;
}

sub get_TMTc_peaks
{
	my ($self,$prec_MH) = @_;
	my $H = $self->get_H_value();	
	my $TMT_peaks = $self->get_TMT_peaks();
	my @TMTc_peaks;
	foreach (@$TMT_peaks)
	{
		my $TMTc_values = $prec_MH- 27.9949146221 - 131.1673996802 + $H + $_ - $TMT_peaks->[0]; 
		push(@TMTc_peaks,$TMTc_values);
	}
	return \@TMTc_peaks;
}
 
1;
