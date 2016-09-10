#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::ProcessingMzXML

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

package Spiders::ProcessingMzXML;

use strict;
use warnings;
use File::Basename;
use Spiders::XMLParser;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_mzXML_file=>undef,
		_dta_path =>undef,
		_converter=>undef,
    };
    bless $self, $class;
	return $self;
}

sub set_mzXML_file
{
	my ($self,$mzXML_file)=@_;
	return $self->{_mzXML_file}=$mzXML_file;
}

sub get_mzXML_file
{
	my $self=shift;
	return $self->{_mzXML_file};
}

sub get_mzXMLfile_dirname
{
	my $self=shift;
	my $mzXMLfile = $self->get_mzXML_file();

	my @suffixlist=(".mzXML");
	my $mzXMLfile_dirname = dirname($mzXMLfile,@suffixlist);

	return $mzXMLfile_dirname;
}

sub set_dta_path
{
	my ($self,$dtapath) = @_;

	if(!(-e $dtapath))
	{
		system(qq(mkdir $dtapath >/dev/null 2>&1));
	}
	$self->{_dta_path}=$dtapath;
}

sub get_dta_path
{
	my ($self) = shift;
	return $self->{_dta_path};
}


sub set_converter
{
	my ($self,$converter)=@_;
	$self->{_converter}=$converter; 
}

sub get_converter
{
	my $self=shift;
	if(!defined($self->{_converter}))
	{
		$self->{_converter} = "/usr/local/bin/MzXML2Search ";	
	}
	return $self->{_converter};
}

sub mzXML2dta
{
	my $self=shift;
	my $mzXMLfile = $self->get_mzXML_file();
	my $dta_path = $self->get_dta_path();

	my $converter = $self->get_converter();

	print "  Converting the mzXML file to dta files\n";

	system(qq($converter $mzXMLfile >/dev/null 2>&1));

}

####### the following subroutine is used to process mzXML without relying on other program

sub readmzXML
{
	my ($self)=@_;
	my $mzXML = $self->get_mzXML_file();

	open (XML, "<$mzXML") || die "can not open the file";
	my $xml = new Spiders::XMLParser();
	my $indexOffset = $xml->get_IndexOffset(*XML); 
	my ($index_array, $last_scan) = $xml->get_IndexArray(*XML, $indexOffset);
	$self->{'_last_scan_num'} = $last_scan;
	$self->{'_total_scan_num'} = scalar (@$index_array);
	$self->{'_index_array'} = $index_array;
	
	my ($runtime) = $xml->get_RT(*XML, $$index_array[$last_scan-1]);
	$self->{'_total_runtime'} = $runtime;
	return *XML;
}

sub get_last_scan_num
{
	my $self=shift;
	return $self->{'_last_scan_num'};
}

sub get_total_scan_num
{
	my $self=shift;
	return $self->{'_total_scan_num'};
}

sub get_index_array
{
	my $self = shift;
	return $self->{'_index_array'};
}

sub get_total_runtime
{
	my $self=shift;
	return $self->{'_total_runtime'};
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

sub get_ms1_rt
{
	my ($self,$scan)=@_;
	my $mshash = $self->{'_mshash'};

	return $mshash->{'surveyhash'}->{$scan}->{'rt'};
}



sub generate_hash_dta
{
	my ($self, $ms_hash, $msms_hash, $mz_array, $param_hash) = @_;
	
	my $xml = new Spiders::XMLParser();
	my $parameter = $self->get_parameter();
	my $win = $parameter->{'isolation_window'}/2;


	my %no_prec;

	my $XML = $self->readmzXML();
	my $index_array = $self->get_index_array();
	my $last_scan = $self->get_last_scan_num();

	my ($number, $scan_order, $prev_index, $prev_survey) = (0,0,0,0);
	my $found = 0;
	my $cycle = 0;

	foreach (@$index_array)
	{
		my $index = $_;
		my ($deleted);
		$number++;
		my $mslevel=0;


		next if ($number>$parameter->{'last_scan_extraction'});			
		print "\r  Gathering scan information: $number of $last_scan scans          ";
		my ($rt) = $xml->get_RT(*XML, $index);

		($$msms_hash{$number}{'xml_mz'}, $$msms_hash{$number}{'xml_int'}, $$msms_hash{$number}{'xml_act'},$mslevel) = $xml->get_PrecursorMZINTACT(*XML, $index);
        if($mslevel == 2)
        {
			next if ($number<$parameter->{'first_scan_extraction'});		
  #  if ($xml->get_MSLevel(*XML, $index) != 1){ # msms scan
			next if ($prev_survey == 0);
			$$msms_hash{$number}{'survey'} = $prev_survey;

     #$$msms_hash{$number}{'xml_int'} = $xml->get_PrecursorIntensity(*XML, $index);
			my ($mz, $int, $empty) = $self->strongest_prec(*XML, $$msms_hash{$number}, $prev_index, $win);
			next if ($empty == 0);
			if ($mz == 0)
			{
				$$msms_hash{$number}{'unverified'} = 1;
				print "  The precursor ion for scan $number is not present in MS ...\n";
				print Dumper($$msms_hash{$number});
				print Dumper($$msms_hash{$number});
				exit;
			}
			                        # Add in msms values
			my @peaks_array;
			$xml->get_Peaks(*XML, \@peaks_array, $index);
			my $values = scalar(@peaks_array);
			my $peak_num = $values/2;
#			next if ($peak_num < 5);

			## Modified by JCho on 12/22/2014 ######################
			if ($$parameter{'TMT_data'} == 0 && $peak_num < 5) {
				next;
			} elsif ($$parameter{'TMT_data'} == 1) {
				my $numPeaks_TMTdata = 0;
				for (my $i = 0; $i < $values; $i++) {
					$i++;
					if ($peaks_array[$i - 1] > 150) {
						$numPeaks_TMTdata++;
					}
				}
				next if ($numPeaks_TMTdata < 5);
			}
			########################################################

			$$msms_hash{$number}{'prec_int'} = $int;
			$$msms_hash{$number}{'prec_mz'} = sprintf("%.6f", $mz);
			my ($sum, $sum2) = (0,0);
			my @msms_mz;
			my @msms_int;	
			for (my $i=0;$i<$values;$i++)
			{
				
				$i++;
				push(@msms_mz,$peaks_array[$i-1]);
				push(@msms_int,$peaks_array[$i]);
				my $temp = $peaks_array[$i];
				$sum += $temp; 
				$sum2 += $temp**2;
			}
			$$msms_hash{$number}{'msms_mz'} = \@msms_mz;
			$$msms_hash{$number}{'msms_int'} = \@msms_int;
			
			my $total = $sum;
			$$msms_hash{$number}{'peak_num'} = $peak_num;
			$$msms_hash{$number}{'intensity'} = $total;
			$$msms_hash{$number}{'rt'} = $rt;
			
			my $scannum = $number;
		
			#$self->generate_dta_file($scannum,$mz,\@peaks_array);
		}
		elsif($mslevel == 1)
		{ # ms scan
			$self->create_surveyhash(*XML, $param_hash, \%{$$ms_hash{'surveyhash'}{$number}}, $index, $scan_order,\@{$$mz_array[$number]});
			$$ms_hash{'orderhash'}{$scan_order} = $number;
			$$ms_hash{'surveyhash'}{$number}{'rt'} = $rt;
			$prev_survey = $number; $prev_index = $index;   $scan_order++;
		}
		elsif($mslevel == 3)
		{
			my @peaks_array;
			$xml->get_Peaks(*XML, \@peaks_array, $index);
			my $values = scalar(@peaks_array);

##################### generate the dta files ############################
			my $ms3_prec_mz = $$msms_hash{$number}{'xml_mz'};
			my $scannum = $number;
			$self->generate_ms3_file($scannum,$ms3_prec_mz,\@peaks_array);
		}
	}	

	$self->{'_mshash'} = $ms_hash;
	$self->{'_msmshash'} = $msms_hash;
}

sub generate_dta_file
{
	my ($self,$scannum,$mz,$peaks_array) = @_;
	my $H = 1.007276466812;
	my $dta_path = $self->get_dta_path();
	my $parameter = $self->get_parameter();	
	my $prec_window = $parameter->{'prec_window'}/2;
	
	my $mzXML_filename = $self->get_mzXML_file();

	my @suffixlist=(".mzXML");
	my $Rawfile_basename = basename($mzXML_filename,@suffixlist);
	
#### generating dta file #########################
	if(defined ($parameter->{'bypass_precursor_preprocessing'}) and $parameter->{'bypass_precursor_preprocessing'} == 2)
	{
#		my @dta = glob("$curr_dir/$dta_path/*.dta");
		my $charge_assignment = $parameter->{'manual_charge_assignment'};
		$charge_assignment =~ s/\s*//g;
		
		my @charge_array = split(/\,/,$charge_assignment);

		foreach my $charge (@charge_array)
		{
			my $file = "$Rawfile_basename\.$scannum\.$scannum" . ".$charge.dta";
			my $values = scalar(@$peaks_array);
		## if the MS2 values <5 then we will not export the dta files	
			return 0 if ($values<=10);
			open(FILE,">$dta_path/$file") || die "can not open the file\n";

			print FILE ($mz-$H)*$charge+$H," ",$charge,"\n";

			for (my $i=0;$i<$values;$i++)
			{
				$i++;
				next if($$peaks_array[$i-1]<($mz+$prec_window) and $$peaks_array[$i-1]>($mz-$prec_window));
				
				print FILE $$peaks_array[$i-1]," ",$$peaks_array[$i],"\n"; 
			}
			close(FILE);
		}
	}
	else
	{

		my $charge = 2;
		my $file = "$Rawfile_basename\.$scannum\.$scannum" . ".$charge.dta";
		my $values = scalar(@$peaks_array);
	## if the MS2 values <5 then we will not export the dta files	
		return 0 if ($values<=10);
		open(FILE,">$dta_path/$file") || die "can not open the file\n";

		print FILE ($mz-$H)*2+$H," ",2,"\n";

		for (my $i=0;$i<$values;$i++)
		{
			$i++;
			next if($$peaks_array[$i-1]<($mz+$prec_window) and $$peaks_array[$i-1]>($mz-$prec_window));
			
			print FILE $$peaks_array[$i-1]," ",$$peaks_array[$i],"\n"; 
		}
		close(FILE);
	}	
}

sub generate_ms3_file
{
    my ($self, $scannum,$mz,$peaks_array) = @_;
	my $dta_path = $self->get_dta_path();
	my $Rawfile_basename = basename($dta_path);
	$Rawfile_basename =~ s/\..*//;
	#  my $charge = 1;
	my $file = "$Rawfile_basename\.$scannum\.$scannum" . ".ms3";
    my $values = scalar(@$peaks_array);


    open(FILE,">$dta_path/$file") || die "can not open the file\n";
    print FILE $mz," ",1,"\n";

    for (my $i=0;$i<$values;$i++)
    {
            $i++;
            print FILE $$peaks_array[$i-1]," ",$$peaks_array[$i],"\n";
    }
}

sub strongest_prec{
	my ($self, $XML, $hash, $prev, $win) = @_;
	my $xml = new Spiders::XMLParser();
	my ($MZ, $INT) = (0, 0);

	# Find strongest precursor peak within +/- .003
	my $mz = $$hash{'xml_mz'};
	my ($low, $high) = ($mz-$win-.10, $mz+$win+.10); 
	#print "$mz mz: low=$low high=$high";
	my @peaks_array; 
	$xml->get_Peaks(*$XML, \@peaks_array, $prev);
	my $values = scalar(@peaks_array);
	my ($intensity, $sumint) = (0, 0);
	for (my $i=0;$i<$values;$i++)
	{
		my $mz = $peaks_array[$i];
		my $int = $peaks_array[$i+1];
		if ($mz < $low)
		{ 
			$i++;
			next;
		}
		last if ($mz > $high);
		$i++;
		$sumint+=$int;
		if ($int>$intensity)
		{
			$intensity = $int; 
			($MZ, $INT)  = ($peaks_array[$i-1], $intensity);	
		}
	}
	return ($MZ, $INT, $sumint); 
}


sub create_surveyhash{
	my ($self, $XML, $paramhash, $mshash, $index, $scan_order, $mz_array) = @_;
################### xusheng ###############################	

	my $parameter = $self->get_parameter();
	my $xml = new Spiders::XMLParser();
	my $win = $parameter->{'preproc_ms1_consolidate'};

#	my $win = 0.19; 
	if (!defined($win)){$win=0};
	# Get all peaks from survey scan and calculate background
	my ($sum, $sum2) = (0,0);	
	my @peaks_array;
	$xml->get_Peaks($XML, \@peaks_array, $index);
	my $values = scalar(@peaks_array);
	my $peak_num = $values/2;
	my %hash;	my @array;
	for (my $i=0;$i<$values;$i+=2)
	{
		my $origmz = sprintf("%.6f", $peaks_array[$i]);
		my $mz = $origmz; $mz =~ s/\..*//;
		my $int = $peaks_array[$i+1];
		$array[$mz]{$origmz} = $int;
		$hash{$origmz} = $int;
	}
	if ($win != 0)
	{
		my %temp = %hash;
		my $largest = 0;
		for my $origmz (sort {$hash{$b}<=>$hash{$a}} keys %hash)
		{
			next if (!defined($temp{$origmz}));
			my ($low, $high) = ($origmz-$win, $origmz+$win); 
			my ($lowint, $highint) = ($low, $high);
			for my $testmz (keys %{$array[$lowint]})
			{
				next if ($testmz > $high || $testmz < $low);
				next if ($testmz == $origmz);
				delete $temp{$testmz};
				#print "     removed $testmz\n";
			} 
			if ($lowint != $highint)
			{
				for my $testmz (keys %{$array[$highint]})
				{
					next if ($testmz > $high);
					next if ($testmz < $low || $testmz == $origmz);
					delete $temp{$testmz};
					#print "     removed $testmz\n";
				} 
			}
		}
		%hash = %temp;
	}
	while (my ($origmz, $int) = each %hash)
	{
		my $mz = $origmz;
		$mz =~ s/\..*//;
		$$mz_array[$mz]{$origmz} = $int;
		$sum += $int;
		$sum2 += $int**2;
	}
	$peak_num = scalar(keys %hash);

	# Delete scans without peaks
	if ($values <= 0){ return (0);}

	$$mshash{'scan_order'} = $scan_order; 
	$$mshash{'peak_num'} = $peak_num; 
	$$mshash{'tic'} = $sum;
}

1;
