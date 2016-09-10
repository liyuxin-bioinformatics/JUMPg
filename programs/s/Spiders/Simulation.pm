#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Simulation

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
package Spiders::Simulation;


use strict;
use warnings;
use File::Basename;
use Spiders::Dta;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
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

sub set_dta_file
{
	my ($self,$dta)=@_;
	$self->{'_dta_file'}=$dta;
}

sub set_param
{
	my ($self,$params)=@_;
	$self->{'_param'}=$params;
}

sub get_dta_file
{
	my ($self)=@_;
	return $self->{'_dta_file'};
}

sub get_param
{
	my ($self)=@_;
	return $self->{'_param'};
}

sub set_sim_path
{
	my ($self,$dir)=@_;
	if(!(-e $dir))
	{
		system(qq(mkdir $dir >/dev/null 2>&1));
	}
	$self->{_sim_path}=$dir;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{_dta_path};
}

sub get_sim_path
{
	my $self=shift;
	return $self->{_sim_path};
}


sub Simulation_Multiplexing_Two_DataSets
{
	#####################
	# Merge files in two different directories using the same file name
	#
	###################
	my $self=shift;

	my $Directory_A = $ARGV[0];
	my $Directory_B = $ARGV[1];

	my @file_A=glob("/home/xwang4/Multiplexing_data/ratbrain01ug_muliplex2a/$Directory_A/$Directory_A.1/*.dta");
	foreach my $file (@file_A)
	{
		print $file,"\n";
		my $fileout=basename($file);
		my $dtafile=$fileout;
		$dtafile =~ s/([A-Za-z0-9\_\-]+)\.((\d+)\.(\d+))\.(\d+).dta\Z/$2/;
		my $file1=glob("/home/xwang4/Multiplexing_data/ratbrain01ug_muliplex2a/$Directory_B/$Directory_B.1/*.$dtafile.*.dta");
		next if (!(defined $file1));
		my $fileout1=basename($file1);
		copy($file, "/home/xwang4/Multiplexing_data/ratbrain01ug_muliplex2a/Orig_Multiplex/$fileout") or die "Copy Faild: $!";
		copy($file1, "/home/xwang4/Multiplexing_data/ratbrain01ug_muliplex2a/Orig_Multiplex1/$fileout1") or die "Copy Faild: $!";
		my %hash_intensity;
		open(FILE,"$file") || die "can not open the file: $file\n";
		my $header=<FILE>;
		while(<FILE>)
		{
			my @data=split(/\s+/,$_);
			$hash_intensity{$data[0]}=$data[1];
		}
		open(FILE1,"$file1") || die "can not open the file: $file\n";
		my $header2=<FILE1>;
		while(<FILE1>)
		{
			my @data=split(/\s+/,$_);
			$hash_intensity{$data[0]}=$data[1];
		}
		close(FILE);
		close(FILE1);
		open(OUTPUT,">Sim_Multiplex/$fileout");
		print OUTPUT $header;
		foreach my $mz (sort {${a}<=>${b}} keys %hash_intensity)
		{
			print OUTPUT $mz," ",$hash_intensity{$mz},"\n";
		}	
		open(OUTPUT1,">Sim_Multiplex1/$fileout1");
		print OUTPUT1 $header2;
		foreach my $mz (sort {${a}<=>${b}} keys %hash_intensity)
		{
				print OUTPUT1 $mz," ",$hash_intensity{$mz},"\n";
		}
	}
}

sub sort_dtafile_by_scan
{
	my ($self,$dta_file) = @_;

	my %dtafile_hash;
	my @sorted_file;
	foreach my $file (@$dta_file) {
		if($file =~ /([A-Za-z0-9\_\-]+)\.((\d+)\.(\d+))\.(\d+).dta\Z/)
		{
			$dtafile_hash{$3}=$file;
		}	
	}

	foreach my $scan (sort {$a<=>$b} keys %dtafile_hash) {
		push (@sorted_file,$dtafile_hash{$scan});
	}

	return \@sorted_file;
}

sub Simulation_Multiplexing_One_DataSet
{
	my $self=shift;
	my $dta_path = $self->get_dta_path();
	my $simulated_path = $self->get_sim_path();
	my $select_num = 60;

	my @file=glob("$dta_path/*.dta");


	my $sorted_file_ref = $self->sort_dtafile_by_scan(\@file);

	my @sorted_file = @$sorted_file_ref;

	for(my $i=0; $i<$#sorted_file;$i++)
	{
		my $fileout = $sorted_file[$i];
		my $fileout1 = $sorted_file[$i+1];
		print $fileout,"\t",$fileout1,"\n";

		$fileout = fileparse($fileout);
		$fileout1 = fileparse($fileout1);
#		$dtafile =~ s/([A-Za-z0-9\_\-]+)\.((\d+)\.(\d+))\.(\d+).dta\Z/$2/;

		my %hash_intensity;
		open(FILE,"$sorted_file[$i]") || die "can not open the file: $sorted_file[$i]\n";
		my $header=<FILE>;
		while(<FILE>)
		{
			my @data=split(/\s+/,$_);
			$hash_intensity{$data[0]}=$data[1];
		}

### normalize the intensity in order to combine the data in the same order of magnitude

		$self->normalize_intensity(\%hash_intensity,$select_num);
# read the second file
		my %hash_intensity1;
		open(FILE1,"$sorted_file[$i+1]") || die "can not open the file: $sorted_file[$i+1]\n";
		my $header2=<FILE1>;
		while(<FILE1>)
		{
			my @data=split(/\s+/,$_);
			$hash_intensity1{$data[0]}=$data[1];
		}
		$self->normalize_intensity(\%hash_intensity1,$select_num);
# merge two hashes
		@hash_intensity{ keys %hash_intensity1 } = values %hash_intensity1;

		close(FILE);
		close(FILE1);

		open(OUTPUT,">$simulated_path/$fileout") || die "can not open the file for output: $simulated_path/$fileout";
		print OUTPUT $header;
		foreach my $mz (sort {$a<=>$b} keys %hash_intensity)
		{
			print OUTPUT $mz," ",$hash_intensity{$mz},"\n";
		}	
		open(OUTPUT1,">$simulated_path/$fileout1") || die "can not open the file for output: $simulated_path/$fileout";
		print OUTPUT1 $header2;
		foreach my $mz (sort {$a<=>$b} keys %hash_intensity)
		{
			print OUTPUT1 $mz," ",$hash_intensity{$mz},"\n";
		}
		$i++;
	}
}

sub normalize_intensity
{
	my ($self,$mz_hash,$select_num) = @_;

	my $mean = 0;
	my $i=0;
	my $sum=0;
	foreach my $mz (sort {$$mz_hash{$a}<=>$$mz_hash{$b}} keys %$mz_hash)
	{
		$i++;
		if($i<=$select_num)
		{
			$sum += $mz_hash->{$mz};
		}
	}
	my $average = $sum/$select_num;

	foreach my $mz (keys %$mz_hash)
	{
		$mz_hash->{$mz} = $mz_hash->{$mz}/$average;
	}	
}

sub simulation
{
	my ($self) = @_;
	my $dtafile = $self->get_dta_file();
	my $param = $self->get_param();
	my $dta = new Spiders::Dta();
	$dta->set_dta_file($dtafile);
	$dta->process_dtafile();
	my $prec_mz = $dta->get_prec_mz();
	my $mz_int_hash = $dta->get_mz_int_hash();
	if($param->{'sim_MS1'} != 0)
	{
#		my $shift_prec_mz = $prec_mz + $param->{'sim_MS1'} * $prec_mz / 1000000;
		my $shift_prec_mz = $prec_mz + $param->{'sim_MS1'};
		$dta->set_prec_mz($shift_prec_mz);
	}

	if($param->{'sim_MS2'} != 0)
	{
		my $sim_mz_int_hash;
		foreach my $mz (keys %{$mz_int_hash})
		{
			my $shift_mz = $mz - $param->{'sim_MS2'} / 2 + rand($param->{'sim_MS2'});
			$sim_mz_int_hash->{$shift_mz} = $mz_int_hash->{$mz}; 
		}
		$dta->set_mz_int_hash($sim_mz_int_hash);
	}
	
	if($param->{'sim_MS1'} != 0 || $param->{'sim_MS2'} != 0)
	{
		$dta->write_dta_file();
	}
}
1;
