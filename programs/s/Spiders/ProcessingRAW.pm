#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::ProcessingRAW

######### Simulation ##########################################
#                                                             #
#       **************************************************    #  
#       **** Coverter RAW file to MzXML		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::ProcessingRAW;

use strict;
use warnings;
use File::Basename;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();
 
sub new{
	my ($class,%arg)=@_;
    my $self = {
        _raw_file => undef,
		_mzXML_file=>undef,
		_converter=> undef,
    };
    bless $self, $class;
	return $self;
}

sub set_raw_file
{
	my ($self,$raw_file)=@_;
	$self->{_raw_file}=$raw_file;
}

sub get_raw_file
{
	my $self=shift;
	return $self->{_raw_file};
}

sub get_rawfile_basename
{
	my $self=shift;
	my $Rawfile = $self->get_raw_file();

	my @suffixlist=(".raw",".RAW");
	my $Rawfile_basename = basename($Rawfile,@suffixlist);

	return $Rawfile_basename;
}


sub get_rawfile_dirname
{
	my $self=shift;
	my $Rawfile = $self->get_raw_file();

	my @suffixlist=(".raw",".RAW");
	my $Rawfile_dirname = dirname($Rawfile,@suffixlist);

	return $Rawfile_dirname;
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
		$self->{_converter} = "wine /usr/local/bin/ReAdW.exe --mzXML -c";	
	}
	return $self->{_converter};
}

sub raw2mzXML
{
	my $self=shift;
	my $Rawfile = $self->get_raw_file();
	my $Dirname=$self->get_rawfile_dirname();
	my $Rawfile_basename = $self->get_rawfile_basename();

	my $converter = $self->get_converter();

	my $mzXML="$Dirname/$Rawfile_basename.mzXML";

	if(-e $mzXML)
	{
		print "  MZXML file ($mzXML) has already existed\n";
	}
	else
	{
		print "  Converting the RAW file to mzXML format\n";
		system(qq($converter "$Rawfile" $mzXML >/dev/null 2>&1));
	}
	$self->{_mzXML_file}=$mzXML;
	return $mzXML;
}

sub get_mzXML_file
{
	my $self=shift;
	return $self->{_mzXML_file};
}

1;
