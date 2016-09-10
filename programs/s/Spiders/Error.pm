#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Error
######### Sequest ###########################################
#                                                           # 
#       **************************************************  #    
#       **** Deisotope program for MS2		          ****  #    
#       ****					                      ****  #    
#       ****Copyright (C) 20212 - Xusheng Wang	      ****  #    
#       ****all rights reserved.		              ****  #    
#       ****xusheng.wang@stjude.org		              ****  #    
#       ****					                      ****  #    
#       ****					                      ****  #    
#       **************************************************  #   
#############################################################

use strict;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

package Spiders::Error;

sub new{
    my ($class) = @_;
    my $self = {
    };
    bless $self, $class;
    return $self;
}

sub no_rawfile_error{
	my $self=shift;
	print "\nERROR: No raw input file(s)\n\n\n";
	$self->usage();
}

sub no_dta_path
{
	my $self=shift;
	print "please specify the dta output path\n";
}

sub no_sim_path
{
	my $self=shift;
	print "please specify the simulation output path\n";
}

sub usage{
	my $self=shift;
	my $progname = $0;
print <<"EOF";
Multiplexed Mass spectrometry program
Usage: $progname <options> rawfile.raw
   -help               print this measly help
   -verbose            lots of messages
EOF
exit 1;
}

1;
