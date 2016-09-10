#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::MassAccuracy
######### MassAccuracy  #######################################
#                                                             #
#       **************************************************    #  
#       **** Binomial module for Tag		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2013 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::MassAccuracy;
        
use strict;
use File::Basename;
use Storable;

use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

=head program example

use File::Basename;
use Storable;
use Spiders::MassAccuracy;

my $path = $ARGV[0];
my $hash_dir = dirname($path)."/.hashes/";
my $ms2_hash_ref = retrieve($hash_dir."./origmsms_hash");

my $massaccu = new Spiders::MassAccuracy();
$massaccu->setmsmshash($ms2_hash_ref);
$massaccu->setmasstolerance(15);
$massaccu->setmasstolerance_frac(0.4);
$massaccu->calculate_dynamic_tolerance();

=cut

sub new{
	my ($class, %arg) = @_;
	my $self = {};
	bless $self, $class;
	return $self
}

sub setmshash
{
	my ($self,$mshash)=@_;
	$self->{'_mshash'}=$mshash;
}

sub getmshash
{
	my ($self)=@_;
	return $self->{'_mshash'};
}

sub setmzarray
{
	my ($self,$mzarray)=@_;
	$self->{'_mzarray'}=$mzarray;
}

sub getmzarray
{
	my ($self)=@_;
	return $self->{'_mzarray'};
}

sub setmsmshash
{
	my ($self,$msmshash)=@_;
	$self->{'_msmshash'}=$msmshash;
}

sub getmsmshash
{
	my ($self)=@_;
	return $self->{'_msmshash'};
}

sub setmasstolerance
{
	my ($self,$tolerance)=@_;
	$self->{'_tolerance'}=$tolerance;
} 

sub getmasstolerance
{
	my ($self)=@_;
	return $self->{'_tolerance'};
} 

sub setmasstolerance_frac
{
	my ($self,$tolerance_frac)=@_;
	$self->{'_tolerance_frac'}=$tolerance_frac;
} 

sub getmasstolerance_frac
{
	my ($self) = @_;
	return $self->{'_tolerance_frac'};
} 


sub calculate_dynamic_tolerance
{
	my $self = shift;
	my $dynamics_tolerance;	
	my $start_mass_tolerance = $self->getmasstolerance();
	my $percentage = $self->getmasstolerance_frac();	
	my $ms2_hash_ref = $self->getmsmshash();

	########## define the dynamic mass tolerance ###################
	my $end_mass_tolerance =  $percentage * $start_mass_tolerance;
	my $reduction_per_window = ($start_mass_tolerance - $end_mass_tolerance) / 500;
	my %ms2 = %{$ms2_hash_ref};

	
	foreach my $scan (keys %ms2)
	{
		if(!defined($ms2{$scan}{'prec_int'}))
		{
			delete $ms2{$scan};
		}
	}	
	my $total_number = scalar keys (%ms2);
	my $moving_window_size = int ($total_number/500)+1;
#	print $total_number,"\n";
	my $window=0;
	my $i=0;
	foreach my $scan (sort {$ms2{$a}{'prec_int'} <=> $ms2{$b}{'prec_int'}} keys %ms2)
	{
		$i++;		
		if(($i % $moving_window_size) == 0)
		{
			$dynamics_tolerance->{$scan} = $start_mass_tolerance - $reduction_per_window * $window;
	#		print $i,"\t",$moving_window_size,"\t",$window,"\t",$dynamics_tolerance->{$scan},"\n";
			$window++;			
		}
		else
		{
			$dynamics_tolerance->{$scan} = $start_mass_tolerance - $reduction_per_window * $window;		
		}
	}	
	return $dynamics_tolerance;
}

1;
