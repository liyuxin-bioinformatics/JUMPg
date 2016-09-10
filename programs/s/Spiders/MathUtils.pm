#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::MathUtils

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

package Spiders::MathUtils;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}

sub min { 
	my ($self,$hash) = @_;

	my $min=1000000;
	foreach my $key (keys %$hash ) {
		$min = $hash->{$key} if $hash->{$key} < $min; 
	}
	return $min;
}

sub max { 
	my ($self,$hash) = @_;
	my $max=0;
	foreach my $key (keys %$hash ) { 
		$max = $hash->{$key} if $hash->{$key} > $max; 
	}
	return $max;
}
=head
sub min_array { 
	my ($self,$array) = @_;

	my $min=$array->[0];
	foreach my $key (keys %$hash ) {
		$min = $hash->{$key} if $hash->{$key} < $min; 
	}
	return $min;
}

sub max_array { 
	my ($self,$array) = @_;
	my $max=0;
	foreach my $key (keys %$hash ) { 
		$max = $hash->{$key} if $hash->{$key} > $max; 
	}
	return $max;
}
=cut
sub combinations
{
    my ($self, $hash, $comb_num) = @_;
	
    my %new_hash;

    for (my $i = 0 ; $i < $comb_num ; $i++)
    {
        if ($i == 0)
        {
            foreach my $key (keys(%$hash))
            {
                $new_hash{$key} = $$hash{$key};
            }
        }
        else
        {
            foreach my $hash_key (keys(%new_hash))
            {
                foreach my $key (keys(%$hash))
                {
                    my $key_join = "$hash_key"."$key";
					my @alpha=sort(split(//,$key_join));
					my $key_comb=join('',@alpha);
                    $new_hash{$key_comb} = $$hash{$key} + $new_hash{$hash_key};
                }
            }
        }
    }
    return (%new_hash);
}



1;

