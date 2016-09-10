#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Sequest

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

package Spiders::Sequest;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_dta_path=>undef,
		_sequest=>undef,
    };
    bless $self, $class;
	return $self;
}

sub set_sequest_path
{
	my ($self,$sequest_path)=@_;
	$self->{_sequest_path}=$sequest_path;
}

sub get_sequest_path
{
	my $self=shift;
	return $self->{_sequest_path};
}

sub copy_sequest_parameter
{
	my ($self,$parameter_file)=@_;
	my $dta_path = $self->get_dta_path();
	system(qq(cp $parameter_file $dta_path));
}

sub set_sequest_engine
{
	my ($self,$engine)=@_;
	$self->{_sequest}=$engine;
}

sub get_sequest_engine
{
	my $self=shift;
	my $sequest_path = $self->get_sequest_path();	
	if(!defined($self->{_sequest}))
	{
		$self->{_sequest}="$sequest_path/sequest28single";
	}
	return $self->{_sequest};
}

sub Run_Sequest
{
	my $self=shift;
	my $dta_dir = $self->get_dta_path();

	print "Running Sequest\n";
    opendir( DIR, "$dta_dir" );
    my @array = ();
    @array = grep( /\.dta/, readdir(DIR) );
    my $job_num = int( $#array / 100 ) + 1;

    for ( my $i = 0 ; $i < $job_num ; $i++ ) {
        open( JOB, ">$dta_dir/job_$i.sh" )
          || die "can not open the job files\n";
        my $dta_file_temp = "";
        for ( my $j = 0 ; $j < 100 ; $j++ ) {
            if ( ( $i * 100 + $j ) <= $#array ) {
                $dta_file_temp .= " $array[$i*100+$j]";
            }
        }

        print JOB "#BSUB -P prot\n";
        print JOB "#BSUB -q normal\n";
        print JOB "#BSUB -eo ./$i.e\n";
        print JOB "#BSUB -oo ./$i.o\n";
        print JOB "$dta_dir/sequest28single $dta_file_temp\n";
    }

    my $job_list;
    for ( my $i = 0 ; $i < $job_num ; $i++ ) {
        my $command_line = "cd $dta_dir && bsub <job_" . "$i.sh ";
        my $job = qx[$command_line];
        chomp $job;
        my $job_id = 0;
        if ( $job =~ /Job \<(\d*)\> is/ ) {
            $job_id = $1;
        }
        $job_list->{$job} = 1;
    }

    my $job_info = 1;
    my ($username) = getpwuid($<);

    while ($job_info) {
        sleep(10);
        my $command_line = "bjobs -u $username";
        my $job_status   = qx{$command_line 2>&1};

        my @job_status_array = split( /\n/, $job_status );

        if ( $job_status =~ /No unfinished job found/ ) {
            $job_info = 0;
        }
        elsif ( ( scalar(@job_status_array) ) == 0 ) {
            $job_info = 0;
        }
        elsif ( $job_status_array[1] =~ /PEND/ ) {
            print "\r cluster is busy, please be patient!";
            sleep(5);
        }
        else {
            for ( my $i = 0 ; $i < $#job_status_array ; $i++ ) {
                my @data_array = split( /\s/, $job_status_array[$i] );
                next if ( $data_array[0] =~ /\d+/ );
                if ( defined( $job_list->{ $data_array[0] } ) ) {
                    $job_info = 0;
                }
            }
            sleep(1);
        }
    }
    print "Sequest search is done!!\n";
}

sub Run_Sequest_standalone
{
	my ($self,$dta_file)=@_;
	my $sequest=$self->get_sequest_engine();
    system(qq($sequest $dta_file));
}

sub check_results()
{
	my $self=shift;
	my $Result_dir = $self->get_dta_path();

	print "Checking search results:\n";
	opendir(DIR,$Result_dir); 
	my @out = grep {/\.out\Z/} readdir(DIR);
	seekdir (DIR,0); my @dta = grep {/\.dta\Z/} readdir(DIR);
	    closedir(DIR);
	if (scalar(@dta) != scalar(@out))
	{
		printf "     THIS SEARCH IS MIGHT BE INCOMPLETE(%d dta files and %d out files)!!!!!\n", scalar(@dta), scalar(@out);
		if(scalar(@out)/scalar(@dta)<0.9)
		{
			print "YOUR SEARCH IS NEEDED TO RERUN AND PLEASE CONTACT BIOINFORMATICS GROUP!!\n";
		}
 	}
	elsif(scalar(@dta) == 0 || scalar(@out) == 0)
	{
		printf "     THERE ARE NO OUTPUT FILES FOR THIS RUN !!!!\n";
	}
	else
	{
		printf "     This search has finished (%d dta files and %d out files)\n", scalar(@dta), scalar(@out);
	}

	print "\n\n";
}

1;
