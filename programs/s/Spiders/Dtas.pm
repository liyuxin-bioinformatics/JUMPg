#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Dtas

package Spiders::Dtas;


use strict;
use warnings;

sub new
{
	my ($class) = @_;
	my $self = {};
	bless ($self,$class);
	return $self;
}

#-------------------------------------------------------------------------------------------
1;

sub add_dta
{
	my ($self,$dtafile)=@_;

	open(IN,"$dtafile") || die "Cannot open dta file: $dtafile!!!\n";
	
	my $head=<IN>;
	my @t=split(/\s/,$head);
	$self->{$dtafile}->{'mass'}=$t[0];
	$self->{$dtafile}->{'charge'}=$t[1];

	my @int;my $k=0;
	while(<IN>)
	{
		$k++;
		@t=split(/\s/,$_);
		$self->{$dtafile}->{'peak'}->{$k}->{'mz'}=$t[0];
		$self->{$dtafile}->{'peak'}->{$k}->{'int'}=$t[1];
	}
	$self->{$dtafile}->{'peaknum'}=$k;

	close IN;
}

sub print_dtas
{
	my ($self,$dtasfile)=@_;

	open(OUT,">>$dtasfile");

	foreach my $dtafile (keys %{$self})
	{
		my @t=split(/\//,$dtafile); 
		print OUT "$t[$#t] ",$self->{$dtafile}->{'mass'}," ",$self->{$dtafile}->{'charge'},"\n";
		for (my $i=1; $i<=$self->{$dtafile}->{'peaknum'}; $i++) {print OUT $self->{$dtafile}->{'peak'}->{$i}->{'mz'}," ";} print OUT "\n";
		for (my $i=1; $i<=$self->{$dtafile}->{'peaknum'}; $i++) {print OUT $self->{$dtafile}->{'peak'}->{$i}->{'int'}," ";} print OUT "\n";
	}

	close OUT;
}
