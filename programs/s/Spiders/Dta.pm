#!/usr/bin/perl

## example code
## Release date: 05/01/2015
## Release version: version 12.0.0
## Module name: Dta

package Spiders::Dta;

use strict;
use warnings;
use File::Basename;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
require Exporter;
@ISA	 = qw(Exporter);
@EXPORT      = qw(JumpExe);
 
sub new{
    my ($class,%arg) = @_;
    my $self = {
        _dta_file => undef,
    };
    bless $self, $class;
    return $self;
}
		
sub set_dta_file
{
	my ($self,$dtafile) = @_;
	return $self->{'_dta_file'} = $dtafile;	
}

sub get_dta_file
{
	my $self = shift;
	return $self->{'_dta_file'};
}

sub parse_dtafile_name
{
	my $self=shift;
	my $dtafile = $self->get_dta_file();
	$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
	my ($scan, $specscan, $specscan2, $charge) = ($3, $3, $4, $5); 
	$self->{'_scan'} = $scan;
	$self->{'_charge'} = $charge;
}

sub get_scan_num
{
	my $self = shift;
	return $self->{'_scan'};
}

sub get_charge
{
	my $self = shift;
	return $self->{'_charge'};
}

sub process_dtafile
{
	my $self = shift;
	my $dtafile = $self->get_dta_file();

# define the mz hash for storing the mz and intensity
    my %mz_hash;
# open each dta file to load the data 
    open(DTAFILE,$dtafile) || die "can not open the dta file: $dtafile";
# get the precusor mz and charge 
    my $prec_mz_charge = <DTAFILE>;
	chomp $prec_mz_charge;
    my ($prec_mz,$prec_charge) = split(/\s+/,$prec_mz_charge);
	my @mz_array;
	my @int_array;
    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
		my @data =split(/\s+/,$_);
		$mz_hash{$data[0]} = $data[1];
		push (@mz_array,$data[0]);
		push (@int_array,$data[1]);
	}
	close(DTAFILE);
	
	$self->{'_charge'} = $prec_charge;
	$self->{'_precMZ'} = $prec_mz;
	$self->{'_mz_int_hash'} = \%mz_hash;
	$self->{'_mz_array'}=\@mz_array;
	$self->{'_int_array'} = \@int_array;
	
}

sub write_dta_file
{
	my ($self) = @_;
	my $dtafile = $self->get_dta_file();
	my $mz_int_hash = $self->get_mz_int_hash();
	
	if(-e $dtafile)
	{
		system(qq(rm $dtafile));
	}
	open (DTA, ">$dtafile") || die "can not open the dta file: $dtafile\n";
	
	print DTA $self->get_prec_mz()," ",$self->get_charge(),"\n";

	foreach my $mz (sort {$a<=>$b} keys %{$mz_int_hash})
	{
		my $intensity = sprintf("%.6f",$mz_int_hash->{$mz});
		$mz = sprintf("%.6f",$mz);
		print DTA $mz," ",$intensity,"\n";
	}
	close(DTA);
}

sub set_prec_mz
{
	my ($self,$prec_mz) = @_;
	$self->{'_precMZ'} = $prec_mz;
}

sub set_charge
{
	my ($self,$charge) = @_;
	$self->{'_charge'} = $charge;
}

sub set_mz_int_hash
{
	my ($self,$mz_int_hash) = @_;
	$self->{'_mz_int_hash'} = $mz_int_hash;
}

sub get_prec_mz
{
	my $self = shift;
	return $self->{'_precMZ'};
}

sub get_mz_int_hash
{
	my $self = shift;
	return $self->{'_mz_int_hash'};
}

sub get_mz_array
{
	my $self = shift;
	return $self->{'_mz_array'};
}

sub get_int_array
{
	my $self = shift;
	return $self->{'_int_array'};
}
sub JumpExe{

        my $self = shift;
	my $peak = shift;
	my $DIS = shift;
	my $tol = shift;
        my $dtafile = $self -> get_dta_file();
        my %inthash;
	my %dishash;
   open(JumpFILE,$dtafile) || die "can not open the dta file: $dtafile";
   chomp(my @array=<JumpFILE>);
	my @rearray = reverse @array;
	my @index;
	my @tmp = split/\s+/,$array[0];
	my $precursor_mass = $tmp[0];
	my $TMTc = $precursor_mass - 154.1226406;
	my $TMTc_mass_shift;
	my $precursor_mass_cutoff = $tmp[0]/$tmp[1];
	$inthash{$precursor_mass}=$tmp[1];
	foreach my $i(0..$#rearray-1)
    {
	    my @data = split/\s+/,$rearray[$i];
		my @pre_data = split/\s+/,$rearray[$i+1];
		my $distance;
		push(@index,$data[0]);
		if($i == $#rearray-1 ){
			#the father node
			$distance = 0;
		}
		else{
			$distance = $data[0]-$pre_data[0];
		}
		if($i == 0){$dishash{$data[0]}=0};
        	$inthash{$data[0]} = $data[1];
			$dishash{$pre_data[0]} = $distance;
	}
	push(@index,$precursor_mass);
    close JumpFILE;
	
	#### Execute
	### split two distance parameters:
	my @disval=split/,/,$DIS;
	my $Tol_dis2 = $disval[1]+$tol;
	my $Tol_dis2_pre = $disval[1]-$tol;
	my $Tol_dis1 = $disval[0]+$tol;
	my @ok_array;
	my $count=0;
	my @Tmp_sec_array;
	for(my $j=0;$j<=$#index;$j++)
	{	
		if($dishash{$index[$j]}<=$Tol_dis2 && $index[$j]>$precursor_mass_cutoff)
		{
			if($dishash{$index[$j]}>=$Tol_dis2_pre)
			{
				$count+=1;
				if($count==1){push(@Tmp_sec_array,$index[$j-1]);$count+=1;}
				push(@Tmp_sec_array,$index[$j]);
			}	
		}
		if($j == 0){push(@Tmp_sec_array,$index[$j]);$count+=1;next;}
		if($dishash{$index[$j]}>$Tol_dis2 || $dishash{$index[$j]}< $Tol_dis2_pre || $index[$j] <=$precursor_mass_cutoff)
		{
			if($count>=$peak)
			{
				###two
				if($dishash{$index[$j-1]}<=$Tol_dis1 && $dishash{$index[$j-1]}>=$disval[0])
				{
				push(@ok_array,@Tmp_sec_array);
				@Tmp_sec_array=();
				$count = 0;
				}
				else
				{	
					pop(@Tmp_sec_array);
					my $raw_count = $count;
					$count=$count-1;
					my $n=$#Tmp_sec_array;
					foreach my $retmp($n..0)
					{
						if($dishash{$Tmp_sec_array[$retmp]}<=$Tol_dis1 && $dishash{$Tmp_sec_array[$retmp]}>=$disval[0]){last;}
						else{$count=$count-1;pop(@Tmp_sec_array);}
					}
					if($count>=$peak)
					{
						push(@ok_array,@Tmp_sec_array);
						@Tmp_sec_array=();
						my $back_time = $raw_count - $count;
						$count = 0;
						$j=$j-$back_time;
					}
				}
				
			}
			else
			{
				@Tmp_sec_array = ();
				$count = 0;	
			}
			if($index[$j]<=$precursor_mass_cutoff){last;}
			next;
		}
	}
	#shift @ok_array;
	my %consecutive_VALUE;
	foreach(@ok_array)
	{
		$consecutive_VALUE{$_}=1;
		if(!defined($TMTc_mass_shift))
		{
			my $temp_shift = ($TMTc - $_)/$TMTc*1000000;
			if(abs($temp_shift)<50)
			{
				$TMTc_mass_shift = $temp_shift;
				print "masscorrection: ",$dtafile,"\t",$precursor_mass,"\t",$TMTc_mass_shift,"\n";
			}
		}
	}
    my @reindex = reverse @index;
	$self->{'_int_hash'}=\%inthash;
	$self->{'_dis_hash'}=\%dishash;
	$self->{'_index_array'}=\@reindex;
	$self->{'_rm_value'}=\%consecutive_VALUE;
	
}
1;

