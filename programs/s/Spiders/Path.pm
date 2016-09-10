#!/usr/bin/perl -w

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Path

use strict;
package Spiders::Path;

#The purpose of this class is the encapsulise a sequest job, either about
# to be run or already running, and provide methods to control the process.
#It can be instantiated with either a .raw file path or by specifying a pid.

use File::Copy;
use Cwd 'abs_path';


#################################################################
#CLASS METHODS
#################################################################
sub new {
  #dat path needs to be the full path
  my($class,$initpath) = @_;
  
  return 0 if(!defined($initpath));
  my $self = {};
  bless $self;
  
  chomp($self->{PWD} = `pwd`);

  if(my $paths = localize_dat_($initpath))
  {    
      $self->{DATPATH} = $paths->{datpath};
      $self->{BASEDIR} = abs_path($paths->{basedir}) || 0;
      $self->{SUBDIRS} = $paths->{subdirs} || 0;
      $self->{BASENAME} = $paths->{basename};
      
      return $self if($self->{DATPATH});
  }
  else
  {
	return undef;
  }
}


sub make_directory_list {
  my ($self) = @_;

  #find the longest word so that we can do the formatting right below
  my $max = 0;
  foreach (keys %{$self->subdirs()}) {
    my $length = length($_);
    $max = ($max > $length) ? $max : $length;
  }

  #calculate the percent done and make labels for the dir. selection
  #list.
  my @keys = map {$_->[1]} sort {$a->[0] <=> $b->[0]}
  map { [($_ =~ /(\d+)$/o) ? $1 : 0,$_]} keys %{$self->subdirs()};

   return \@keys;
   
}


#################################################################
#PUBLIC INSTANCE METHODS
#################################################################

#Accessor methods


sub currentdir {
  my($self) = shift @_;
  my $return = $self->{CURRENTDIR};
  $self->{CURRENTDIR} = shift if(@_);
  return $return;
}

sub datpath {
  my($self) = shift @_;
  my $return = $self->{DATPATH};
  $self->{DATPATH} = shift if(@_);
  return $return;
}

sub basedir {
  my($self) = shift @_;
  my $return = $self->{BASEDIR};
  $self->{BASEDIR} = shift if(@_);
  return $return;
}

sub subdirs {
  my($self) = shift @_;
  my $return = $self->{SUBDIRS};
  $self->{SUBDIRS} = shift if(@_);
  return $return;
}

sub basename {
  my($self) = shift @_;
  my $return = $self->{BASENAME};
  $self->{BASENAME} = shift if(@_);
  return $return;
}


sub sequest {
  my($self) = shift @_;
  my $return = $self->{SEQUEST};
  $self->{SEQUEST} = shift if(@_);
  return $return;
}

sub extractor {
  my($self) = shift @_;
  my $return = $self->{EXTRACTOR};
  $self->{EXTRACTOR} = shift if(@_);
  return $return;
}



#takes a datpath and constructs a basedir and moves the datfile into it.
#also responsible for adjusting all of the instance data pointing to the 
#dat file, etc.
sub construct_basedir 
{
  my ($self) = @_;
  return if($self->{BASEDIR});
  
  #if no basedir already exists, make make
  my $basedir = construct_basedir_($self->{DATPATH});
  $self->{BASEDIR} = $basedir;
  
  if($self->{DATPATH} =~ /\/?([^\/]+\.dat)$/){
    $self->{DATPATH} = "$basedir/$1";
  }
    
  $self->{SUBDIRS} = {};
  return $basedir;
}

#adds a subdirectory for a given datfile with the suffix passed in
sub add_subdir 
{
  my ($self,$newdir) = @_;
  
  #self->{subdirs} will be defined but empty if there is a valid base
  return if(!$self->{SUBDIRS} ||!$self->{BASEDIR});
    
  #the subdir returned is without the basedir localization
  my $subdir = add_subdir_("$self->{BASEDIR}/$newdir");
  $self->{SUBDIRS}->{$subdir} = $subdir;
  $self->{CURRENTDIR} = $subdir;
  return $subdir;
}

#-----------------------------------------------------------------

sub get_next_subdir {
  my ($self,$dbase)= @_;
  return undef unless(defined($self->{SUBDIRS}) && 
		      defined($self->{BASEDIR}) &&
		      defined($self->{BASENAME}));

  #find the next number in the series.
  #the convention now for naming is: datprefix.database.num
  my @numbered_dirs = grep {/\d+$/} keys %{$self->{SUBDIRS}};
  my @nums = sort {$b <=> $a} map {/\.(\d+)$/;$1 || 0;} @numbered_dirs;
  my $number = (@nums > 0) ? shift(@nums) + 1 : 1;
  
#  return "$self->{BASENAME}.$dbase.$number";
  return "$self->{BASENAME}.$number";
}

#---------------------------------------------------------------------------


#-----------------------------------------------------------------

#################################################################
#PRIVATE INSTANCE METHODS
#################################################################

 
#if the dat file has a base directory, return a list of matching subdirs
#if no base directory exists, return 0
sub get_updir_ {
  chomp(my $pwd = `pwd`);
  if($pwd =~ /([^\/]*)$/o){return $1;}
  return '';
}
     
#This method doesn't try to create anything, just returns the current status
#of things. All of this is rolled into essentially one method to reduce
#redundant parsing of the datfile path and system calls to pwd at the cost
#of reducing some of the generality of finding these thigns out.
sub localize_dat_ {
  my ($path) = @_;
  my %paths = (datpath => "",basedir => "", subdirs => "");
  
  #strip off the trailing / and break up the path into a list
  $path =~ s/\/$//o;
  my @parray = split('/',$path);
  
  #first step is to resolve the given path into the path of the datfile
  my $l = $parray[-1];
  my ($datpath,$dir);

  if(-d $path){
    my $t = "$path/" . $l . ".raw";
    my $t1 = "$path/" . $l . ".RAW";
	my $t2 = "$path/" . $l . ".mzXML";
    if(-e $t)
	{
      $datpath=$t;
      $dir = $path;
    }
    elsif(-e $t1){
      $datpath=$t1;
      $dir = $path;
    }
    elsif(-e $t2){
      $datpath=$t2;
      $dir = $path;
    }	
    else {
      return undef;
    }
  }
  $paths{datpath} = $datpath;
  
  #once we've found the dat file, we can determine the basename
  (my $basename = $datpath) =~ s/^.+\///o;
  $basename =~ s/\.raw$//o;
  $basename =~ s/\.RAW$//o;
  $basename =~ s/\.mzXML$//o;
  my $basedir;
  my $updir = get_updir_();
  
  if(defined($dir) && $dir =~ /$basename$/){
    $basedir = $dir;
  }
  elsif($basename eq $updir){
    $basedir = "../$updir";
  }
  
  $paths{basedir} = $basedir;
  $paths{basename} = $basename;
  $paths{subdirs} = ($basedir) ? resolve_subdirs_($basedir,$basename) : 0;
  return \%paths;
}

#---------------------------------------------------------------------------

#simple method for finding the subdirectories
sub resolve_subdirs_ {

  my($basedir,$prefix) = @_;

  my %subdirs;
  if(opendir(DIR,$basedir))
  {
	map {$subdirs{$_} = "$basedir/$_"} grep {($_=~/$prefix\..+/) && (-d "$basedir/$_")} readdir(DIR);    
	
    close DIR;
  }
  
  return \%subdirs;
}

#-----------------------------------------------------------------
#This is a simple method to pull off the basename of the dat file
#make the corresponding directory and move the dat file into it.
sub construct_basedir_ {
  #this assumes we've checked, and we're not already using subdirs
  #path is the path to the dat file.
  my($datpath) = @_;
  
  my $basedir = $datpath;
  if($basedir =~ s/\.raw//o){
print "########################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##################################\n";
    system(qq(mkdir $basedir >/dev/null 2>&1));
    move($datpath,$basedir);
    return $basedir;
  }
  
  return;
}

#-----------------------------------------------------------------
sub add_subdir_ {
  my($newdir) = @_;
  return 0 unless(defined($newdir));
  my $subdir;
 ##################added by xusheng #####################################
 # because I did not found create a fold for input file
  if(! -e $newdir)
  {
	system(qq(mkdir $newdir  >/dev/null 2>&1));
  }
  else
  {
    print "  Do you wish to overwrite (y or n)?:";
	my $overwrite=<STDIN>;
	chomp $overwrite;
	if($overwrite eq 'y')
	{
      	system(qq(rm $newdir >/dev/null 2>&1));
		system(qq(mkdir $newdir  >/dev/null 2>&1));
    }
	else
	{
		exit;
	}
  }
#  move($datpath,$basedir);
#     ####################################################################
  $subdir->{$newdir} =$newdir;
  return $subdir;
}

#-----------------------------------------------------------------

sub ask 
{
  my ($class,$message,$default) = @_;
  local $| = 1;
  
  my ($result,$choice);
  print "$message [$default]: ";
  chomp($choice = <STDIN>);
  if($choice eq ''){
    return $default;
  }
  else {print '\n';}
  return $choice;
}

sub choose_dir_entry {
  my($class,$dir_entries,$msg,$default) = @_;
  local $|=1;
  my $choice;

  while(1){

    my $num = 1;
    my @entries = @$dir_entries;
	my $updated_file;
	if($entries[$#entries] =~ /(.*)\.(\d+)/)
	{
		my $num = $2 + 1;
		$updated_file = $1 . ".$num"; 
	}
#	print "$msg [$entries[$#entries]]:";
	print "$msg [$updated_file]:";
	
    my $key = <STDIN>;
    chomp($key);

    if($key eq ""){
      #$choice = $entries[$#entries];
	  $choice = $updated_file;
      last;
    }
    elsif(defined($key)){
      $choice = $key;
      last;
    }
  }

  return $choice;
}



sub choose_default_entry {
  my($class,$dir_entries,$msg,$default) = @_;
  local $|=1;
  my $choice;

  while(1){

    my $num = 1;
    my @entries = @$dir_entries;
	my $updated_file;
	if($entries[$#entries] =~ /(.*)\.(\d+)/)
	{
		my $num = $2 + 1;
		$updated_file = $1 . ".$num"; 
	}
	print "$msg [$updated_file]:";
	
    my $key = '';
 
    if($key eq ""){      
	  $choice = $updated_file;
      last;
    }
    elsif(defined($key)){
      $choice = $key;
      last;
    }
  }

  return $choice;
}


1;
