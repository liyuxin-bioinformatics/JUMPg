#!/usr/bin/perl -w

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: idsum2::CL

package idsum2::CL;
use strict;
#use Term::ReadKey;
  
###########################################################################
#This is pretty much straight out of the perl cookbook. pg. 123, Tom Christiansen and
#nathaniel torkington
sub printcolumns {

  my($class,$list) = @_;
  my $maxlength =1;
  my ($item,$cols,$rows,$maxlen);
  my($xpixel,$ypixel,$mask);    
  
  ($cols,$rows,$xpixel,$ypixel) = GetTerminalSize();
  
  $maxlen = 1;
  foreach (@$list) {
    my $mylen;
    $maxlen = $mylen if(($mylen = length) > $maxlen);
  }
  
  $maxlen++;
  
  $cols = int($cols/$maxlen) || 1;
  $cols = 1 if($#$list < 10);
  $rows = int(($#$list+$cols) / $cols);
  
  #pre-create mask for faster computation
  $mask = sprintf("%%-%ds ", $maxlen-1);
    
  for($item = 0;$item < $rows*$cols; $item++){
    my $target = ($item % $cols) * $rows + int($item/$cols);
    my $piece = sprintf($mask,$target < @$list ? $list->[$target] : "");
    $piece =~ s/\s+$// if(($item+1) % $cols == 0); #don't blank pad to EOL
    print $piece;
    print "\n" if (($item+1) % $cols ==0);
  }
  
  print "\n" if(($item+1) % $cols ==0);
}

###########################################################################

sub choose_dir_entry {
  my($class,$dir_entries,$msg,$default) = @_;
  local $|=1;
  my $choice;
  
  while(1){
  
    my $num = 1;
    my @entries = map{"[".$num++ . "] $_      "} @$dir_entries;
    printcolumns($class,\@entries);
        
    if(defined($msg)) {print $msg;}
    else {print "Choose one";} 
    if($default){print " [$default]: ";}
    else {print ": ";}
    
    my $key = <STDIN>;
    chomp($key);
        
    if($key eq "" && $default){
      $choice = $default;
      last;
    }
    if(--$key <= $#$dir_entries){
      $choice = $dir_entries->[$key];
      last;
    }
  }
  
  return $choice;
}

#---------------------------------------------------------------------------
#simple tool to ask a yes or no question. Requires Term::ReadKey
=head
sub ask_bool {
  my ($class,$message) = @_;
  local $|=1;
  ReadMode 'cbreak';
  my ($result,$key);
  $result = 0;
  while(1){
    print "$message [y/n]: ";
    $key = ReadKey(0);
    print $key;
    ReadMode 'normal';
    print "\n" if!($key eq "\n");
    if($key=~/([yYnN])/o){
      $result = ($1 =~ /[yY]/);
      last;
    }
  }
  return ($result);
}
=cut
#---------------------------------------------------------------------------
sub ask {
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

sub ask1 {
  my ($class,$message,$default) = @_;
  local $| = 1;
  my ($result,$choice);
    return $default;
}

#---------------------------------------------------------------------------
#like printcolumns, prints a list out across the screen, but sorts it so that 
#lines go across rather than down.

sub printrows {
  my ($class,$array) = @_;
  my @lengths = ();
  my ($a,$l,$i,$test);
  foreach $a (@$array){
    for($i=0;$i<=$#$a;$i++){
      $test = $lengths[$i] || 0;
      $l = length($a->[$i]);
      $lengths[$i] = ($test > $l) ? $test : $l;
    }
  }
  
  
  foreach $a (@$array){
    for($i=0;$i<=$#$a;$i++){
      my $maxlen = $lengths[$i];
      my $mask = sprintf("%%-%ds ",$maxlen-1);
      printf($mask,$a->[$i]);
    }
    print "\n";
  }
  
  print "\n";
}


1;
