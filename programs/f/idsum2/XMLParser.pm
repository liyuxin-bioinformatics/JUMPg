#!/usr/bin/perl -wT

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: idsum2::XMLParser

package idsum2::XMLParser;
use strict;
#use GD;
use MIME::Base64;
use Data::Dumper;#DMD June 28, 2007

my %residMass =(A => 71.0788,
		 B => 114.5962,
                 C => 103.1388,
                 D => 115.0886,
                 E => 129.1155,
                 F => 147.1176,
                 G => 57.0519,
                 H => 137.1411,
                 I => 113.1594,
                 K => 128.1741,
                 L => 113.1594,
                 M => 131.1926,
                 N => 114.1038,
                 O => 114.1472,
                 P => 97.1167,
                 Q => 128.1307,
                 R => 156.1875,
                 S => 87.0782,
                 T => 101.1051,
                 V => 99.1326,
                 W => 186.2132,
                 X => 113.1594,
                 Y => 163.1760,
                 Z => 128.6231);
sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);
  
  return $self;
}

#-----------------------------------------------------------------

sub get_runtime{
	shift @_;
	(*XML) = @_;
	seek (XML, 0, 0);
	
	my $runtime;
	while (<XML>){
		next if (!/endTime=/);
		$runtime = $_;
		$runtime =~ s/^\s+endTime=\"PT([\d\.]+)S\".*/$1/;
		last;
	}
	
	return $runtime;
} 
sub get_scannum{
	shift @_;
	(*XML, my $scan_index) = @_;
	seek (XML, $scan_index, 0);
	my $scan_num;
	while (<XML>){
		next if (!/<scan\snum=\"\d+\"/);
		$_ =~ s///g;
		chomp;
		$scan_num = $_;
		$scan_num =~ s/<scan\snum=\"(\d+)\"//;
		$scan_num = $1;
		last;
	}
	return $scan_num;
}

sub get_endpointrts{
	shift @_;
	(*XML) = @_;
	seek (XML, 0, 0);
	my ($start, $end);
	while (<XML>){
		if (/startTime/ || /endTime/){
			chomp;
			$_ =~ s///g;
			if (/startTime/){
				$start = $_;
				$start =~ s/^\s+startTime=\"PT([0-9\.]+)S\"/$1/;
				$start *= 60;
			} else {
				$end = $_;
				$end =~ s/^\s+endTime=\"PT([0-9\.]+)S\">/$1/;
				$end *= 60;
				return ($start, $end);
			}
		}
	}
}

sub get_IndexOffset{
	shift @_;
	(*XML) = @_;
	
	seek (XML, -120, 2);
	my ($indexOffset);
	while (<XML>){
	  next if (!/<indexOffset>/);
	  chomp;
	  $indexOffset = $_;
	  $indexOffset =~ s/\s+<indexOffset>(\d+)<\/indexOffset>.*/$1/o;
	  last;
	}
	return $indexOffset;
}

sub get_IndexArray{
	shift @_;
	(*XML, my $indexOffset) = @_;
	my @index_Array;
	my $lastscan = 0;

	seek (XML, $indexOffset, 0);
	while (<XML>){
	  next if (/^\s+<index/);
	  last if (/^\s+<\/index>.*/);
	  chomp;
		next if (/scan/);
	  $lastscan++;
	  my $index = $_;
	  #$index =~ s/[\s\t]+\<offset id=\"d+"\s>(\d+)<\/offset>.*/$1/o;
	  $index =~ s/[\s\t]+\<offset id="\d+"[\s]*>(\d+)<\/offset>.*/$1/o;
	  push (@index_Array, $index);
	}
	return (\@index_Array, $lastscan);
}

sub get_EntireScanInfo{
	shift @_;
  (*XML, my $scan_index) = @_;
  seek (XML, $scan_index, 0);
  while(<XML>){
		last if (/\<\/scan\>/);
		print "$_\n";
	}
}

sub get_Peaksinfo{
	shift @_;
	(*XML, my $scan_index) = @_;
	my ($peaks_count, $peaks_line);	

	seek (XML, $scan_index, 0);
	while(<XML>){
	  next if (!/peaksCount/);
          chomp;
          $peaks_count = $_;
          $peaks_count =~ s/\s+peaksCount="(\d+)".*/$1/o;
	  last;
	}
	while(<XML>){
	  if (/<peaks precision=\"32\">/ || /pairOrder=\"m\/z-int\">/){
      chomp;
      $peaks_line = $_;
	  	if (/<peaks precision="32">/){
       	$peaks_line =~ s/\s+<peaks precision="32">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} else {
       	$peaks_line =~ s/\s+pairOrder="m\/z-int">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			}
	  	last;
		} else {
			next;
		}
	}
	
	return ($peaks_count, $peaks_line);
}

sub get_MSLevel{
	shift @_;
	(*XML, my $scan_index) = @_;
	my $msLevel;
	
	seek (XML, $scan_index, 0);
	while (<XML>){
	  next if (!/msLevel=/);
	  chomp;
	  $msLevel = $_;
	  $msLevel =~ s/\s+msLevel="(\d)".*/$1/o;
	  last;
	}
	
	return $msLevel;
}
sub get_PrecursorMZINTACT{
  shift @_;
  (*XML, my $scan_index) = @_;
  my $prec_mz; my $prec_int = 0; my $prec_act = "CID";
  my $mslevel=2;                                                                                                                            
  seek (XML, $scan_index, 0);
  while (<XML>){

###### changed by yanji
    next if (!/filterLine/);

    if ((/filterLine/)) {
      chomp;
        if($_ =~ m/ms(.*)(\s\d+\.\d+)\@([a-z]+).*(\s\d+\.\d+)\@([a-z]+)/)
        {
                ($mslevel,$prec_mz,$prec_act) = ($1, $2, $3);
                $mslevel =~ s/\s+//g;
                $prec_mz =~ s/\s+//g;
                $prec_act =~ tr/a-z/A-Z/;
                last;
        }
	elsif($_ =~ m/ms(.*)(\s\d+\.\d+)\@([a-z]+)/)
      	{
     	 	($mslevel,$prec_mz,$prec_act) = ($1, $2, $3);
 	 	$mslevel =~ s/\s+//g;
      	 	$prec_mz =~ s/\s+//g;
      	 	$prec_act =~ tr/a-z/A-Z/;
      	 	last;
	}
	elsif($_ =~ / ms \[\d+\./)
	{
		$mslevel=1;
		$prec_mz = 0;
		last;
	}
    }
  }

#    next if (!/<precursorMz/);
#    if (/<\/precursorMz>/){
#      chomp;
#      $prec_mz = $_;
#      if ($prec_mz =~ /activationMethod/) {
#				if ($prec_mz =~ /precursorCharge/){
#        	$prec_mz =~ s/\s+<precursorMz precursorIntensity="([e\d\.\-\+]+)" precursorCharge="[\d]+" activationMethod="([A-Z]+)" >([e\d\.\+]+)<\/precursorMz>.*//;
#					($prec_mz, $prec_act, $prec_int) = ($3, $2, $1);
#				} else { #ADDED Feb 12 2010 for non-data dependent scans
				#print "\n$prec_mz\n";
#        	$prec_mz =~ s/\s+<precursorMz precursorIntensity="([e\d\.\-\+]+)" activationMethod="([A-Z]+)" >([e\d\.\+]+)<\/precursorMz>.*//;
#					($prec_mz, $prec_act, $prec_int) = ($3, $2, $1);
#				}
			#print "$prec_mz, $prec_act, $prec_int\n";
#      }elsif ($prec_mz =~ /precursorCharge/){
#        $prec_mz =~ s/\s+<precursorMz precursorIntensity="([e\d\.\-\+]+)" precursorCharge="[\d]+">([e\d\.\+]+)<\/precursorMz>.*//;
#        ($prec_mz, $prec_int) = ($2, $1);
#			} else {
#       $prec_mz =~ s/\s+<precursorMz precursorIntensity="([e\d\.\-\+]+)">([e\d\.\+]+)<\/precursorMz>.*//;
#       ($prec_mz, $prec_int) = ($2, $1);
#    }


#      last;
#    } 
#    else {
#      while (<XML>){
#        chomp;
#        exit if (!/collisionEnergy=/);
#        $prec_mz = $_;
#        $prec_mz =~ s/\s+collisionEnergy="[\d\.]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
#        last;
#      }
#      last;
#    }
#  }



#  if(defined $prec_mz) {
#    print "YYYYYYYYYYYanji $prec_mz\n";
#  } else {
#    print "JJJJJJJJJJJJJ None\n";	
#  }  
# {

# seek (XML, $scan_index, 0);
#    while (<XML>){
#        chomp;
	
#      if(/filterLine/)
#      {
#	$prec_mz = $_;
#        if($prec_mz=~ /Full ms2 (.*)\@cid/)
#	{
#       		 $prec_mz = $1;
#	}
#       }
#    }
#  }

  return ($prec_mz, $prec_int, $prec_act,$mslevel);
}
sub get_Precursor{
  shift @_;
  (*XML, my $scan_index) = @_;
  my $prec_mz;
                                                                                                                                                             
  seek (XML, $scan_index, 0);
  while (<XML>){
  	next if (!/<precursorMz/);
		if (/<\/precursorMz>/){
			chomp;
			$prec_mz = $_;
			#print "$prec_mz\n";
			if ($prec_mz =~ /activationMethod/){
				$prec_mz =~ s/\s+<precursorMz precursorIntensity="[e\d\.\-\+]+" precursorCharge="[\d]+" activationMethod="[A-Z]+"\s+>([e\d\.\+]+)<\/precursorMz>.*/$1/;
			
			}elsif ($prec_mz =~ /precursorCharge/){
				$prec_mz =~ s/\s+<precursorMz precursorIntensity="[e\d\.\-\+]+" precursorCharge="[\d]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
			} else {
				$prec_mz =~ s/\s+<precursorMz precursorIntensity="[e\d\.\-\+]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
			}
			last;
		} else {
			while (<XML>){
				chomp;
				exit if (!/collisionEnergy=/);
				$prec_mz = $_;
				$prec_mz =~ s/\s+collisionEnergy="[\d\.]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
				last;
			}
			last;
		}
  }
	#print "$prec_mz\n";exit;
                                                                                                                                                             
  return $prec_mz;
}

sub get_PrecursorIntensity{
  shift @_;
  (*XML, my $scan_index) = @_;
  my $prec_intensity;
                                                                                                                                                             
  seek (XML, $scan_index, 0);
  while (<XML>){
    next if (!/<precursorMz/);
		chomp;
		$prec_intensity = $_;
		$prec_intensity =~ s/[\s\t]+<precursorMz\sprecursorIntensity="([e\d\.\-\+]+)".*/$1/;
		last;
	}              
  return $prec_intensity;
}

sub get_Precursorinfo{
	shift @_;
	(*XML, my $scan_index) = @_;

	my $ce;
  seek (XML, $scan_index, 0);
  while (<XML>){
		next if (!/collisionEnergy=/);
		chomp;
		$ce = $_;
		$ce =~ s/\s+collisionEnergy="(\d+)".*/$1/;
		last;
	}
	
	my $prec_mz;
	while (<XML>){
		next if (!/<precursorMz precursorIntensity=/);
		chomp;
		$prec_mz = $_;
		print "$prec_mz\n";
		$prec_mz =~ s/\s+<precursorMz precursorIntensity="[e\d\.\-\+]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
		print "$prec_mz\n";exit;
		last;
	}
	
	return ($ce, $prec_mz);
}


sub get_BasePeakinfo{
  shift @_;
  (*XML, my $scan_index) = @_;
  my ($mz, $intensity, $rt);
                                                                                                                                                             
  seek (XML, $scan_index, 0);
  while (<XML>){
    if (/basePeakIntensity=/){
      chomp($intensity = $_);
      $intensity =~ s/.*basePeakIntensity=\"([\d\.]+)\".*/$1/o;
      last;
    }
    if (/retentionTime=/){
      chomp($rt = $_);
      $rt =~ s/.*retentionTime="PT([\d\.]+)S".*/$1/o;
    }
    next if (!/basePeakMz=/);
    chomp;
    $mz = $_;
    $mz =~ s/.*basePeakMz="([\d\.]+)".*/$1/o;
  }
                                                                                                                                                             
  return ($mz, $intensity, $rt);
}


sub get_BasePeakIntensity{
  shift @_;
  (*XML, my $scan_index) = @_;
  my $intensity;
  seek (XML, $scan_index, 0);
  while (<XML>){
    if (/basePeakIntensity=/){
      chomp($intensity = $_);
      $intensity =~ s/.*basePeakIntensity=\"([e\d\.\+]+)\".*/$1/o;
      last;
    }
  }
	if ($intensity =~ /e/){
		$intensity =~ s/([\d\.]+)e\+[0]+(\d+)/$1/;
		$intensity *= 10**$2;
	}
  return ($intensity);
}

sub get_RT{
	shift @_;
	(*XML, my $scan_index) = @_;
	my $rttime;
####### debug by xusheng ###
#	print $scan_index,"aa\n";
############################
	seek (XML, $scan_index, 0);
	while (<XML>){
	  next if (!/retentionTime=/);
	  chomp;
	  $rttime = $_;
	  $rttime =~ s/.*retentionTime="PT([\d+\.]+)S".*/$1/o;
	  last;
	}
	
	return $rttime;
}

sub get_PeaksCount{
	shift @_;
	(*XML, my $scan_index) = @_;
	my ($peaksCount);
	
	seek (XML, $scan_index, 0);
	while (<XML>){
	  next if (!/peaksCount/);
	  chomp;
	  $peaksCount = $_;
	  $peaksCount =~ s/\s+peaksCount="(\d+)".*/$1/o;
	  last;
	}
	
	return $peaksCount;
}

sub get_PeaksString{
	shift @_;
	(*XML, my $scan_index) = @_;
        my $peaks_line;
                                                                                                                                                             
        seek (XML, $scan_index, 0);
        while (<XML>){
          next if (!/<peaks precision="32">/);
          chomp;
          $peaks_line = $_;
          $peaks_line =~ s/\s+<peaks precision="32">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
          last;
        }
	return $peaks_line;
}

sub get_Peaks{
	shift @_;
	(*XML, my $peak_array, my $scan_index) = @_;
	my ($peaks, $peaks_line);
	
	seek (XML, $scan_index, 0);
	#print "$scan_index\n";
	while(<XML>){
		#print "$_\n";
		#if (/m\/z-int/){print "$_\n\n\n";	}
	  if (/<peaks\sprecision="32"/ || /pairOrder="m\/z-int"[\s]*>/){
      chomp;
			next if ($_ =~ /<peaks precision="32"\Z/);
      $peaks_line = $_;
	  	if (/<peaks precision="32">/){
       	$peaks_line =~ s/\s+<peaks precision="32"[\s\w\W\d\=\"]+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} else {
       	$peaks_line =~ s/\s+pairOrder="m\/z-int"[\s]*>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
				#print "$peaks_line\n";exit;
			}
	  	last;
		} elsif (/compressedLen/){
      chomp;
      $peaks_line = $_;
     	$peaks_line =~ s/\s+compressedLen="[\d]"\s+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			last;
		} else {
			next;
		}
	}
	#Base64 Decode
	#print "$peaks_line\n";exit;
	$peaks = decode_base64($peaks_line);
	my @hostOrder32 = unpack("N*", $peaks);
	for (@hostOrder32){
		my $float = unpack("f", pack("I", $_));
		push (@$peak_array, $float);
	}
	#print Dumper($peak_array);exit;
}

sub MS_surveytype{ #DMD Jan 28, 2008
	shift @_;
	(*XML, my $index) = @_;
	my ($smz, $emz) = get_mzrange("", *XML, $index);
	my ($peaksum) = get_PeaksCount("", *XML, $index);
	my $sourceval = $peaksum/($emz-$smz);
	return ($sourceval);
}

sub hybrid_check{
	shift @_;
	(*XML) = @_;
	my $hybrid = 0;
	my $indexOffset = get_IndexOffset("",  *XML);
  my ($index_array, $last_scan) = get_IndexArray("", *XML, $indexOffset);
	my ($smz, $emz, $smz2, $emz2);
	my ($peaksum, $peaksum2) = (0,0);
	#open (OUT, ">test_orbi.txt");
	#open (OUT2, ">test_ltq.txt");
	my $testnum = 0; my $round = 0;
	for (my $i=0; $i<$last_scan; $i++){
		my $index = $$index_array[$i];
		next if (get_MSLevel("", *XML, $index) != 1);
		next if (!defined($$index_array[$i+1]));
		my $next = $$index_array[$i+1];	next if (get_MSLevel("", *XML, $next) != 2);
		next if (!defined($$index_array[$i-1]));
		my $prev = $$index_array[$i-1]; return $hybrid if (get_MSLevel("", *XML, $prev) == 2);
		$round++;
		my ($smz, $emz) = get_mzrange("", *XML, $index); my ($peaksum) = get_PeaksCount("", *XML, $index);
		my ($smz2, $emz2) = get_mzrange("", *XML, $prev); my ($peaksum2) = get_PeaksCount("", *XML, $prev);
		my $sourceval2 = ($peaksum)/($emz-$smz);
		my $sourceval1 = ($peaksum2)/($emz2-$smz2);
	#	print OUT2 "ltq: $index $smz $emz $peaksum $sourceval2\n";
#		print OUT "orbi: $prev $smz2 $emz2 $peaksum2 $sourceval1\n";
		$testnum += 1 if ($sourceval1 < .15 && $sourceval2 > .15);
		#last if $round == 5;
	}
	#print "$testnum $round\n";
	$testnum /= $round;
	$hybrid = sprintf("%.0f", $testnum);
	#print "$testnum\n";exit;
	return $hybrid;
	#exit;
}

sub MS_Source{ #DMD June 28, 2007
	shift @_;
	(*XML) = @_;
	my $indexOffset = get_IndexOffset("",  *XML);
  my ($index_array, $last_scan) = get_IndexArray("", *XML, $indexOffset);
	my ($smz, $emz);
	my ($num, $peaksum) = (0,0);
	my ($minpeak, $maxpeak) = (10000, 0);
	for (my $i=0; $i<$last_scan; $i++){
		my $index = $$index_array[$i];
		next if (get_MSLevel("", *XML, $index) != 1);
		next if (!defined($$index_array[$i+1]));
		my $next = $$index_array[$i+1];
		next if (get_MSLevel("", *XML, $next) != 2);
		my $scan = get_scannum("", *XML, $index);
		if (!defined($smz)) {($smz, $emz) = get_mzrange("", *XML, $index);}
		my ($peakcount) = get_PeaksCount("", *XML, $index);
		$minpeak = $peakcount if ($peakcount<$minpeak);
		$maxpeak = $peakcount if ($peakcount>$maxpeak);
		$num++; $peaksum+=$peakcount;
	}
	my $sourceval = ($peaksum/$num)/($emz-$smz);
	#printf "Average # of peaks in Survey Scans = %d (scan range = %d)\n", $peaksum/$num, $emz-$smz;
	#printf "sum=$peaksum num=$num range=%d average=%d min=$minpeak max=$maxpeak final=$sourceval\n", $emz-$smz, $peaksum/$num;
	return ($sourceval);
}

sub get_mzrange{ #DMD June 28, 2007
	shift @_;
	(*XML, my $scan_index) = @_;
  my ($smz, $emz);

  seek (XML, $scan_index, 0);
  while (<XML>){
    next if (!/startMz=/ && !/endMz=/);
    chomp;
		if (/startMz=/){
    	$smz = $_; $smz =~ s/.*startMz=\"([\d+\.]+)\".*/$1/o;
		} else {
    	$emz = $_; $emz =~ s/.*endMz=\"([\d+\.]+)\".*/$1/o;
    	last;
		}
  }

  return ($smz, $emz);
}
