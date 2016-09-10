#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::MassCorrection

##################################################################
##	Mass table calculated by Junmin Peng on 10/03/2014	##
##################################################################
## 	TMT reporters and y1-ions
##	 	126 = 126.1277259380
## 		127N = 127.1247608314
##		127C = 127.1310807758
##		128N = 128.1281156692
##		128C = 128.1344356136
##		129N = 129.1314705070
##		129C = 129.1377904514
##		130N = 130.1348253448
##		130C = 130.1411452892
##		131 = 131.1381801826
##		K = 376.2757362992
##		R = 175.1189521741
##	y1-ions of non-TMT
##		K = 147.1128041645
##		R = 175.1189521741
##	ESI MS1-ion
##		(Si(CH3)2O))6 + H+ = 445.1200245337

package Spiders::MassCorrection;
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new {
	my ($class,%arg)=@_;
    my $self = {};
    bless $self, $class;
	return $self;
}

sub massCorrection {
	##########################
	## Initialization	##
	##########################
	
	my ($self, $msHash, $ms2Hash, $mzArray, $params) = @_;
	my @options = split(/,/, $$params{'mass_correction'});
	my $option = shift (@options);
	my $manualMassShift = shift (@options);
	if (defined $manualMassShift) {
		$manualMassShift =~ s/\s//;
	}
	my $tolPpm = 50;				## Window-width for each reference ion
	my $thresholdPercentage = 0.05;			## Threshold ratio of the number of spectra used for the mass-shift calculation to the total number of spectra	 
							## e.g. $thresholdPercentage = 0.2 indicates that at least 20% of spectra (MS1 or MS2)
							## should be used to calculate mass shifts
									
	##################################
	## Mass-shift correction module	##
	##################################
	
	if ($option == 0) {
		##########################
		## No correction	##
		##########################
		print "  No mass-shift correction\n";
		return ($ms2Hash, $mzArray);
	} elsif ($option == 1) {
		##################################################
		## Mass-shift correction based on MS1 spectra	##
		##################################################
		my @referenceMass = (445.1200245337);
		my $referenceIndex = int($referenceMass[0]);
		my @massShifts;
		my $nScans = 0;
		for (my $scanNumber = 0; $scanNumber < scalar(@{$mzArray}); $scanNumber++) {
			if (defined $$mzArray[$scanNumber]) {
				##################################
				## Mass-shift calculation	##
				##################################
				$nScans++;
				my @masses = keys %{$$mzArray[$scanNumber][$referenceIndex]};
				my @intensities = values %{$$mzArray[$scanNumber][$referenceIndex]};
				my @sortedIndex = sort {$masses[$a] <=> $masses[$b]} 0..$#masses;
				@masses = @masses[@sortedIndex];
				@intensities = @intensities[@sortedIndex];
				my $observedMass = (getObservedMasses(\@masses, \@intensities, \@referenceMass, $tolPpm))[0];
				## Mass shift (ppm-scale for each reference mass)
				## "ppm" is better than "Da" since it is a relative scale
				## e.g. 1 Da is a large correction for small masses, but
				## 		it may be negligible factor for large masses (>1000)
				if ($observedMass ne "") {
					my $ppm = ($observedMass - $referenceMass[0]) / $referenceMass[0] * 1e6;
					push (@massShifts, $ppm);
				}
			}
		}
		
		##########################################################################################
		## Filter the calculated mass shifts							##
		## 1. More than 50% of spectra should be used to calculate mass shifts			##
		## 2. Remove 20% highest/lowest mass shifts to get a more reliable correction factor	## 
		##########################################################################################
		my $nMassShifts = scalar(@massShifts);
		my $measuredPercentage = $nMassShifts / $nScans;	
		## More than 5% of spectra should be used to calculate mass shifts
		## Otherwise, mass-shift correction will not be performed
		if ($measuredPercentage >= $thresholdPercentage) {
			printf("  %.2f%% of spectra (%d of %d) is used to calculate mass shifts\n", $measuredPercentage * 100, $nMassShifts, $nScans);
			@massShifts = sort {$a <=> $b} @massShifts;
			## Filter the 20% of highest/lowest mass shifts
			my $numFiltered = int(0.2 * scalar(@massShifts));
			splice (@massShifts, 0, $numFiltered);					## Remove the lowest mass shifts
			splice (@massShifts, scalar(@massShifts) - $numFiltered, $numFiltered);	## Remove the highest mass shifts
		} else {
			printf("  Only %.2f%% of spectra (%d of %d) is used to calculate mass shifts\n", 
					$measuredPercentage * 100, scalar(@massShifts), $nScans);
			printf("  Mass correction will not be performed when less than %.2f%% of spectra is used\n", $thresholdPercentage * 100);
			print "  Please choose another option like,\n";
			print "    \"0\" (no correction) or\n";
			print "    \"3\" (manual correction with a specified mass-shift value)\n";
			print "  for the \"mass_correction\" parameter and run again\n";
			exit;
		}
		
		##########################################################################
		## Calculate a global "correction factor" by taking the mean value of	##
		## mass shifts calculated from MS1 scans (i.e. mean of @massShifts)	##
		##########################################################################
		my $meanMassShift = mean(@massShifts);
		my $stdMassShift = stdev(@massShifts);
		printf ("  Calculated mass-shift: mean = %.5f ppm and SD = %.5f ppm\n", $meanMassShift, $stdMassShift); 
		my $massCorrectionFactor = mean(@massShifts);
		
		##########################################################################
		## Mass-shift correction of all MS1 spectra using the correction factor	##
		##########################################################################
		for (my $scanNumber = 0; $scanNumber < $nScans; $scanNumber++) {
			## Correction of MS1 masses
			if (defined $$mzArray[$scanNumber]) {
				for (my $i = 0; $i < scalar(@{$$mzArray[$scanNumber]}); $i++) {
					next if (!defined $$mzArray[$scanNumber][$i]);
					## Keys of %{$$mzArray[$scanNumber][$i]} (correspond to mz values) need to be changed
					%{$$mzArray[$scanNumber][$i]} = map {
						my $newKey = $_ / (1 + $massCorrectionFactor / 1e6);
						$newKey => $$mzArray[$scanNumber][$i]{$_}
					} keys (%{$$mzArray[$scanNumber][$i]});
				}
			}
		}
		print "  Mass-shift correction has been finished for MS1 spectra\n";
		return ($ms2Hash, $mzArray);
	} elsif ($option == 2) {
		##################################################
		## Mass-shift correction based on MS2 spectra	##
		##################################################
		
		## Modified by JCho on 12/29/2014
		#my @referenceMasses = (126.1277259380, 147.1128041645, 175.1189521741, 376.2757362992);
		my @referenceMasses = (126.1277259380, 147.1128041645, 175.1189521741);
		if ($$params{'add_K_Lysine'} > 0) {
			my $modK = 147.1128041645 + $$params{'add_K_Lysine'};
			push (@referenceMasses, $modK);
		}
		if ($$params{'add_R_Arginine'} > 0) {
			my $modR = 175.1189521741 + $$params{'add_R_Arginine'};
			push (@referenceMasses, $modR);
		}
		## End of modification on 12/29/2014

		my $nReferences = scalar(@referenceMasses);
		my @massShifts;
		my $nScans = 0;
		## Calculate mass shifts between observed reference masses and theoretical reference masses
		foreach my $scanNumber (keys %{$ms2Hash}) {
			next if (!defined $$ms2Hash{$scanNumber}{'prec_mz'});
			$nScans++;
			my @observedMasses = getObservedMasses(\@{$$ms2Hash{$scanNumber}{'msms_mz'}}, \@{$$ms2Hash{$scanNumber}{'msms_int'}}, \@referenceMasses, $tolPpm);
			for (my $i = 0; $i < $nReferences; $i++) {
				if ($observedMasses[$i] ne "") {
					my $ppm = ($observedMasses[$i] - $referenceMasses[$i]) / $referenceMasses[$i] * 1e6;
					push (@{$massShifts[$i]}, $ppm);
				}
			}
		}
		
		##########################################################################
		## Choose a reference (ion) with the largest number of data points	##
		##########################################################################
		my $ind = 0;
		my $nMassShifts = 0;
		for (my $i = 0; $i < $nReferences; $i++) {
			if (defined $massShifts[$i] && scalar(@{$massShifts[$i]}) > $nMassShifts) {
				$ind = $i;
				$nMassShifts = scalar(@{$massShifts[$i]});
			}
		}
		if ($nMassShifts == 0) {
			print "  Any of reference ions (TMTreporter-126, y1-ions of K and R) is not found\n\n";
			print "  Please choose another option like,\n";
                        print "    \"0\" (no correction) or\n";
                        print "    \"3\" (manual correction with a specified mass-shift value)\n";
                        print "  for the \"mass_correction\" parameter and run again\n";
			exit; 
		} else {
			@massShifts = @{$massShifts[$ind]};
		}
		
		##########################################################################################
		## Filter the calculated mass shifts							##
		## 1. More than 50% of spectra should be used to calculate mass shifts			##
		## 2. Remove 20% highest/lowest mass shifts to get a more reliable correction factor	## 
		##########################################################################################
		$nMassShifts = scalar(@massShifts);
		my $measuredPercentage = $nMassShifts / $nScans;
		## More than 50% of spectra should be used to calculate mass shifts
		## Otherwise, mass-shift correction will not be performed
		if ($measuredPercentage >= $thresholdPercentage) {
			printf("  %.2f%% of spectra (%d of %d effective MS2) is used to calculate mass shifts\n",
					$measuredPercentage * 100, $nMassShifts, $nScans);
			@massShifts = sort {$a <=> $b} @massShifts;
			## Filter the 20% of highest/lowest mass shifts
			my $numFiltered = int(0.2 * scalar(@massShifts));
			splice (@massShifts, 0, $numFiltered);									## Remove the lowest mass shifts
			splice (@massShifts, scalar(@massShifts) - $numFiltered, $numFiltered);	## Remove the highest mass shifts
		} else {
			printf("  Only %.2f%% of spectra (%d of %d effective MS2) is used to calculate mass shifts\n", 
					$measuredPercentage * 100, scalar(@massShifts), $nScans);
			printf("  Mass correction will not be performed when less than %.2f%% of spectra is used\n", $thresholdPercentage * 100);
			print "  Please choose another option like,\n";
			print "    \"0\" (no correction) or\n";
			print "    \"3\" (manual correction with a specified mass-shift value)\n";
			print "  for the \"mass_correction\" parameter and run again\n";
			exit;
		}
		
		##################################################################################
		## Calculate a global "correction factor" by taking the mean value of		##
		## mass shifts calculated from MS2 or MS3 scans (i.e. mean of @massShifts)	##
		##################################################################################
		my $meanMassShift = mean(@massShifts);
		my $stdMassShift = stdev(@massShifts);
		printf ("  Calculated mass-shift: mean = %.5f ppm and SD = %.5f ppm\n", $meanMassShift, $stdMassShift); 
		my $massCorrectionFactor = mean(@massShifts);
		
		##################################################################################
		## Mass-shift correction of all MS1/MS2 spectra using the correction factor	##
		##################################################################################
		foreach my $scanNumber (keys %{$ms2Hash}) {
			## Correction of MS2 masses
			if (defined $$ms2Hash{$scanNumber}{'prec_mz'}) {
				for (my $i = 0; $i < scalar(@{$$ms2Hash{$scanNumber}{'msms_mz'}}); $i++) {
					$$ms2Hash{$scanNumber}{'msms_mz'}[$i] = $$ms2Hash{$scanNumber}{'msms_mz'}[$i] / (1 + $massCorrectionFactor / 1e6);
				}
			}
			## Correction of MS1 masses
			if (defined $$mzArray[$scanNumber]) {
				for (my $i = 0; $i < scalar(@{$$mzArray[$scanNumber]}); $i++) {
					next if (!defined $$mzArray[$scanNumber][$i]);
					## Keys of %{$$mzArray[$scanNumber][$i]} (correspond to mz values) need to be changed
					%{$$mzArray[$scanNumber][$i]} = map {
						my $newKey = $_ / (1 + $massCorrectionFactor / 1e6);
						$newKey => $$mzArray[$scanNumber][$i]{$_}
					} keys (%{$$mzArray[$scanNumber][$i]});
				}
			}
		}
		print "  Mass-shift correction has been finished for MS1 and MS2 spectra\n";
		return ($ms2Hash, $mzArray);
	} elsif ($option == 3) {
		##################################################
		## Manual mass-correction using the user input	##
		##################################################
		
		print "  Mass-shift (user specified) is $manualMassShift\n";
		for (my $scanNumber = 0; $scanNumber < scalar(@{$mzArray}); $scanNumber++) {
			## Correction of MS1 masses
			if (defined $$mzArray[$scanNumber]) {
				for (my $i = 0; $i < scalar(@{$$mzArray[$scanNumber]}); $i++) {
					next if (!defined $$mzArray[$scanNumber][$i]);
					## Keys of %{$$mzArray[$scanNumber][$i]} (correspond to mz values) need to be changed
					%{$$mzArray[$scanNumber][$i]} = map {
						my $newKey = $_ / (1 + $manualMassShift / 1e6);
						$newKey => $$mzArray[$scanNumber][$i]{$_}
					} keys (%{$$mzArray[$scanNumber][$i]});
				}
			}
		}
		print "  Mass-shift correction has been finished for MS1 spectra\n";
		return ($ms2Hash, $mzArray);				
	} else {
		print "  Stopped\n\n";
		print "  Set a correct parameter for the mass-shift correction in your .params file and run again\n";
		print "  Please see the line of \"mass_correction = \"\n";
		print "  mass_correction = 1, for AUTOMATIC MS1-based mass-shift correction (It will correct MS1 spectra)\n";
		print "  mass_correction = 2, for AUTOMATIC MS2-based mass-shift correction (It will correect MS1 and MS2 spectra)\n";
		print "  mass_correction = 3, for MANUAL MS1-based mass-shift correction (It will correct MS1 spectra)\n";
		print "    If you choose a MANUAL mass-shift correction, \n";
		print "    you have to specify the mass-shift value in a ppm-scale followed by , (comma) as following\n";
		print "    e.g. mass_correction = 3, 5   => 5 ppm will be subtracted from all MS1 spectra\n";
		print "                                     The corrected m/z values will be 5-ppm smaller than the original values\n";
		print "         mass_correction = 3, -10 => -10 ppm will be subtracted from all MS1 spectra\n";
		print "                                     The corrected m/z values will be 10-ppm larger than the original values\n";
		exit;
	}
}

##################
## Subroutines	##
##################

sub stdev {
	my (@x) = @_;
	my $n = scalar(@x);
	my $stdev;	
	if ($n > 1) {
		my $mean = mean(@x);
		my $sum = 0;
		foreach (@x) {
			$sum = $sum + ($_ - $mean) ** 2;
		}
		$stdev = sqrt($sum / ($n - 1));
	} else {
		$stdev = 0;
	}
	return ($stdev);
}

sub mean {
	my (@x) = @_;
	my $n = scalar(@x);
	my $sum = 0;
	foreach (@x) {
		$sum = $sum + $_;
	}
	my $mean = $sum / $n;
	return ($mean);
}

sub getObservedMasses {
	my ($massRef, $intensityRef, $referenceMasses, $tolPpm) = @_;
	my @masses = @{$massRef};
	my @intensities = @{$intensityRef};
	my @obsMasses = ();
	my $nPeaks = scalar(@masses);
	my $nReferences = scalar(@{$referenceMasses});
	for (my $i = 0; $i < $nReferences; $i++) {
		my $obsMass = "";
		my $obsIntensity = 0;	## This parameter can control the selection of a low-intensity peak
		my $uL = $$referenceMasses[$i] + $$referenceMasses[$i] * $tolPpm / 1e6;
		my $lL = $$referenceMasses[$i] - $$referenceMasses[$i] * $tolPpm / 1e6;
		for (my $j = 0; $j < $nPeaks; $j++) {
			if ($masses[$j] >= $lL && $masses[$j] <= $uL) {
				if ($intensities[$j] > $obsIntensity) {
					$obsMass = $masses[$j];
					$obsIntensity = $intensities[$j];
				}
			}
			last if ($masses[$j] > $uL);
		}
		push (@obsMasses, $obsMass);
	}
	return (@obsMasses);
}

1;
