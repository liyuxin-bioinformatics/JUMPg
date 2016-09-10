#!/usr/bin/perl

use strict;
use warnings;
use Clone 'clone';

if (scalar(@ARGV)!=1)
{
        die "usage: perl eventSummary.pl alg.updated2.pep\n";
}

#----------------------------------------------------
my $nearPepRange=2000;
#----------------------------------------------------

my (%events,%eventType);
open(IN,"$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>) {
	chomp;
        my ($outfile,$XCorr,$dCn,$peptide,$scanCount,$readCount,$readID,$ucscID,$mismatch,$gap,$readLength,$alignLength,$qStart,$qEnd,$levels,$sStart,$sEnd,$sStrand,$b5,$b3,$strand,$frame,$sType)=split /\t/,$_;

	next if ($sType eq 'CDS' and $frame == 0); # CDS in frame

	my $eventKey = "$sType\|$ucscID\|$strand\|$frame";
	$events{$sType}{$eventKey}{$peptide}=$_;
}
close IN;

# assign events
foreach my $sType (keys %events) {
	if ($sType eq 'intronic' || $sType eq 'ncGene') {
		foreach my $eventKey (keys %{$events{$sType}}) {
			$eventType{$sType}{$eventKey}=clone($events{$sType}{$eventKey});
			#$eventType{$sType}{$eventKey}=$events{$sType}{$eventKey};
		}
	} elsif ($sType eq 'mRNA') {
		foreach my $eventKey (keys %{$events{$sType}}) {
			$eventType{'UTR'}{$eventKey}=clone($events{$sType}{$eventKey});
			#$eventType{'UTR'}{$eventKey}=$events{$sType}{$eventKey};
		}
	} elsif ($sType eq 'CDS') {
		foreach my $eventKey (keys %{$events{$sType}}) {
			my ($sType,$ucscID,$strand,$frame)=split /\|/,$eventKey;
			if ($strand eq '+' && $frame != 0) {
				$eventType{'frameShift'}{$eventKey}=clone($events{$sType}{$eventKey});
			} elsif ($strand eq '-') {
				$eventType{'antiSenCDS'}{$eventKey}=clone($events{$sType}{$eventKey});
			}
		}
	} elsif ($sType eq 'intergenic') {
		foreach my $eventKey (keys %{$events{$sType}}) {
			my @pepPos;
			foreach my $peptide  (keys %{$events{$sType}{$eventKey}}) {
	        		my ($outfile,$XCorr,$dCn,$peptide,$scanCount,$readCount,$readID,$ucscID,$mismatch,$gap,$readLength,$alignLength,$qStart,$qEnd,$levels,$sStart,$sEnd,$sStrand,$b5,$b3,$strand,$frame,$sType)=split /\t/,$events{$sType}{$eventKey}{$peptide};
				push @pepPos, $b5; # start genomic pos
			}
			@pepPos=sort {$a <=> $b} @pepPos; # sort pep pos 
			my $k=0; my @range; # record range of each segment 
			$range[$k][0]=my $crtPos=$pepPos[0]; # start of the 1st segment
			for (my $i=1;$i<=$#pepPos;$i++) {
				if ( $pepPos[$i] < $crtPos + $nearPepRange ) { # same segment
					$crtPos = $pepPos[$i] ; 
				} else { # new segment
					$range[$k][1] = $crtPos + $nearPepRange; # end of this segment
					$k++; # begins a new segment
					$range[$k][0]=$crtPos=$pepPos[$i]; # start of new segment
				}
			}
			$range[$k][1] = $crtPos + $nearPepRange; # end of last segment
			foreach my $peptide (keys %{$events{$sType}{$eventKey}}) {
				for (my $i=0; $i<=$k; $i++) {
	        			my ($outfile,$XCorr,$dCn,$peptide,$scanCount,$readCount,$readID,$ucscID,$mismatch,$gap,$readLength,$alignLength,$qStart,$qEnd,$levels,$sStart,$sEnd,$sStrand,$b5,$b3,$strand,$frame,$sType)=split /\t/,$events{$sType}{$eventKey}{$peptide};
					if ($range[$i][0] <= $b5 and $b5 <= $range[$i][1]) {
						$eventType{'intergenic'}{"$eventKey|$range[$i][0]\-$range[$i][1]"}{$peptide}=$events{$sType}{$eventKey}{$peptide};
						last;
					}
				}

			}


		}
	}
}

# print
print "eventType\tevenyKey\tpeptideN\tpeptides\n";
foreach my $sType (keys %eventType) {
	foreach my $eventKey (keys %{$eventType{$sType}}) {
		print "$sType\t$eventKey\t";
		my $k=0;
		foreach my $peptide (keys %{$eventType{$sType}{$eventKey}}) {
			print "$peptide,";
			$k++;
		}
		print "\t$k\n";
	}
}
