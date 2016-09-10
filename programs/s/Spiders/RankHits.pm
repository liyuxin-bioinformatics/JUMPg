#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::RankHits

package Spiders::RankHits;

######### Deisotope ##########################################
#                                                             #
#       **************************************************    #  
#       **** Postprocessing for search		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2012 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################


require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw();
use vars qw($VERSION @ISA @EXPORT);
use strict; 

$VERSION     = 1.01;

use List::Util qw(first);
use File::Basename;

############## Example ##################
## Get database info
## use Spiders::RankHits;
## use Spiders::Params;
##
## my $Rankhit = new Spiders::RankHits();
## my $p = Spiders::Params->new('-path'=>$parameter);
## my $params=$p->parse_param();
## $Rankhit->set_parameter($params);
## my ($miscleavage_modification_coeff,$modification_coeff,$miscleavage_coeff) = $Rankhit->get_db_mis_mod();
## my $mainhash	= parse_spOut_files_v5($folder);
## my $panelty_mainhash = $Rankhit->get_miscleavage_mod_spout($mainhash);
## my $best_weight = $Rankhit->calculate_best_weight($panelty_mainhash);
## $Rankhit->rerank_Jscore($panelty_mainhash,$miscleavage_modification_coeff,$modification_coeff,$miscleavage_coeff);

sub new {
  my ($class,%args) = @_;
  my $self = {};
  bless $self;
 
  return $self;
}


sub get_mis_cleavage_number
{
	my ($self,$peptide) = @_;
	
	my $mis = ($peptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;
	return $mis;
} 

sub get_mod_site_number
{
	my ($self,$peptide) = @_;
	my $mods = ($peptide =~ s/([\@\#\*\^\~\$]+)/$1/g) || 0;
	return $mods;
}

sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}



sub parse_spOut_files_v5{
	my ($self,$folder)= @_;
	my $mainhash;
	my $maxHitCondiered = 10;
	
	if(!-e $folder){ die "Error opening folder $folder\n";}
	
	my @files = glob("$folder/*.spout");
	
	foreach my $spout (@files){

		print "\r  Reading file $spout                   ";		
		open(INPUT,$spout) or die "Could onpen $spout. Error:!$\n";
		my @content = <INPUT>;		
		$spout = basename($spout);
		$spout =~ s/\.spout//;		
		
		if(scalar(@content)>0)
		{ #if the file is not empty
	#		($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"},$mainhash->{$spout}->{"outhash"}->{"Tag"},$mainhash->{$spout}->{"outhash"}->{"Tag_E_value"},$mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}) = grep(/Precursor mass\s*=/ || /Tag\s*=/ || /Tag E_value\s*=/ || /Tag Side Mass\s+=/,@content);

			($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}) = grep(/Precursor mass\s*=/,@content);
			($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}) = grep(/Tag number\s*=/,@content);
								
			chomp($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}); 
			
             if ($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} =~ /^Precursor mass\s*=\s*([0-9\.]+)[\s\t]+Percentage of precursor peak intensity\s*=\s*([\d\.\%]+)/)
            {
                        $mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} =$1;
                        $mainhash->{$spout}->{"outhash"}->{"PPIpct"} =$2;
            }
            if ($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"} =~ /^MS2 signal noise ratio\s*=\s*([0-9\.]+)[\,\s\t]+Tag number\s*=\s*([\d]+)/)
            {
                    $mainhash->{$spout}->{"outhash"}->{"MS2SN"}=$1;
                    $mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}=$2;
            }				
			
		#	$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];

		#	$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/\s+/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];
		#	chomp($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}); $mainhash->{$spout}->{"outhash"}->{"TotalTagNum"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}))[1];
	#		chomp($mainhash->{$spout}->{"outhash"}->{"Tag_E_value"}); $mainhash->{$spout}->{"outhash"}->{"Tag_E_value"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"Tag_E_value"}))[1];
	#		chomp($mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}); $mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}))[1];
			
			my $startIndex = first {$content[$_]=~/-+\n/} 0..$#content;

			foreach my $i ($startIndex+1..$#content)
			{	
				
				last if ($content[$i]=~/^\s*$/);
				last if ($content[$i]=~/All identified peptide/);
				chomp($content[$i]); #Remove return line
				$content[$i]=~ s/^\s+//; #Remove leading spaces
				my @vals = split /\s+/,$content[$i];
				if(scalar(@vals)==1)
				{
					next;
				}
				if(scalar(@vals)==10)
				{
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"TagRank"} =0;
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[8];	
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Weight_E_value"} = $vals[6];					
				}
			#	print $spout,"\t",scalar(@vals),"\n";
				elsif(scalar(@vals)==11){
				#	die "Error parsing file $spout at l:qine: $i\n Some values are missing\n";
				
				
					if($vals[0]!~/\d+/){
						next;
					#	die "Error parsing file $spout at line $i\n The Oder field should be a numeric value\n";
					}
				#	last if ($vals[0]>$maxHitCondiered);
				#	next if ($vals[9]=~/\@/);
					
					$self->checkParameters($i,@vals);
					#print $spout,$vals[8],"\n"  if ($vals[8]<0.1);
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"MH"}= $vals[1];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"lsideMass"}= $vals[2];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"rsideMass"}= $vals[3];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide_E_value"}= $vals[4];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag"} = $vals[5];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"TagRank"} = $vals[6];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag_num"} = $vals[7];				
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Weight_E_value"} = $vals[8];

								
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[9];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide"}= $vals[10];

					my $intpep = $vals[10];
					$intpep =~ s/[A-Z\-]\.([A-Z\@\#\*]+)\.[A-Z\-]/$1/;
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"intpep"}= $intpep;
									
			#		$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"modification"}= $vals[7];
				}
			}
		}	 
		close($spout);		
	}

	return $mainhash;
}

sub get_score_scans
{
	my ($self,$mainhash,$cutoff_score,$dta_path) = @_;
	my @research_scans;
	foreach my $spout (keys %$mainhash)
	{
		my $tagseq = $mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"};
		$tagseq =~ s/[^A-Za-z]//g;
		next if(length($tagseq) >3);
		next if($mainhash->{$spout}->{"outhash"}->{"MS2SN"}<5);
		next if($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"}>$cutoff_score);
		$spout = $spout . ".dta";
		push (@research_scans,"$dta_path/$spout");
	}		
	return \@research_scans;
}

sub calculate_random_matching_score
{
	my ($self,$mainhash,$cutoff_score)= @_;
	my ($sum, $average, $n) = (0,0,0);
	foreach my $spout (keys %$mainhash)
	{
		next if($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"}>$cutoff_score);
		
		my $protein = $mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"protein"};
		if($protein =~ /\#\#Decoy/)
		{
			$n++;
			$sum += $mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"};
		}
	}
	$average = $sum / $n;
	return $average;
}

sub rerank_Jscore
{
	my ($self,$newline,$weight)= @_;
	my $mainhash;
	my $maxHitCondiered = 10;
	my $newlines;
	my $all_peptide_info;
	my $peptide_Jscore;
	my $phash = $self->get_db_mis_mod();	
	print "\n  Reranking with weight: $weight\n";
	foreach my $spout (keys %$newline)
	{
		print "\r  Parsing file $spout  ";
		my @content = @{$newline->{$spout}};
        my @content_copy = @content;		
		my ($startIndex1, $startIndex2) = (0,0);
        if (scalar(@content) > 0) 
		{
#			print "Before:",@content,"\n";
            # my $startIndex = first {$content[$_] =~ /-+\n/} 0 .. $#content;
            ($startIndex1, $startIndex2) =
                (grep {$content[$_] =~ /^-{2,}\n/} 0 .. $#content)[0, 1];

			($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}) = grep(/Precursor mass\s*=/,@content);
			($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}) = grep(/Tag number\s*=/,@content);
			
			chomp($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}); $mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];

			$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/\s+/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];
			chomp($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}); $mainhash->{$spout}->{"outhash"}->{"TotalTagNum"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}))[1];
			
			my $startIndex = first {$content[$_]=~/-+\n/} 0..$#content;

			foreach my $i ($startIndex+1..$#content)
			{	
				
############## not re-rank those PSMs whose score is larger than 				
								
				last if ($content[$i]=~/^\s*$/);
				last if ($content[$i]=~/All identified peptide/);
				chomp($content[$i]); #Remove return line
				$content[$i]=~ s/^\s+//; #Remove leading spaces
				my @vals = split /\s+/,$content[$i];
				if(scalar(@vals)==1)
				{
					next;
				}

				elsif(scalar(@vals)==11){				
					if($vals[0]!~/\d+/){
						die "Error parsing file $spout at line $i\n The Oder field should be a numeric value\n";
					}
					
					$self->checkParameters($i,@vals);
#					print $spout,$vals[8],"\n"  if ($vals[8]<0.1);
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"MH"}= $vals[1];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"lsideMass"}= $vals[2];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"rsideMass"}= $vals[3];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide_E_value"}= $vals[4];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag"} = $vals[5];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"TagRank"} = $vals[6];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag_num"} = $vals[7];				
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Weight_E_value"} = $vals[8];

								
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[9];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide"}= $vals[10];

					my $intpep = $vals[10];
					$intpep =~ s/[A-Z\-]\.([A-Z\@\#\*]+)\.[A-Z\-]/$1/;
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"intpep"}= $intpep;

				}
			}
		
            # parsing part 3 of input file
            # ============================

            my $peptide;
            for my $i (($startIndex2 + 1) .. $#content) 
			{
                last if ($content[$i] =~ /^\s*$/);
                chomp($content[$i]); # Remove return line
                $content[$i] =~ s/^\s+//; # Remove leading spaces
                my @vals = split /\s+/, $content[$i];
                if ($#vals == 11) {
                    my @data = (my ($MH, $lsideMass, $rsideMass,
                        $peptide_E_value, $Tag, $TagEvalue, $TagRank,
                        $reference), $peptide) =
                        @vals[1, 2, 3, 4, 5, 6, 7, 10, 11];
						

					my $intpep = $peptide;
                    $intpep =~ s/[A-Z\-]\.([A-Z\@\#\*]+)\.[A-Z\-]/$1/;
						
                    my $mis_cleavage_num = $self->get_mis_cleavage_number($intpep);
					
                    my $mod_site_num = $self->get_mod_site_number($intpep);
########################################because multiple same #########################################
					$all_peptide_info->{$spout}->{$peptide}->{'Num_of_Tags'}++;

#					my $Jscore = $peptide_E_value + $TagEvalue - $weight * ($phash->{"mis$mis_cleavage_num"}); 
					my $Jscore = $peptide_E_value + $TagEvalue - $weight * ($phash->{"mis$mis_cleavage_num"} + $phash->{"mod$mod_site_num"}) * ($peptide_E_value + $TagEvalue);
#					print $phash->{"mis$mis_cleavage_num"},"\t",$phash->{"mod$mod_site_num"},"\t",$Jscore,"\n";
					if($Jscore <= 0)
					{
						$Jscore = 0.5;
					}
					if(!defined($all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'}))
					{
						$all_peptide_info->{$spout}->{$peptide}->{'data'} = \@data;
						$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'} = $Jscore;
						$peptide_Jscore->{$spout}->{$peptide} = $Jscore;						
					}
					elsif($Jscore>$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'})
					{

						$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'} = $Jscore;
						$peptide_Jscore->{$spout}->{$peptide} = $Jscore;
					}

		####### do not put the duplicates references ##################			
					unless($all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference})
					{					
						push @{$all_peptide_info->{$spout}->{$peptide}->{'ref'}}, $reference;
						$all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference} = 1;
					}
					
                } elsif (@vals == 1) {
                    my $reference = $vals[0];
		####### do not put the duplicates references ##################			
					unless($all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference})
					{
						push @{$all_peptide_info->{$spout}->{$peptide}->{'ref'}}, $reference;
						$all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference} = 1;						
					}
                } 
				else {
                    printf STDERR "Warning: input line %d not parsed.\n", $i + 1;
                }
            }
		}

		my @top_three = (
                map {$_->[0]}
                sort {$b->[1] <=> $a->[1]}
                map {[$_, $peptide_Jscore->{$spout}->{$_}]}
                keys %{$peptide_Jscore->{$spout}})[0 .. 2];
				
		my @newlines_part2;
		my $order = 1;
		for my $peptide (@top_three) {
			next if(!defined($peptide));
			my ($MH, $lsideMass, $rsideMass, $peptide_E_value, $Tag,
                   $TagEvalue, $TagRank, $reference, $peptide) =
                   @{$all_peptide_info->{$spout}->{$peptide}->{'data'}};
#				print $spout,"\t",$peptide,"\t",$peptide_Jscore->{$spout}->{$peptide},"\t",$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'},"\n";
			my $Jscore = $all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'};
			my $TagNum =  $all_peptide_info->{$spout}->{$peptide}->{'Num_of_Tags'};
			my @reference = @{$all_peptide_info->{$spout}->{$peptide}->{'ref'}};
			my $first_ref = shift @reference;

			my $part2_line1 = sprintf("%-2d  %-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-2d   %-3s    %-0.4f  %-20s    %-20s\n",$order,$MH,$lsideMass,$rsideMass,$peptide_E_value,$Tag,$TagRank,$TagNum,$Jscore,$first_ref,$peptide);
			push @newlines_part2, $part2_line1, map {"\t\t\t\t\t\t\t\t\t$_\n"} @reference;
			$order++;

			
		@{$newlines->{$spout}}=();	
		push @{$newlines->{$spout}}, @content_copy[0 .. $startIndex1],
               @newlines_part2, "\n", @content_copy[($startIndex2 - 2) .. $#content];
		}

	}
	
	return ($mainhash,$newlines);
}

sub get_miscleavage_mod_spout
{
	my ($self,$mainhash) = @_;
	my $panelty_mainhash;
	foreach my $spout (keys %$mainhash)
	{
		my $intpep = $mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"intpep"};
		my $mods = $self->get_mod_site_number($intpep);
		my $miscleavage = $self->get_mis_cleavage_number($intpep);
		if($mods or $miscleavage)
		{
			$panelty_mainhash->{$spout} = $mainhash->{$spout};
		}		
	}
	return $panelty_mainhash;
}

sub calculate_best_weight
{
	my ($self, $folder, $mainhash_all, $panelty_mainhash, $miscleavage_modification_coeff, $modification_coeff, $miscleavage_coeff) = @_;
	my ($weight, $best_weight, $best_sum_target) = (0,0,0);
	my $best_newline;
	my $mainhash;
#	my $mainhash_all = $self->parse_spOut_files_v5($folder);
	my ($sum_target,$sum_decoy,$cutoff_score) = $self->calculate_FDR($mainhash_all,0);
	my $newline = $self->parse_spOut_data($folder,$panelty_mainhash,$weight,$cutoff_score);
	my $phash = $self->get_db_mis_mod();
	my $mis_key = "mis" . $phash->{max_mis_cleavage};
	my $mod_key = "mod" . $phash->{max_modification};
	my $max_mis_space = $phash->{$mis_key};
	my $max_mod_space = $phash->{$mod_key};


	while($weight<2)
	{
		($mainhash,$newline) = $self->rerank_Jscore($newline,$weight,$phash);
		foreach my $spout (keys %$mainhash)
		{
			$mainhash_all->{$spout} = $mainhash->{$spout};
		}
		
		($sum_target,$sum_decoy,$cutoff_score) = $self->calculate_FDR($mainhash_all,0.01);
			
		if($sum_target > $best_sum_target)
		{
			$best_weight = $weight;
			$best_newline = $newline;
			$best_sum_target = $sum_target;
		}
		print $weight,"\t",$sum_target,"\t",$best_sum_target,"\n";		
		$weight += 0.5 / (int($max_mis_space+$max_mod_space-2) * 10);
	#	$weight += 0.1;
	}

	return $best_weight,$best_newline;
}




sub parse_spOut_data{
    my ($self,$folder,$panelty_mainhash,$weight,$cutoff_score)= @_;

    my $maxHitCondiered = 10;
    my $all_peptide_info;   # part 3 of input file
    my $newlines;           # ready-to-print new data of each file
	my $mainhash;
    my $peptide_Jscore;

	print "\n  Extracting PSMs that needs re-ranking\n";
    foreach my $spout (keys %$panelty_mainhash) {
        print "\r  Parsing file $spout                ";
		$spout .= ".spout";
        open(INPUT, "$folder/$spout") or die "Could open $spout. Error:$!\n";
        my @content = <INPUT>;
        my @content_copy = @content;
        close INPUT;
        $spout = basename($spout);
        $spout =~ s/\.spout//;

		my ($startIndex1, $startIndex2) = (0,0);
        if (scalar(@content) > 0) {
         
            # my $startIndex = first {$content[$_] =~ /-+\n/} 0 .. $#content;
            ($startIndex1, $startIndex2) =
                (grep {$content[$_] =~ /^-{2,}\n/} 0 .. $#content)[0, 1];

            for my $i ($startIndex1 + 1 .. $#content) 
			{
                last if ($content[$i] =~ /^\s*$/);
                last if ($content[$i] =~ /All identified peptide/);
                chomp($content[$i]); # Remove return line
                $content[$i] =~ s/^\s+//; # Remove leading spaces
				
				chomp($content[$i]); #Remove return line
				$content[$i]=~ s/^\s+//; #Remove leading spaces
				my @vals = split /\s+/,$content[$i];
				if(scalar(@vals)==1)
				{
					next;
				}

				elsif(scalar(@vals)==11)
				{								
					if($vals[0]!~/\d+/){
						next;
						#die "Error parsing file $spout at line $i\n The Oder field should be a numeric value\n";
					}
					
					$self->checkParameters($i,@vals);
					print $spout,$vals[8],"\n"  if ($vals[8]<0.1);
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"MH"}= $vals[1];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"lsideMass"}= $vals[2];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"rsideMass"}= $vals[3];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide_E_value"}= $vals[4];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag"} = $vals[5];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"TagRank"} = $vals[6];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag_num"} = $vals[7];				
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Weight_E_value"} = $vals[8];
								
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[9];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide"}= $vals[10];

				}	
            }
			next if($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"} > $cutoff_score);
			
            # parsing part 3 of input file
            # ============================

            my $peptide;
            for my $i (($startIndex2 + 1) .. $#content) {
                last if ($content[$i] =~ /^\s*$/);
                chomp($content[$i]); # Remove return line
                $content[$i] =~ s/^\s+//; # Remove leading spaces
                my @vals = split /\s+/, $content[$i];
                if ($#vals == 11) {
                    my @data = (my ($MH, $lsideMass, $rsideMass,
                        $peptide_E_value, $Tag, $TagEvalue, $TagRank,
                        $reference), $peptide) =
                        @vals[1, 2, 3, 4, 5, 6, 7, 10, 11];
						

					my $intpep = $peptide;
                    $intpep =~ s/[A-Z\-]\.([A-Z\@\#\*]+)\.[A-Z\-]/$1/;
						
                    my $mis_cleavage_num = $self->get_mis_cleavage_number($intpep);
					
                    my $mod_site_num = $self->get_mod_site_number($intpep);
########################################because multiple same #########################################
					$all_peptide_info->{$spout}->{$peptide}->{'Num_of_Tags'}++;

					my $Jscore = $peptide_E_value + $TagEvalue;

					$all_peptide_info->{$spout}->{$peptide}->{'data'} = \@data;
						
					if(!defined($all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'}))
					{

						$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'} = $Jscore;
						$peptide_Jscore->{$spout}->{$peptide} = $Jscore;						
					}
					elsif($Jscore>$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'})
					{

						$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'} = $Jscore;
						$peptide_Jscore->{$spout}->{$peptide} = $Jscore;
					}
		####### do not put the duplicates references ##################			
					unless($all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference})
					{					
						push @{$all_peptide_info->{$spout}->{$peptide}->{'ref'}}, $reference;
						$all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference} = 1;
					}
					
                } elsif (@vals == 1) {
                    my $reference = $vals[0];
		####### do not put the duplicates references ##################			
					unless($all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference})
					{
						push @{$all_peptide_info->{$spout}->{$peptide}->{'ref'}}, $reference;
						$all_peptide_info->{$spout}->{$peptide}->{'exist'}->{$reference} = 1;						
					}
                } 
				else {
                    printf STDERR "Warning: input line %d not parsed.\n", $i + 1;
                }
            }
		}
		
         my @top_three = (
                map {$_->[0]}
                sort {$b->[1] <=> $a->[1]}
                map {[$_, $peptide_Jscore->{$spout}->{$_}]}
                keys %{$peptide_Jscore->{$spout}})[0 .. 2];
				
        my @newlines_part2;
        my $order = 1;

        for my $peptide (@top_three) {
			next if (!defined($peptide));

            my ($MH, $lsideMass, $rsideMass, $peptide_E_value, $Tag,
                    $TagEvalue, $TagRank, $reference, $peptide) =
                    @{$all_peptide_info->{$spout}->{$peptide}->{'data'}};
#				print $spout,"\t",$peptide,"\t",$peptide_Jscore->{$spout}->{$peptide},"\t",$all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'},"\n";
			my $Jscore = $all_peptide_info->{$spout}->{$peptide}->{'highest_Jscore'};
			my $TagNum =  $all_peptide_info->{$spout}->{$peptide}->{'Num_of_Tags'};
            my @reference = @{$all_peptide_info->{$spout}->{$peptide}->{'ref'}};
            my $first_ref = shift @reference;

			my $part2_line1 = sprintf("%-2d  %-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-2d   %-3s    %-0.4f  %-20s    %-20s\n",$order,$MH,$lsideMass,$rsideMass,$peptide_E_value,$Tag,$TagRank,$TagNum,$Jscore,$first_ref,$peptide);
            push @newlines_part2, $part2_line1, map {"\t\t\t\t\t\t\t\t\t$_\n"} @reference;
            $order++;
        }
			
			
        push @{$newlines->{$spout}}, @content_copy[0 .. $startIndex1],
                @newlines_part2, "\n", @content_copy[($startIndex2 - 2) .. $#content];
        
    }
    return $newlines;
}


sub get_db_mis_mod
{
	my ($self) = @_; 
	my $param = $self->get_parameter();
	my $database = $param->{'database_name'};

	my $phash;
	$database =~ s/.mdx//;
	open(SDX,"$database.sdx") || die "can not open $database.sdx.\n";
	while(<SDX>)
	{
		my $line = $_;
		$line =~ s/\s*([;\#].+)$//;
     
		if($line =~ /^(.+?)\s*=\s*(.+)$/)
		{
			my ($key,$data) = ($1,$2);
			$data =~ s/\s+$//o;
			$phash->{$key} = $data;
      }     		
	}
	close(SDX);
	return $phash;	
}
=head
sub calculate_mass_shift
{
	my ($mass, $expmass, $charge) = ($$outhash{'MH'}, $$outhash{'expMH'}, $$outhash{'charge'});
	my $realdiff = $mass-$expmass; 
    my $ppm = ($diff/$mass/$charge)*(10**6);
	
    if ($$outhash{'XCorr'} >= $XCorr)
	{
		next if ($$outhash{'protein'} =~ /Decoy__/);
		%{$goodhash{$outfile}} = %$outhash;
	}
	my ($sum, $sum2, $n) = (0,0,scalar(keys %goodhash));
	while (my ($outfile, $hash) = each %goodhash)
	{
		my $ppm = $$hash{'ppm'};
		$sum += $ppm; 
		$sum2 += $ppm**2;
	}
	my $premean= $sum/$n;
	my $prestdev = $utils->calc_stdev($sum, $sum2, $n);	
    if ($n < $bf){
		while (my ($outfile, $outhash) = each %$runhash){
			$$outhash{status} = 1;
		}
        print "Skipping accurate mass filtering because of not enough good scans.\n";
        return;
    }	
	my @vgarray;
	for (my $cycle=1; $cycle<=2; $cycle++)
	{
		my @array;
		($sum, $sum2, $n) = (0,0,0);
		for my $outfile (sort {$goodhash{$a}{'scan'}<=>$goodhash{$b}{'scan'}} keys %goodhash)
		{
			my $ppm = $goodhash{$outfile}{'ppm'};
			if (($ppm < $premean - 2*$prestdev) || ($ppm > $premean + 2*$prestdev))
			{
				$$runhash{$outfile}{'status'} = -1;
				delete ($goodhash{$outfile});
			}
			else 
			{
				push (@array, $outfile);
				$sum += $ppm;
				$sum2 += $ppm**2;
				$n++;
			}
		}
		$premean = $sum/$n;
		$prestdev = $utils->calc_stdev($sum, $sum2, $n);
		@vgarray = @array;
	}
	my $ppmshift = $sum/$n;
	
	return $ppmshift;	
}
=cut

sub calculate_FDR
{
	my ($self,$mainhash,$FDR_cutoff)=@_;
	
	my $score = $self->score_distribution($mainhash);
	my $cutoff_score=0;
	my $best_ratio = 0;
	my %merge_score;
	foreach (keys %{$score->{"Decoy"}})
	{
			$merge_score{$_}=1;
	}

	foreach (keys %{$score->{"Target"}})
	{
			$merge_score{$_}=1;
	}

	my $sum_target = 0;
	my $sum_decoy = 0;

	foreach (reverse sort {$a<=>$b} keys %merge_score)
	{
			$score->{"Target"}->{$_}=0 if(!defined($score->{"Target"}->{$_}));
			$score->{"Decoy"}->{$_}=0 if(!defined($score->{"Decoy"}->{$_}));
			next if ($_==0);
			my $total = $sum_target;
			$total = 1 if $total == 0;
			my $FDR = ($sum_decoy)/($total);
			if($FDR>$FDR_cutoff)
			{		
				last;			
			}			
			$sum_target = $score->{"Target"}->{$_}+$sum_target;
			$sum_decoy = $score->{"Decoy"}->{$_}+$sum_decoy;
			$cutoff_score = $_;	
			next if ($score->{"Target"}->{$_} == 0);

	}
	return ($sum_target,$sum_decoy,$cutoff_score);
}

sub score_distribution
{
	my ($self,$mainhash) = @_;
	my $score;
	foreach my $file(keys %{$mainhash}){
		next if(!defined($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"}));
		my $pep_score = sprintf("%.01f",($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"}+0.002));
		
		if($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"reference"} =~ /Decoy/)
		{
	
			$score->{"Decoy"}->{$pep_score}++;
		}
		else
		{
			$score->{"Target"}->{$pep_score}++;			
		}
	}
	return $score;
}

sub checkParameters{
	my ($self,$i,@line) = @_;
	
	if($line[1]!~/\d+/ ){
		die "Error parsing file at line $i\n The MH field should be a numeric value\n";
	}
	
	if($line[2]!~/\d+/){
		die "Error parsing file at line $i\n The lSideMass field should be a numeric value\n";
	}
	
	if($line[3]!~/\d+/){
		die "Error parsing file  at line $i\n The rSideMass field should be a numeric value\n";
	}
	
	if($line[4]!~/\d+/){
		die "Error parsing file  at line $i\n The peptide E value  field should be a numeric value\n";
	}
}




1;

=head
my $tryptic = $utils->istryptic($peptide);

sub istryptic{  #returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
        shift @_;
        my ($Sequence) = @_;

        my $Nterm = 0;
        my $Cterm = 0;
        if ($Sequence =~ m/^[RK]\.[A-OQ-Z].*/ || $Sequence =~ m/^[\-]\..*/) {$Nterm = 1}; # Fixed bug DMD 3/20/2007
        if (($Sequence =~ m/.\.[A-Z\#\@\*]*[KR][\#\@\*]*\.[A-OQ-Z]\Z/) || ($Sequence =~ m/.\.[A-Z\#\@\*]*\.\-\Z/)) {$Cterm = 1};
        if ($Nterm && $Cterm) {
                return 3;
        } elsif ($Nterm || $Cterm){
                return 1 if $Nterm;
                return 2 if $Cterm;
        } else {
                return 0;
        }
}


=cut

