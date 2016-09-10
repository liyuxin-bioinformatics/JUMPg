#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Digestion
package Spiders::Digestion;
######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Coverter MzXML to dta files	          ****    #
#       ****					                      ****    #
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #
#       ****all rights reserved.		              ****    #
#       ****xusheng.wang@stjude.org		              ****    #
#       ****					                      ****    #
#       ****					                      ****    #
#       **************************************************    #
###############################################################

require Exporter;
use Spiders::DatabaseModif;
use Spiders::Params;

$VERSION='0.02';
@ISA=qw(Exporter);
@EXPORT=qw(enzDigest readFasta writeFasta writeFreq);

########### version 0.02 ####################################################
# To reduce the memory usage, save the masshash, peptidehash and proteinhash#
# on the disc. Keep the filepoint on the memory                             #
#                                                                           #
########### Version 12.1.0 ####################################################
# add side amino acid to peptide for example K.FGHLOPVR.L

sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
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


sub enzDigest{
	my ($self, $sqs, $paraHash) = @_;
	my $params = $self->get_parameter();
	my $enzyme = $paraHash->{'enzyme_info'};
	my $p = new Spiders::Params();
	my ($modif,$largest_modif) =$p->get_dynamic_modifications($params);	
	my $masshash;
	my $peptidehash;
	my $proteinhash;
	
##### $j is peptideID #######
	my $j = 0;

#### $i is proteinID #############
	for my $i (0..$#{$sqs}){
		my $id= $sqs->[$i]->{'id'};
		print "building index for protein: $id\n";
		my $desc = $sqs->[$i]->{'desc'};		
######### Version 12.1.0 ##############################

		$sqs->[$i]->{'seq'} = "-\." . $sqs->[$i]->{'seq'};
		$sqs->[$i]->{'seq'} = $sqs->[$i]->{'seq'} . "\.-";
######################################################
		
		my $peps = $self->digest($enzyme, $sqs->[$i]->{'seq'}, $paraHash->{'mis_cleavage'});
		
		foreach my $pepsequence (sort @{$peps})
		{
			next if (length($pepsequence)<7);
			my @pepseq=();
			push (@pepseq,$pepsequence);
			if ($params->{'digestion'} eq 'partial')
			{
				my $temp_pep = $self->generatePeps($pepsequence);
				push (@pepseq, @$temp_pep);
			}
		
			foreach $p (@pepseq)
			{			
				next if (!defined ($p));
	# changed by xusheng; remove the unknown amino acids			
				next if ($p =~ /X/);
				next if ($p =~ /U/);			
	#			$j++;
				next if(length($p)<7);
				
				$mw = $self->MW($p);
				next if ($mw < $paraHash->{'min_mass'} || $mw > $paraHash->{'max_mass'});
		################## partial tryptic ##############					

		#################################################
		
				$mw = $mw * 1000;
				$mw = int($mw); 
	#### remove these lines after add modified peptides 		
	#			push @{$masshash->{$mw}},$j;
			
	#			$peptidehash->{$j}->{'seq'}=$p;
	#			$peptidehash->{$j}->{'proteinID'}=$i;
	############################################################			
	############# generate modified peptide #####################			
				if((scalar keys %$modif)>0)
				{
					my $modif_masshash;
					undef $modif_masshash;
					my @temp_peptide;
					undef @temp_peptide;
	########## version 1.07 ######################
					my @pep_array = (split(/\./,$p));
					my $p = $pep_array[1];					
					push (@temp_peptide,$p);
					$modif_masshash->{$mw} = \@temp_peptide;

					$self->generate_modif_peptide($modif_masshash);
					
					foreach my $mass ( keys(%$modif_masshash) ) 
					{
				###### if the mass of modified peptide larger than 6000, next
						next if ($mass > ($paraHash->{'max_mass'} * 1000));	
						
						foreach my $pep (@{$modif_masshash->{$mass}})
						{
							my $mass_mw = int($mass); 
							$j++;
							push @{$masshash->{$mass_mw}},$j;
							$peptidehash->{$j}->{'seq'}="$pep_array[0]\.$pep\.$pep_array[2]";
							$peptidehash->{$j}->{'proteinID'}=$i;
						}
						undef @{$modif_masshash->{$mass}}; 						
					}
					undef %$modif_masshash;
					
				}
	################# no dynamic modif ##########################			
				else
				{
					$j++;
					push @{$masshash->{$mw}},$j;
					$peptidehash->{$j}->{'seq'}=$p;
					$peptidehash->{$j}->{'proteinID'}=$i;				
				}
	#########################################################
			}	
		}
		undef $peps;
		$proteinhash->{$i}->{'id'} = $id;
		$proteinhash->{$i}->{'desc'} = $desc;
	}
	return ($masshash,$peptidehash,$proteinhash);
	
}

sub generate_modif_peptide
{
	my ($self,$pephash) = @_;
	my $cmdatabase = new Spiders::DatabaseModif();
#	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");
	my %mod_symbol = ('M'=>"@",'S'=>"#",'T'=>"%",'Y'=>"*",'G'=>"^",'K'=>"&",'D'=>'?','A'=>'~','Q'=>'!','P'=>"(",'E'=>")",'E'=>"{",'V'=>"}","V"=>"[","H"=>"]","C"=>":","F"=>",","I"=>';',"L"=>',',"R"=>"<","N"=>">","W"=>"'");	
	my $params = $self->get_parameter();
#	my $combNUM = $params->{'max_modif_num'};
	my $p = new Spiders::Params();
	my ($modif,$largest_modif) =$p->get_dynamic_modifications($params);

	my %new_modif=();
	my $i=0;
	foreach $modif_aa (sort keys %$modif)
	{
#		$new_modif{$modif_aa} = $modif_aa . $mod_symbol[$i];
		$new_modif{$modif_aa} = $modif_aa . $mod_symbol{$modif_aa};
		$i++;		
	}
	$cmdatabase->GenerateModifSeq($modif,$params,$pephash,\%new_modif);
}

sub rmRedundancy{

	my ($self, $ps) = @_;
	my $h = {};
	foreach (@{$ps}){
		$h->{$_->{'id'}} = $_;
	}
	
	my $hh = [];
	foreach (sort keys %$h){
		push @{$hh},$h->{$_};
	}
	return $hh;
}


sub writeFreq{
	my ($self, $f, $freq) = @_;
	open OUT, '>'.$f or die;
	print OUT "Mass\tFrequency\n";
	foreach (keys %$freq){
		print OUT $_,"\t",$freq->{$_},"\n";
	}
	close OUT;
}



sub digest{

	my ($self, $enzyme, $seq, $NumOfMisCleavage) = @_;
	my @enzyme_array = split(/\s+/,$enzyme);
	my @sites = split(/\s*/,$enzyme_array[1]);
	my $uncleavage = $enzyme_array[2];
####################3 skip KP ##############
	for my $i (0..$#sites){
		$seq =~ s/($sites[$i])(?!$uncleavage)/$1 /g;
	}
	$seq =~ s/ \Z//;
	
	my @a = split(/ /,$seq);	
	my @pep;
	for(my $i=0;$i<=$#a;$i++)
	{

		my $frag="";
		if($i==0)
		{
			$frag = $a[$i] . "\.";
			$frag .= substr($a[$i+1],0,1); 	
		}
		elsif($i==$#a)
		{
			$frag = "\." . $a[$i];
			next if $a[$i] eq ".-";
			$frag = substr($a[$i-1],length($a[$i-1])-1,1) . $frag;					
		}
		else
		{
			$frag = substr($a[$i-1],length($a[$i-1])-1,1) . "\.$a[$i]";
			if($i==($#a-1) and $a[$i+1] eq ".-")
			{
				$frag = "$frag" . $a[$i+1];
			}
			else
			{
				$frag = "$frag\." . substr($a[$i+1],0,1);
			}
		}

		push(@pep,$frag);
	}
	undef @a;
	@a = @pep;
	undef @pep;
	
	if ($#a > -1){

		my @f = @a;
		
		for my $i (0..$#a){		
			my $ub = ($#a > $i+$NumOfMisCleavage) ? $i+$NumOfMisCleavage : $#a;
			for my $j ($i+1..$ub){
				my @b = @a;
				my @comb = splice(@b,$i,$j-$i+1);
				my $comb_pep = "";
				for(my $k=0;$k<=$#comb;$k++)
				{

					my @data = split(/\./,$comb[$k]);
					if($k==0)
					{
						$comb_pep = "$data[0]\." . $data[1];
					}
					elsif($k==$#comb)
					{
						$comb_pep .= "$data[1]\.$data[2]";
					}
					else
					{
						$comb_pep .= $data[1];
					}
				}
				
				push @f, $comb_pep;
				 
			}
		}
		my $h = {};
		my $flag = 1;
## Version 12.1.0		
#		foreach my $frag (sort @f){
		for(my $i=0;$i<=$#f;$i++)
		{
			my $frag = $f[$i];

			for my $i (0..$#sites){
				my $tmpFrag = $frag;
				$tmpFrag = (split(/\./,$tmpFrag))[1];
				my $site_uncleavage = $sites[$i] . $uncleavage;
				my $lc_site_uncleavage = lc($site_uncleavage);
				$tmpFrag =~ s/$site_uncleavage/$lc_site_uncleavage/;
				my $num = grep /$sites[$i]/,split(//,$tmpFrag);
			#	print $frag,"\t",$num,"\n";
				$flag *= ($num > $NumOfMisCleavage+1) ? 0 : 1;
			}
			
			if ($flag){
		#		print $frag,"\n";
				$h ->{$frag}++;
			}
			$flag = 1;
		}
		
		return [sort keys %$h];
	}else{
		return [@a];
	}
}

sub MW{

	my ($self, $sq) = @_;
	my $sq = (split(/\./,$sq))[1];
	my $parameter = $self->get_parameter();
	
	my $aaMW = {	
				'G' => 57.02146372057,	
				'D' => 115.02694302383,
				'A' => 71.03711378471,
				'Q' => 128.05857750528,
				'S' => 87.03202840427,
				'K' => 128.094963014,
				'P' => 97.05276384885,
				'E' => 129.04259308797,
				'V' => 99.06841391299,
				'M' => 131.04048491299,
				'T' => 101.04767846841,
				'H' => 137.05891185845,
				'C' => 103.00918478471,
				'F' => 147.06841391299,
				'I' => 113.08406397713,
				'L' => 113.08406397713,
				'R' => 156.1011110236,
				'N' => 114.04292744114,	
				'Y' => 163.06332853255,
				'W' => 186.07931294986,				
				};
### add the static modification ################## 		
	foreach my $param (keys %$parameter)
	{
		if($param =~/add\_(\w+)\_/ and $parameter->{$param}>0)
		{
			$aaMW->{$1} += $parameter->{$param};
		}
	}	
			
	my $h = {};
	my @aas = split(//,$sq);
	foreach (@aas){
		$h->{$_}++;
	}
	my $mw = 0;
	if($parameter->{"add_Nterm_peptide"}>0)
	{
		$mw += $parameter->{"add_Nterm_peptide"};
	}
	if($parameter->{"add_Cterm_peptide"}>0)
	{
		$mw += $parameter->{"add_Cterm_peptide"};
	}
	
	foreach $aa (keys %$h){
		if(!defined($aaMW->{$aa}))
		{
#			print "\ryour sequences contain unkown amino acid: $aa!\n";
			$mw=0;
			last;
		}
		$mw += $h->{$aa} * $aaMW->{$aa};
	}

	### added the H20 and H+
	$mw += 19.017841150512;
#	print $sq,"\t",$mw,"\n";
	return $mw;
}

sub readFasta{

	my ($self, $f) = @_;
	
	my $parsedSqs = [];
	open F, $f or die;
	my @lines = <F>;
	close F;
	
	my @sqs = split("\n\>",join('',@lines));
	$sqs[0] =~ s/^\>//;
	
	
	for my $i (0..$#sqs){
		my $h = {};
		
		my @a = split("\n",$sqs[$i]);

		
		$h->{'id'} = substr($a[0], 0, index($a[0],' '));
		$h->{'desc'} = substr($a[0], index($a[0],' ')+1);
		shift @a;
		$h->{'seq'} = join('',@a);
		push @{$parsedSqs}, $h;
	}
	return $parsedSqs;
}

sub writeFasta{
	my ($self, $f, $sqs) = @_;
	open OUT, '>'.$f or die;
	for my $i (0..$#{$sqs}){
		print OUT '>'.$sqs->[$i]->{'id'},' ',$sqs->[$i]->{'desc'},"\n";
		print OUT $sqs->[$i]->{'seq'},"\n";
	}
	close OUT;
}

sub generatePeps{

	my ($self, $sq) = @_;
######### Have not been fully tested ############
	my ($start,$sq,$end) = (split(/\./,$sq));	

	my @a = split(//,$sq);
	my @b = split(//,$sq);
	
	my @g;
	for my $i (0..$#a){
		my @t= @a;
		$a[$i+1] = $end if($i==$#a);
		my $p_frag = "$start\." . join('',splice(@t,0,$i+1)) . "\.$a[$i+1]";
		push @g, $p_frag;
	}
	for( my $i=$#b; $i>=0; $i--){
		my @t= @b;

		my $p_frag = "$b[$i]\." . join('',splice(@t,$i+1,$#b)) . "\.$end";	
#		push @g, join('',splice(@t,0,$i+1));
		push @g, $p_frag;
	}
	
	return [@g];

}






#############################################

1;

__END__
