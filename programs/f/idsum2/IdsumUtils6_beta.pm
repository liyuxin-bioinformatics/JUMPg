#!/usr/bin/perl -wT

## Release date: 11/01/2015
## Release version: version 12.1.0

package idsum2::IdsumUtils6_beta;
use strict;
use idsum2::CommonUtils;
use Data::Dumper;
use Storable;
use idsum2::CL;
use Cwd;
use FileHandle;
use idsum2::pepXML_parser;
use POSIX;
#use Clone qw(clone);
#use File::Sync qw(fsync sync);

###### added by yanji
use File::Basename;
use idsum2::XMLParser;
use Statistics::R;

my $xml = idsum2::XMLParser->new();

#my $log_file = "idsum.log";
my $log_file = "jump_f.log";
###### end of addition

my $utils = idsum2::CommonUtils->new();
my $H = 1.00728;
my %dixon; $utils->get_critical_values(\%dixon);
my ($M0,$M1,$M2,$M3,$M4,$MX) = (0,0,0,0,0,0);

###### added by yanji
#my $max_mis_del = 0;
#my $max_mod_del = 0;

sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);
                                                                                                                                                             
  return $self;
}

#####################################################################
# PARSING FUNCTIONS 																								#
#####################################################################
sub parse_params{
  shift @_;
	my ($paramhash, $paramfile, $listarray) = @_;
	
	my @combinations = ("FT1","FT2","FT3","FT4","PT1","PT2","PT3","PT4","NT1","NT2","NT3","NT4");
	my @masstypes = ("MH","MH, MH+1","MH, MH+1, MH+2","MH, MH+1, MH+2, MH+3","MH, MH+1, MH+2, MH+3, MH+4", "M-1, MH, MH+1, MH+2", "MH-2, MH-1, MH, MH+1, MH+2");
	open (IN, "<$paramfile");
	my ($group, $groupsetup, $grouptxt) = ("",0,0);
	my $groupnum = 0;
	while (<IN>)
	{
		next if (/^\#/ || /^\s*\Z/ || /^\s*\#\.*/);
		#next if (/^\#/ or /^\s*\Z/);

		if (/\=/)
		{
			if (/\#/){
				$_ =~ s/\#.*//g;
			}
			$_ =~ s/\s+//g;
			my ($key, $value)  = split("\=", $_);
			if ($key =~ /12combinations/){
				my @bools = split("", $value);
				my $num = 0;
				for my $entry (@bools){
					$$paramhash{$key}{$combinations[$num]} = $entry;
					$num++;
				}
			} else {
				if ($key =~ /mass_consideration/){
					$$paramhash{$key}{'value'} = $value;
					$$paramhash{$key}{'type'} = $masstypes[$value-1];
				} else {
					$$paramhash{$key} = $value;
				}
			}
		}
		else 
		{ # Groups
			#next if (/^\s*\Z/);
			#next if (/^\#/ || /^\s*\Z/ || /^\s*\#\.*/);
			$_ =~ s/\s+//g;
			if (/\:/){
				($group, my $member)  = split(/\:/, $_);
				$groupnum++;	$$paramhash{'grouporder'}{$group} = $groupnum;
				push (@$listarray, $group) if (defined($listarray));
				$member =~ s/\/$//;		# rm ending '/' for input path
				push (@{$$paramhash{'groups'}{$group}}, $member);
			} else {
				s/\/$//;			# rm ending '/' for input path
				push (@{$$paramhash{'groups'}{$group}}, $_);
			}
		}
	}
	close IN;

	#for new version of idsum2.params hard Code params
	hardCode_params($paramhash);
}

sub old_parse_params{
  shift @_;
	my ($paramhash, $paramfile, $listarray) = @_;
	
	my @combinations = ("FT1","FT2","FT3","FT4","PT1","PT2","PT3","PT4","NT1","NT2","NT3","NT4");
	my @masstypes = ("MH","MH, MH+1","MH, MH+1, MH+2","MH, MH+1, MH+2, MH+3","MH, MH+1, MH+2, MH+3, MH+4", "M-1, MH, MH+1, MH+2", "MH-2, MH-1, MH, MH+1, MH+2");
	open (IN, "<$paramfile");
	my ($group, $groupsetup, $grouptxt) = ("",0,0);
	my $groupnum = 0;
	while (<IN>)
	{
		$groupsetup = 0 if (/^\#\sMinimum/);
		if ($groupsetup == 0)
		{
			if (/Grouping/){
				$groupsetup = 1; next;
			}
			next if (/^\#/ || /^\s*\Z/ || /^\s*\#\.*/);
			if (/\#/){
				$_ =~ s/\#.*//g;
			}
			$_ =~ s/\s+//g;
			my ($key, $value)  = split("\=", $_);
			if ($key =~ /12combinations/){
				my @bools = split("", $value);
				my $num = 0;
				for my $entry (@bools){
					$$paramhash{$key}{$combinations[$num]} = $entry;
					$num++;
				}
			} else {
				if ($key =~ /mass_consideration/){
					$$paramhash{$key}{'value'} = $value;
					$$paramhash{$key}{'type'} = $masstypes[$value-1];
				} else {
					$$paramhash{$key} = $value;
				}
			}
		} 
		else 
		{ # Groups
			next if (/^\s*\Z/);
			next if (/^\#/ || /^\s*\Z/ || /^\s*\#\.*/);
			$_ =~ s/\s+//g;
			$grouptxt = 1 if (/\.txt\Z/);
			if (/\:/){
				($group, my $member)  = split(/\:/, $_);
				$groupnum++;	$$paramhash{'grouporder'}{$group} = $groupnum;
				push (@$listarray, $group) if (defined($listarray));
				push (@{$$paramhash{'groups'}{$group}}, $member);
			} else {
				push (@{$$paramhash{'groups'}{$group}}, $_);
			}
		}
	}
	close IN;
	$$paramhash{'grouptxt'} = $grouptxt;

	#for new version of idsum2.params hard Code params
	hardCode_params($paramhash);
}

sub hardCode_params
{
	my ($paramhash)=@_;

	if (!defined($$paramhash{output_pepXML})) {$$paramhash{output_pepXML}=0;}
	$$paramhash{MS_Source}=1;
	$$paramhash{add_redundancy}=0;
	$$paramhash{combine_filtermethods}=0;
	$$paramhash{norm_XCorr}=0;
	$$paramhash{corereport_ID}=1;
	$$paramhash{corereport_SC_comp1}=0;
	if (defined($$paramhash{initial_outfile_fdr})) { $$paramhash{min_peptide_fpr}=$$paramhash{initial_outfile_fdr}; }
	else {$$paramhash{min_peptide_fpr}=5;}
	#$$paramhash{min_peptide_fpr}=50;

	$$paramhash{trypticity}=2;
	$$paramhash{charge}=2;
=head
	if (!defined($$paramhash{XCorr}))
	{
		if ( $$paramhash{search_engine} eq 'jump' ) { $$paramhash{XCorr}=20; }
		else {	$$paramhash{XCorr}=2.5;}
	}
=cut
	$$paramhash{dCn}=0.00;
	$$paramhash{scan_events}=4;

	$$paramhash{keep_decoy}=1;
	$$paramhash{core_report}=0;
	$$paramhash{publication_table}=0;
	$$paramhash{bypass_grouping}=1;
	$$paramhash{subgroup}=50;
	$$paramhash{abundance_index}=0;
	$$paramhash{PAI_top}='TP';
	$$paramhash{max_peptide_length}=30;
	$$paramhash{min_peptide_hydro}= -24;
	$$paramhash{max_peptide_hydro}=17;
	$$paramhash{enzyme}='trypsin';
}

sub parse_badproteinfile{
	shift @_;
	my ($proteinhash, $file) = @_;
	
	open (IN, "<$file");
	while (<IN>){
		chomp;
		my $protein = $_;
		$protein =~ s/(^[A-Za-z0-9\.\_\:\-\|]+)\s//; $protein = $1;
		$$proteinhash{$protein}++;
		if ($protein =~ /\|/){
			my @proparts = split('\|', $protein);
			my $shortpro = $proparts[scalar(@proparts)-1];
			$$proteinhash{$shortpro}++;
		}
	}
	close IN;
}

sub parse_jumpparams
{
	shift @_;
	my ($seqparamhash, $fraction) = @_;
	my $paramfile = "$fraction\/jump.params";
	unless (-e $paramfile) { $paramfile = 'jump.params'; }

	while (!(-e $paramfile))
	{
		die "\n$paramfile does not exists!!!\n\n";
	}

	open (IN, "<$paramfile");
	while (<IN>)
	{
		next if (/^#/  || /^\s+\Z/);
		chomp; s/[;#].*//; s/\s+\Z//;
		my ($key, $value) = split(" = ", $_);
		if ($key eq 'search_engine')
		{
			$$seqparamhash{$key}=$value;
		}
		elsif ($key eq 'pit_file')
		{
			if ( !defined($$seqparamhash{$key}) or $$seqparamhash{$key} eq '0' )
			{ $$seqparamhash{$key}=$value; }
		}
	}
	close IN;
}

sub parse_seqparams{
	shift @_;
  my ($seqparamhash, $fraction) = @_;

###### added by yanji
  open (LOGFILE, ">>$log_file");
  #select((select(LOGFILE), $|=1)[0]); 
  my $paramfile = "$fraction\/sequest.params";
  unless (-e $paramfile) { $paramfile = 'sequest.params'; }
                                                                                                                                                             
  while (!(-e $paramfile)){
    print "\n$paramfile does not exists!!!";
    print LOGFILE "\n$paramfile does not exists!!!";
    $paramfile = idsum2::CL->ask("\nPlease enter the correct full path to the sequest parameter file.", "");
  }
  open (IN, "<$paramfile");
  while (<IN>){
    last if (/SEQUEST_ENZYME_INFO/);
    next if (/^#/ || /\[SEQUEST\]/ || /^\s+\Z/);
    chomp; s/[;#].*//; s/\s+\Z//;
    my ($key, $value) = split(" = ", $_);
    next if (!defined($value));
    if ($key =~ /^diff_search_options/){
      my @array = split('\s', $value);
      for (my $i=0;$i<scalar(@array);$i+=2){
        my ($aa, $diff) = ($array[$i+1], $array[$i]);
				next if ($diff =~ /^[0\.]+\Z/);
        $$seqparamhash{$fraction}{"dynamicmods"}{$aa} = $diff;
      }
    } elsif ($key =~ /^add/){
			next if ($value == 0);
      if ($key =~ /terminus/ || /term/){
        if ($key =~ /add_N_terminus/ || /add_Nterm_peptide/){
          $$seqparamhash{$fraction}{"staticmods"}{'nterm'} = $value;
        } else { # Cterm mod
          $$seqparamhash{$fraction}{"staticmods"}{'cterm'} = $value;
        }
      } else {
      	$key =~ s/add_([A-Z])_.*/$1/;
      	$$seqparamhash{$fraction}{"staticmods"}{$key} = $value;
			}
    }
  }
  close IN;
}

sub get_dCn {
	shift @_;
	(*IN) = @_;

	seek(IN,0,0);
  my ($dCntemp, $dCntemp2, $dCntemp3, $dCntemp4) = grep(/^\s+2\.\s+\d+\s+\// || /^\s+3\.\s+\d+\s+\// || /^\s+4\.\s+\d+\s+\// || /^\s+5\.\s+\d+\s+\//, <IN>);
	$dCntemp =~ s/\s*\/\s*/\//g; $dCntemp =~ s/(\s)\s+/$1/g; $dCntemp =~ s/^\s+//; $dCntemp =~ s/\s+\Z//;
########## new format ################
  my ($z, $y, $i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp);
	if ($dCn <= 0.01){
		$dCntemp2 =~ s/\s*\/\s*/\//g; $dCntemp2 =~ s/(\s)\s+/$1/g; $dCntemp2 =~ s/^\s+//; $dCntemp2 =~ s/\s+\Z//;
		($z, $y, $i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp2);
		if (defined($dCntemp3) && $dCn <= 0.01){
			$dCntemp3 =~ s/\s*\/\s*/\//g; $dCntemp3 =~ s/(\s)\s+/$1/g; $dCntemp3 =~ s/^\s+//; $dCntemp3 =~ s/\s+\Z//;
			($z, $y, $i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp3);
		}
		if (defined($dCntemp4) && $dCn <= 0.01){
			$dCntemp4 =~ s/\s*\/\s*/\//g; $dCntemp4 =~ s/(\s)\s+/$1/g; $dCntemp4 =~ s/^\s+//; $dCntemp4 =~ s/\s+\Z//;
			($z, $y, $i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp4);
		}
	}
	return $dCn;
}

sub get_intprored{
	shift @_;
	(*IN) = shift @_; 
	my ($hash, $protein, $red) = @_;
	
	seek(IN,0,0);
	my @lines = <IN>;	splice(@lines, -6);
	my @proteinlines = grep(/^\s\s\s\s\s\s[0-9A-Za-z]/, @lines, $red);
	my $num = 0;
	for my $pro (@proteinlines){
		$num++;	$pro =~ s/^\s+([a-zA-Z0-9\.\_\-\|\:]+)[\s\,\;]*.*//; $pro = $1;	$$hash{$pro} = 1;
		last if ($num == $red || $num == 20);
	}

####changed by xusheng ##############
#
#	if ($protein =~ /^(Decoy__|Random__)/i){
     if ($protein =~ /(Decoy__|Random__)/i){
		for my $pro (keys %$protein){

#			next if (/^(Decoy__|Random__)/i);
			next if (/(Decoy__|Random__)/i);
		}
		print Dumper(\@lines);exit;
	} else {
		
	}
	#print Dumper($hash);exit;
}

sub establish_CPR_CPN_hashes
{
        my ($self,$in,$cpr,$cpn)=@_;
	my $pth=getcwd();#print "$pth\/$in\n";
        open(IN,"<$pth\/$in") or die "Cannot open $pth\/$in!!!\n";
        my $line=<IN>;
        while(<IN>)
        {
                my @tmp=split(/\t/,$_);
                if (int($tmp[2])>0) { $$cpr{$tmp[0]}=$tmp[2]*1000/$tmp[1]; }
                else { $$cpn{$tmp[0]}=0;}# print "$tmp[0]";}
        }
        close IN;
}

sub establish_establish_rnaseqReads_hash
{
        my ($self,$in,$hash)=@_;
	my $pth=getcwd();#print "$pth\/$in\n";
        open(IN,"<$pth\/$in") or die "Cannot open $pth\/$in!!!\n";
        while(<IN>)
        {
                my @tmp=split(/\t/,$_);
                $$hash{$tmp[1]}=$tmp[2];
        }
        close IN;
}

sub establish_pep2reads
{
	my ($self,$in,$hash)=@_;
	my $pth=getcwd();
	open(IN,"<$pth\/$in") or die "Cannot open $pth\/$in!!!\n";
	while(<IN>)
	{
		chomp;
		my @tmp=split(/\t/,$_);
		$$hash{$tmp[0]}=$tmp[1];
	}
	close IN;
}

sub _update_max_levels
{
        my ($mxlv,$mxtd,$mxG,$mxRpk,$mxProName,$mxRank,$pr,$cpr,$cpn,$i)=@_;

        my ($lv,$td)=protein_priorty_level($pr,$cpr,$cpn);
	#if ( $td eq 'D' ) { $lv+=0.5; } #decoy penalty

	# mx level change
        if ( $$mxlv>$lv ) 
	{ 
		$$mxlv=$lv; $$mxtd=$td; $$mxProName=$pr; $$mxRank=$i;

		# for specific levels (1 & 4)
		if ($lv==1)
		{
			$$mxRpk=get_RNAseq_value($pr,$cpr);
			$$mxG=$pr; 
		}
	}
	elsif ( $$mxlv==$lv )
	{
		# for same level, check for quantitive variables (rpk, readNumber)
		if ($lv==1)
		{
			my $rpk=get_RNAseq_value($pr,$cpr);
			if ( $rpk>$$mxRpk )
			{
				$$mxProName=$pr; $$mxRpk=$rpk; $$mxG=$pr; $$mxRank=$i;
			}
		}
	}
}

sub protein_priorty_level
{
        my ($pr,$cpr,$cpn)=@_;
        my $lv=0; my $td='';

        if ( $pr =~ m/Decoy__/) { $pr =~ s/\#\#Decoy__//; $td='D';}
        else { $td='T'; }

        if (defined($$cpr{$pr})) {$lv=1;}
        elsif ( $pr =~ m/Con$/ or $pr =~ /^co\|CON/ ) {$lv=2;}
        elsif (defined($$cpn{$pr})) {$lv=3;}
	#4_1102_5476_147340-1|F0|45-87
        #elsif ( $pr =~ /[0-9\_]+\-\d\|[FR]\d\|\d+-\d+$/ )  {$lv=4;}
        elsif ( $pr =~ /[0-9\_]+\-\d\|[FR]\d\|/ )  {$lv=4;}
	else { die "unclassified protein names: $pr !!!"; }

        return ($lv,$td);
}

sub get_RNAseq_value
{
        my ($pr,$cpr)=@_;

        $pr =~ s/\#\#Decoy__//;
        if (defined($$cpr{$pr})) {return $$cpr{$pr};}
        else { die "not exist protein ($pr) in CPR!!!\n"; }
}

sub get_reads_counts
{
        my ($pr,$cpr)=@_;

        $pr =~ s/\#\#Decoy__//;
        if (defined($$cpr{$pr})) {return $$cpr{$pr};}
        else { die "not exist seq ID ($pr) in pep ~ reads counts!!!\n"; }
}

sub peptide_similarity
{
	my ($oldPep,$peptide)=@_;

	# rm flanking AA
	#$oldPep =~ s/^(\w\.)//; $oldPep =~ s/(\.\w)$//;
	#$peptide =~ s/^(\w\.)//; $peptide =~ s/(\.\w)$//;
	$oldPep =~ s/[A-Z\-]\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\.[A-Z\-]/$1/;
	$peptide =~ s/[A-Z\-]\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\.[A-Z\-]/$1/;

	# rm mods
	$oldPep =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;
	$peptide =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;

	# substitute AA with the same mass
	$oldPep =~ s/I/L/g; $peptide =~ s/I/L/g;

	# checking starts
#=head
	# method 1 (Junmin's method): composition based
	my %mt;
	my @a=split(//,$oldPep); my @b=split(//,$peptide);
	foreach my $x (@a) { if (defined($mt{$x}{'1'})) {} else { $mt{$x}{'1'}=$mt{$x}{'2'}=0; } $mt{$x}{'1'}++; }
	foreach my $x (@b) { if (defined($mt{$x}{'2'})) {} else { $mt{$x}{'1'}=$mt{$x}{'2'}=0; }$mt{$x}{'2'}++;}
	my $k=0; foreach my $x (keys %mt) { if ( $mt{$x}{1}<$mt{$x}{2} ) { $k+=$mt{$x}{1}; }else{ $k+=$mt{$x}{2}; } }
	#if ( $k*2/(scalar(@a)+scalar(@b))>0.9 and scalar(@a)!=scalar(@b)) { print "$oldPep,$peptide\n"; }
	return $k*2/(scalar(@a)+scalar(@b));
#=cut

=head
	# method 2 (Yuxin's method): should be equal length; positional check
	if (length($oldPep) != length($peptide)) { return -1; }
	my @a=split(//,$oldPep); my @b=split(//,$peptide); my $k=0;
	for ( my $i=0; $i<scalar(@a); $i++ ) { if ($a[$i] eq $b[$i]) {$k++;} }
	return $k/scalar(@a);
=cut
}

sub orig_pep_for_decoy
{
	my ($paramhash,$peptide)=@_;
	if ( $$paramhash{decoy_strategy} eq 'RS3' )
	{
		$peptide =~ s/[\@\%\&\^\~\$\#\*\^\~\$\.]//g;
		#my $enzyme_info=$$paramhash{'enzyme_info'};
		my $enzyme_info='Trypsin 1 1 KR P';
		my @eArray=split(/\s+/,$enzyme_info);
		my @sites = split(/\s*/,$eArray[3]);
		my $uncleavage = $eArray[4];

		my $seq=$peptide;
		my @nt=split(//,$seq);

		for (my $i=$#nt; $i>=0; $i--)
		{
			for my $j (0..$#sites)
			{
				if ($nt[$i] eq $sites[$j])
				{
					my $pre;
					if ($i<$#nt) { $pre=$i+1;} else { $pre=0;}
					$nt[$i]=$nt[$pre];
					$nt[$pre]=$sites[$j];
					last;
				}
			}
		}
		return reverse(join('',@nt));
	}
	else { die "undefined decoy strategy: $$paramhash{decoy_strategy}!!!\n"; }
}

sub get_read_spt
{
	my ($peptide,$protein,$paramhash,$pep2reads)=@_;

	my $origPep;
	if ( $protein =~ /Decoy/ ){ $origPep=orig_pep_for_decoy($paramhash,$peptide);}
	else { $origPep=$peptide; $origPep =~ s/[\@\%\&\^\~\$\#\*\^\~\$\.]//g;}

	my @aa=split(//,$origPep); pop(@aa); shift(@aa);
	$origPep=join('',@aa);
	if (defined($$pep2reads{$origPep})) { return $$pep2reads{$origPep};}
	else {return 0;}
}

sub check_mod_pep
{
	my ($peptide,$mod)=@_;

	#open(OUT,">>mod_pepride_annotation.txt");

	my @mods=split(//,$mod);
	my $withMod=0;
	foreach my $m (@mods)
	{
		if ($peptide =~ m/$m[\@\%\&\^\~\$\#\*\^\~]/) { $withMod=1; last; }
	}

	#print OUT "$peptide\t";
	#if ($withMod) { print OUT "mod\n"; } else { print OUT "nomod\n"; }
	#close OUT;

	return $withMod;
}

sub parse_outfile
{
  shift @_;
  my ($paramhash,$peptidehash,$runhash,$run,$outfile,$delhash,$blank,$cpr,$cpn,$pep2reads) = @_;
  my %proteins;
  my $fullpath = "$run\/$outfile";

  # check if RNAseq data is vailable
  my $noRNAseq;
  if (!defined($cpr) || !defined($cpn) || !defined($pep2reads)) { $noRNAseq=1; } else { $noRNAseq=0; }
 
  open (IN, $fullpath);
  #$outfile =~ m/\.(\d+)\.\d+\.(\d)\.out$/;
  #my ($scan,$charge) = ($1,$2);
  $outfile =~ m/\.(\d+)\.(\d+)\.(\d)\.out$/;
  my ($scan,$charge) = ("$1.$2",$3);
  my ($sorcerer) =  grep(/SageNResearch/ || /SORCERER/, <IN>);
  if (defined($sorcerer)){ $$paramhash{'sorcererouts'} = 1; }

  #get PSM lines
  seek(IN, 0,0);
  my @hits=grep(/^\s+\d+\.\s+\d+\s+\//,<IN>);
  if (scalar(@hits)==0) { $$blank++; return; }

  #get red protain names
  seek(IN,0,0);
  my $red_count=0;
  my @redProNames=grep(/^\s+\d+\s+[0-9A-Za-z#]/,<IN>);
  close IN;

  #decide which PSM and protein name to use
  my ($mxlv,$mxtd,$mxG,$mxRpk,$dCn,$mxProName,$mxRank,$mxReadSpt)=(10,'','',0,0,'',0,0,0);
  my @tmpProteins; my $oldPep='';
  #for (my $i=0; $i<scalar(@hits); $i++)
  for (my $i=0; $i<scalar(@hits) and $i<10; $i++)
  {
    #remove redundant space
    $hits[$i] =~ s/\s*\/\s*/\//g; $hits[$i] =~ s/(\s)\s+/$1/g; $hits[$i] =~ s/^\s+//; $hits[$i] =~ s/\s+\Z//;
    my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $hits[$i]);
    $dCn=$c;
    my ($peptide, $red)  = ($g, 0);
    if(defined($h)) { $red = $g; $red =~ s/\+//; $peptide = $h; }
    #last if ( $i>=1 and $dCn>0 ); # old critera

    last if ( $i>=1 and $dCn>0.05 ); 
    my $pepSim=1;
    if ( $oldPep ne '' ) { $pepSim=peptide_similarity($oldPep,$peptide); }
	#if ( $pepSim<0 ) {die "uneuqla length for peptide sequences with small dCn!!!\n $outfile,$oldPep,$peptide\n";}
    last if ( $pepSim<0.9 ); 
    last if ( $i>=10 );

    # for 1st protein name
    if ($noRNAseq) { if ($i==0) { $mxProName=$protein; } }
    else {_update_max_levels(\$mxlv,\$mxtd,\$mxG,\$mxRpk,\$mxProName,\$mxRank,$protein,$cpr,$cpn,$i);}

    # for redundant protein names
    my $j;
    for ( $j=1; $j<=$red; $j++)
    {
	if (!defined($redProNames[$red_count])) {die "undefined red proten names lines:\n $outfile,$i,$j,$red,$red_count\n";}
      chomp($redProNames[$red_count]); my $tmpP=$redProNames[$red_count]; $red_count++;
      #$tmpP =~ s/^\s+\d+\s+//;  if ( $tmpP =~ /\s/ ) { $tmpP =~ s/^(.*?)\s+/$1/; }   $tmpProteins[$i][$j]=$tmpP;
      if($tmpP=~/^\s+\d*\s*([\#a-zA-Z0-9\.\_\-\|\:]+)\s*/) { $tmpProteins[$i][$j]=$1;} else {die "Ero:red pro\n";}
      if ($noRNAseq) {}
      else {_update_max_levels(\$mxlv,\$mxtd,\$mxG,\$mxRpk,\$mxProName,\$mxRank,$tmpP,$cpr,$cpn,$i);}
    }
    $tmpProteins[$i][0]=$j-1;

    $oldPep=$peptide;

    if ($noRNAseq) {}else {
    my $readSpt=get_read_spt($peptide,$protein,$paramhash,$pep2reads);
#if ( $outfile eq 'u170k_siRNA_02ug.62439.62439.3.out' ) { print "$outfile,$i,$mxReadSpt,$readSpt,$peptide,$protein\n"; }
    if ($mxlv==4 and $mxReadSpt<$readSpt )
    { $mxReadSpt=$readSpt; $mxRank=$i; $mxProName=$protein; }}
  }
  if ($noRNAseq) { $mxlv=1; $mxRank=0;}
  if ( scalar(@hits)==1 ) { $dCn=1; }
  my $final_rank=$mxRank;
#  for (my $j=1; $j<=$tmpProteins[$final_rank][0]; $j++) {print "$outfile,$final_rank,$j,$tmpProteins[$final_rank][$j]\n";}

  #continue to use old code, but make an interface
  open (IN, $fullpath);
  #seek(IN, 0,0);
  my ($mass, $database) = grep(/\+\smass\s\=/ || /proteins\s\=/, <IN>);
  my $temp=$hits[$final_rank];#interface
  $temp =~ s/\s*\/\s*/\//g; $temp =~ s/(\s)\s+/$1/g; $temp =~ s/^\s+//; $temp =~ s/\s+\Z//;
  $mass =~ s/mass\s+\=\s+([\d\.]+)\s//; $mass = $1;
  $database =~ s/\s+//g; my @temparray = split (/\,/, $database); $database = pop(@temparray);
  if($database=~/hdr/)	{$database=~s/\.hdr$//;	}
  my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $temp);
  $protein=$mxProName;

  if (defined($sorcerer)){    $$paramhash{'sorcererouts'} = 1;  }
  
  if ($dCn < $$paramhash{'min_dCn'} || $XCorr < $$paramhash{'min_XCorr'}){
    $$delhash{'DX'}++;
    if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'DX'}++; } else { $$delhash{'target'}{'DX'}++; }
    return ($database);
  }
  my $peptide = $g; my $red = 0; 

my $ions = $e;
if (!defined($protein)) {die "$outfile,$protein,$mxRank,$noRNAseq,$tmpProteins[0][1]\n";}
  $protein =~ s/\,.*//;
  if (defined($$paramhash{'filter_contaminants'})){
    if ($$paramhash{'filter_contaminants'} == 1){
#      if (defined($$paramhash{'badproteins'}{$protein}) || $protein =~ /keratin/ || $protein=~ /KERATIN/){
     if (defined($$paramhash{'badproteins'}{$protein}) || $protein =~ /keratin/ || $protein=~ /CON_/){ 
       $$delhash{'contaminants'}++;
       if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'contaminants'}++; } else { $$delhash{'target'}{'contaminants'}++; }
        return ($database);
      }
    }
  }
  $proteins{$protein} = 1;
  my $target_decoy_mixed = 0;

  if(defined($h))
  {
    #if(defined($h) && $protein !~ /Random__/ && $protein !~ /Decoy__/  && $protein !~/CON_/){
    $red = $g; $red =~ s/\+//; $peptide = $h;

  # check whether peptide have specific dynamic modifications
  #if (defined($$paramhash{keep_only_mod_peptide}) and $$paramhash{keep_only_mod_peptide})
  if (defined($$paramhash{keep_only_mod_peptide}))
  {
	my $withMod=check_mod_pep($peptide,$$paramhash{keep_only_mod_peptide});
	if (!$withMod) 
	{ 
		$$delhash{no_required_mod}++;
		if ($protein =~ /Decoy/) { $$delhash{decoy}{no_required_mod}++; } else { $$delhash{target}{no_required_mod}++; }
		return ($database); 
	}
  }

    for (my $j=1; $j<scalar(@{$tmpProteins[$final_rank]}); $j++) 
    {
      my $pro=$tmpProteins[$final_rank][$j];#print "$outfile,$final_rank,$j,$tmpProteins[$final_rank][$j]\n";

	#Remove the decoy in the redundancy??
        #mark targets in %proteins
        #if ($pro !~ /Random__/ && $pro !~ /Decoy__/ && $pro !~ /CON_/){
        if ($pro !~ /CON_/ && $pro !~ /Decoy/)
        {
                if (defined($$paramhash{'filter_contaminants'})){
                        if ($$paramhash{'filter_contaminants'} == 1){
                                if (!defined($$paramhash{'badproteins'}{$pro}) && $pro !~ /keratin/ && $pro !~ /KERATIN/ && $pro !~ /CON_/){
                                        $proteins{$pro} = 1;
                                }
                        }
                        else {
                                $proteins{$pro} = 1;
				#print "$outfile,$pro\n";
                        }
                }
        }
        if($pro =~ /Decoy/)
        {
                $target_decoy_mixed = 1;
        }
	last if ($j>=20);
    }

  }

  if($protein !~ /Decoy/ && $target_decoy_mixed == 1)
  {
	$$delhash{'sharedecoy'}++;
  }
  my $testpeptide = $peptide;
  my $mods = ($testpeptide =~ s/([\@\%\&\^\~\$\#\*\^\~]+)/$1/g) || 0;
  $testpeptide =~ s/.\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\../$1/;
  my $mis = ($testpeptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;
  if (defined($$paramhash{'max_peptide_mod'})){
    if ($mods > $$paramhash{'max_peptide_mod'}){
      #$$max_mod_del++;
      
      $$delhash{'maxmod'}++;
      if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'maxmod'}++; } else { $$delhash{'target'}{'maxmod'}++; }
      return ($database);
    }
  }

  if (defined($$paramhash{'max_peptide_mis'})){
    if ($mis > $$paramhash{'max_peptide_mis'}){
      #$$max_mis_del++;
      #print $mis."mismismis\t".$$paramhash{'max_peptide_mis'}."hashmishashmis\t".$max_mis_del."\tccccccccccc\n";
      $$delhash{'maxmis'}++;
      if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'maxmis'}++; } else { $$delhash{'target'}{'maxmis'}++; }
      return ($database);
    }
  }
  my $tryptic = $utils->istryptic($peptide);#returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
  $tryptic -= 1 if ($tryptic>1);
  my $intpep = $peptide;
  $intpep =~ s/[A-Z\-]\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\.[A-Z\-]/$1/;
  if ($$paramhash{'mix_label'} ne 0){  #???
    my $labels = "\[$$paramhash{'mix_label'}\]";
    if ($intpep =~ /$labels[\*\#\@\%\&\^\~\$\^\~]/ && ($intpep =~ /$labels[A-Z]/ || $intpep =~ /$labels\Z/)){
        $$delhash{'mixlabel'}++;
        if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'mixlabel'}++; } else { $$delhash{'target'}{'mixlabel'}++; }
				#print "$outfile $peptide $XCorr $dCn\n";
        return ($database);
    }
  }
  if (defined($$paramhash{'peptide_mod_removal'})){
    if ($$paramhash{'peptide_mod_removal'} ne 0){
      my @labels = split("", $$paramhash{'peptide_mod_removal'});
      for my $aa (@labels){
        my $label = "\[$aa\]";
        if ($intpep =~ /$label[\*\#\@\%\&\^\~\$\^\~]/){
          $$delhash{'modremoval'}++;
          if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'modremoval'}++; } else { $$delhash{'target'}{'modremoval'}++; }
          return ($database);
        }
      }
    }
  }
  if (defined($$paramhash{'peptide_aa_removal'})){
    if ($$paramhash{'peptide_aa_removal'} ne 0){
      my @labels = split("", $$paramhash{'peptide_aa_removal'});
      for my $aa (@labels){
        my $label = "\[$aa\]";
        if ($intpep =~ /$label/){
          $$delhash{'aaremoval'}++;
          if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'aaremoval'}++; } else { $$delhash{'target'}{'aaremoval'}++; }
          return ($database);
        }
      }
    }
  }
  my $nomod = $intpep; $nomod =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;
  if (length($nomod) < $$paramhash{'min_peptide_length'}){
    $$delhash{'length'}++;
    if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'length'}++; } else { $$delhash{'target'}{'length'}++; }
    return ($database);
  }
  my @array = split("", $nomod);
  my $calcmw = $utils->get_MW(\@array); #returns molecular weight of given sequence
  my $num = 0;
	if (defined($$paramhash{$run}{staticmods}{'nterm'})){
    $calcmw += $$paramhash{$run}{staticmods}{'nterm'};
  } elsif (defined($$paramhash{$run}{staticmods}{'cterm'})){
    $calcmw += $$paramhash{$run}{staticmods}{'cterm'};
  }
  for my $aa (keys %{$$paramhash{$run}{staticmods}}){
    $num += ($intpep =~ s/$aa/$aa/g); #print "$outfile: $peptide, $aa, $num\n";
    $calcmw += $num * $$paramhash{$run}{staticmods}->{$aa};
    $num = 0;
  }
  for my $aa (keys %{$$paramhash{$run}{dynamicmods}}){
    my $test = "\[$aa\]";
    $num += ($intpep =~ s/($test[\@\%\&\^\~\$\#\*\^\~])/$1/g);
    $calcmw += $num * $$paramhash{$run}{dynamicmods}->{$aa};
    $num = 0;
  }

#	if ($expmass != $calcmw) { die "$outfile: $expmass, $calcmw, $peptide\n"; }	
  $expmass = $calcmw;#$expmass predefined; why bothers recalculation? why diff?
  if ($$paramhash{'norm_XCorr'}){
    my $length = length($nomod);
    if ($charge == 2){
      if ($length >= 15){
        $XCorr = log($XCorr)/log(15*2);
      } else {
        $XCorr = log($XCorr)/log($length*2);
      }
    } elsif ($charge == 3){
      if ($length >= 25){
        $XCorr = log($XCorr)/log(25*4);
      } else {
        $XCorr = log($XCorr)/log($length*4);
      }
    }
  }

#print "XCorr:",$XCorr,"\n";
#print "Sp:",$Sp,"\n";
#print "rank:",$rank,"\n";
#print "dCn:",$dCn,"\n";
#print "protein:",$protein,"\n";
#print "MH:",$mass,"\n";

###### added by yanji
  # get scan number from outfile name
#  my @outfile_element = split/\./, $outfile;
#  pop(@outfile_element);
#  pop(@outfile_element);
#  my $scan_number = pop(@outfile_element);

  # get original hash, origmsms_hash
#  my $hash_dir = dirname($run)."\/".".hashes";
#  my $hash_ref = retrieve("$hash_dir/origmsms_hash");
#  my %origmsms_hash = %$hash_ref;

  # get rention time, rt, and intensity from the hash table
#  $scan_number = int($scan_number || 0);
#  my $rt = $origmsms_hash{$scan_number}{'rt'};
#  my $intensity = $origmsms_hash{$scan_number}{'intensity'};
###############added by xusheng on 06/14/2012 ###########
#  my $prec_int = $origmsms_hash{$scan_number}{'prec_int'};

  # save rt and intensity
#  $$runhash{$run}{$outfile}{'rt'} = $rt;
#  $$runhash{$run}{$outfile}{'intensity'} = $intensity;
#  $$runhash{$run}{$outfile}{'prec_int'} = $prec_int;

	
	if (defined($$peptidehash{$intpep})){ ##ADDED DMD NOVEMBER 25, 2008 to remove peptide redundancy (decoys) for peptides with real and decoy matches
		my $decoy = 0;
		my $target = 0;
		for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){
		#	print $intpep,"\t",$pros,"\n";
########## changed by xusheng #####################
#Remove all start test ^ sign from the statement
#			$decoy++ if ($pros =~ /Decoy/ || $pros =~ /Decoy__/);
#
			$target ++ if ($pros !~ /Decoy__/ && $pros !~ /Random__/);
			$decoy++ if ($pros =~ /Decoy__/ || $pros =~ /Random__/);
		}

		if ($protein !~ /Decoy__/ && $protein !~ /Random__/ && $decoy > 0){#current out is target; but previous have decoys
			for my $delout (keys %{$$peptidehash{$intpep}{'outfiles'}}){#deltete all previous outfiles of this peptide
				my $run = $$peptidehash{$intpep}{'outfiles'}{$delout}{'run'};
				delete $$runhash{$run}{$delout};
			}
			delete $$peptidehash{$intpep};#delete previous item of this peptide (for renew?)
		} elsif (($protein =~ /Decoy__/ || $protein =~ /Random__/) && $decoy == 0) {#current decoy; deltete current out
			delete $$runhash{$run}{$outfile};
			return;
		}

                for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){#no error??? delete $$peptidehash{$intpep};???
                        if($target>0 and ($pros =~ /Decoy__/ || $pros =~ /Random__/))
                        {
                                $$delhash{'sharedecoy'}++;
				delete $$peptidehash{$intpep}{'proteins'}{$pros};#rm decoy items;keep target; why not delete outfiles?
                        }
		}

	}





 	$$runhash{$run}{$outfile}{'scan'} = $scan;
        $$runhash{$run}{$outfile}{'run'} = $run;
       $$runhash{$run}{$outfile}{'peptide'} = $peptide;
        $$runhash{$run}{$outfile}{'tryptic'} = $tryptic;
        $$runhash{$run}{$outfile}{'mod'} = $mods;
        $$runhash{$run}{$outfile}{'mis'} = $mis;
        $$runhash{$run}{$outfile}{'red'} = $red;
        $$runhash{$run}{$outfile}{'charge'} = $charge;
        $$runhash{$run}{$outfile}{'XCorr'} = $XCorr;
       $$runhash{$run}{$outfile}{'Sp'} = $Sp;#??
       $$runhash{$run}{$outfile}{'rank'} = $rank;
       $$runhash{$run}{$outfile}{'dCn'} = $dCn;
       $$runhash{$run}{$outfile}{'protein'} = $protein;
       $$runhash{$run}{$outfile}{'ions'} = $ions;
       $$runhash{$run}{$outfile}{'MH'} = $mass;
       $$runhash{$run}{$outfile}{'expMH'} = $expmass;
       $$runhash{$run}{$outfile}{'path'} = $fullpath;
       $$runhash{$run}{$outfile}{'intpep'} = $intpep;
       $$runhash{$run}{$outfile}{'protein_priorty_level'} = int($mxlv);
       $$runhash{$run}{$outfile}{'rpk'} = $mxRpk;
       $$runhash{$run}{$outfile}{'psmRank'} = $final_rank;

	if ($noRNAseq) {} else {
	if ( $protein =~ /Decoy/ ) 
	{ 
		my $origPep=orig_pep_for_decoy($paramhash,$peptide);#print "$outfile,$peptide,$origPep\n";
		my @aa=split(//,$origPep); pop(@aa); shift(@aa);
		$origPep=join('',@aa);
		if (defined($$pep2reads{$origPep})) {$$runhash{$run}{$outfile}{'readSpt'} = $$pep2reads{$origPep};}
		else { $$runhash{$run}{$outfile}{'readSpt'} = 0; }
	}
	else 
	{
		if (defined($$pep2reads{$nomod})) {$$runhash{$run}{$outfile}{'readSpt'} = $$pep2reads{$nomod};}
		else { $$runhash{$run}{$outfile}{'readSpt'} = 0; }
	}}


  $$peptidehash{$intpep}{'orig_peptide'} = $peptide;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'scan'} = $scan;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'run'} = $run;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'peptide'} = $peptide;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'tryptic'} = $tryptic;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'mod'} = $mods;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'mis'} = $mis;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'red'} = $red;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'charge'} = $charge;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'XCorr'} = $XCorr;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'Sp'} = $Sp;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'rank'} = $rank;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'dCn'} = $dCn;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'protein'} = $protein;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'ions'} = $ions;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'MH'} = $mass;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'expMH'} = $expmass;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'path'} = $fullpath;
  $$peptidehash{$intpep}{'tryptic'} = $tryptic;
  $$peptidehash{$intpep}{'mis'} = $mis;
  $$peptidehash{$intpep}{'mod'} = $mods;
  $$peptidehash{$intpep}{'readSpt'} = $$runhash{$run}{$outfile}{'readSpt'};

###### added by yanji, add rention time, rt, and intensity
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'rt'} = $rt;
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'intensity'} = $intensity;
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'prec_int'} = $prec_int;
###### end of addition
  my $target = 0;
  for my $key (keys %proteins){
    if($key !~ /Decoy/)
    {
	$target++;
     }
$key =~ s/\#\#//;
    $$peptidehash{$intpep}{'proteins'}{$key}++;
  }

########## count the scan number rather than protein number #####################
  if($target>0)
  {
	my $flag = 0;
       for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){#rm decoy protein names
           if($pros =~ /Decoy/)
           {
		$flag=1;
                #  $$delhash{'sharedecoy'}++;
                  delete $$peptidehash{$intpep}{'proteins'}{$pros};
            }
      }
	if($flag>0)
	{
      		$$delhash{'sharedecoy'}++;
	}
   }
  
  $$peptidehash{$intpep}{'expMH'} = $expmass;
  #%{$$runhash{$run}{$outfile}} = %{$$peptidehash{$intpep}{'outfiles'}{$outfile}};
  #$$runhash{$run}{$outfile}{'intpep'} = $intpep;

  return ($database);
}

sub create_proteinhash
{
	open (LOGFILE, ">>$log_file");
	#select((select(LOGFILE), $|=1)[0]);
	shift @_;
	my ($proteinhash, $peptidehash, $dbhash) = @_;
	my $peptidenum = 1;
	my $pepnum = scalar(keys %$peptidehash);
	while (my ($peptide, $pephash) = each %$peptidehash)
	{
		
		$peptidenum++;
#		printf "\rCreating proteinhash with Peptide_Hash: %d%% done", ($peptidenum/$pepnum)*100;
        	my $SC = scalar(keys %{$$pephash{'outfiles'}});  #print "$SC ";
		my $red = scalar(keys %{$$pephash{'proteins'}});
		
		# get the highest xcorr for each protein 
		my $peptide_max_xcorr=0; my $bestoutfile='';

		my %outfileshash = %{$$pephash{'outfiles'}};
		foreach my $outfile (keys %outfileshash)
		{
			my $xcorr = $outfileshash{$outfile}{'XCorr'};
			if($peptide_max_xcorr<$xcorr)
			{
				$peptide_max_xcorr = $xcorr;
				$bestoutfile = $outfile;
			}
		}

		#unique peptide?
	        if ($red == 1)
		{
			$$pephash{'unique'} = 1;
		} 
		else 
		{
	            $$pephash{'unique'} = 0;
        	}

		for my $protein (keys %{$$pephash{'proteins'}})
		{
			$protein=~s/\#\#//;
			if (!defined($$proteinhash{$protein}))
			{
				$$proteinhash{$protein}{'max_xcorr'} = 0;	
				$$proteinhash{$protein}{'abundance'} = 0;
				$$proteinhash{$protein}{'total'} = 0;
				$$proteinhash{$protein}{'unique'} = 0;
				$$proteinhash{$protein}{'shared'} = 0;
				$$proteinhash{$protein}{'occurrence'} = 0;
				$$proteinhash{$protein}{'sequence'} = $$dbhash{$protein}{'sequence'};
				if ( length($$proteinhash{$protein}{'sequence'})==0 ){die "Zero length protein:$protein!!!\n";}

			if (!defined($$dbhash{$protein}{'annotation'}))
                        {
                                $$dbhash{$protein}{'annotation'} = "not annotated";
                        }
                        $$proteinhash{$protein}{'annotation'} = $$dbhash{$protein}{'annotation'};
			#Cyclin-dependent kinase 18 (Fragment) OS=Homo sapiens GN=CDK18 PE=1 SV=2
			if ($$proteinhash{$protein}{'annotation'} =~ m/\sGN=([\d\w\-]+)/)
			{
				$$proteinhash{$protein}{'GN'}=$1;
			}
			else { $$proteinhash{$protein}{'GN'}='NA'; }
                        if($$proteinhash{$protein}{'annotation'}=~/Uncharacterized/)
                        {
                                $$proteinhash{$protein}{'Unchar'}=1;
                        }
                        else
                        {
                                $$proteinhash{$protein}{'Unchar'}=0;
                        }
                        $$proteinhash{$protein}{'annotation'} =~ s/\[MASS\=(\d+)\]//;
                        if (defined($1))
                        {
                                #$$proteinhash{$protein}{'MW'} = $1;
				my @seqArray = split("", $$dbhash{$protein}{'sequence'});
				$$proteinhash{$protein}{'MW'} = $utils->get_MW(\@seqArray);#returns molecular weight of protein
                        }
                        else
                        {


				if (!defined($$dbhash{$protein}{'sequence'}))
				{
                                      #print Dumper($$dbhash{$protein}) if (!defined($$dbhash{$protein}{'sequence'}));
                                      print "No protein sequence: $protein\n";
                                      print LOGFILE Dumper($$dbhash{$protein}) if (!defined($$dbhash{$protein}{'sequence'}));
                                      print LOGFILE "$protein\n";
                                }
				if(!defined($$dbhash{$protein}{'sequence'}))
				{
					print $protein,"\n";
					exit;
				}

        	                my @seqArray = split("", $$dbhash{$protein}{'sequence'});
                	        $$proteinhash{$protein}{'MW'} = $utils->get_MW(\@seqArray);#returns molecular weight of protein
			}
      			}
 
		       if($$proteinhash{$protein}{'max_xcorr'}<$peptide_max_xcorr)
		       {
		 		$$proteinhash{$protein}{'max_xcorr'} = $peptide_max_xcorr;
		 		$$proteinhash{$protein}{'best_outfile'} = $bestoutfile;
		 		$$proteinhash{$protein}{'best_peptide'} = $peptide;
	       	       }
	  
	      		$$proteinhash{$protein}{'total'}++;
        	 	$$proteinhash{$protein}{'occurrence'} += $SC;
	        	$$proteinhash{$protein}{'abundance'} = ($$proteinhash{$protein}{'occurrence'}*50*1000/($$proteinhash{$protein}{'MW'}+0.0001));
		        $$proteinhash{$protein}{'peptides'}{$peptide}{'tryptic'} = $$pephash{'tryptic'};

			# numbers of unique and shared peptides
		        if ($red == 1)
			{
		        	$$proteinhash{$protein}{'unique'}++;
			        $$proteinhash{$protein}{'peptides'}{$peptide}{'unique'} = 1;
			} 
			else 
			{
			        $$proteinhash{$protein}{'shared'}++;
			        $$proteinhash{$protein}{'peptides'}{$peptide}{'unique'} = 0;
	      		}
    		}
	 }
	# print "\n";
}


sub count_fpr
{
	open (LOGFILE, ">>$log_file");
	#select((select(LOGFILE), $|=1)[0]);

	shift @_;
	my ($proteinhash, $fprhash, $filter, $skipPrintResults) = @_;

	my ($t, $r, $t3, $t2, $t1, $r3, $r2, $r1) = (0,0,0,0,0,0,0,0);
	for my $protein (keys %$proteinhash)
	{
    		if ($protein =~ /Random__/ || $protein =~ /Decoy__/)
		{
      			if ($$proteinhash{$protein}{'total_nomod'} >= 3){ $r3++;$t3++; } 
			elsif ($$proteinhash{$protein}{'total_nomod'}==2){$r2++;$t2++;} 
			else { $r1++;$t1++; }
		      	$r++;
			# added by xusheng for keeping the decoy record ############
			delete ($$proteinhash{$protein}) if ($filter == 0);
    		} 
		else 
		{
      			if ($$proteinhash{$protein}{'total_nomod'} >= 3){$t3++;} 
			elsif ($$proteinhash{$protein}{'total_nomod'}==2){ $t2++;} 
			else { $t1++; }
		}
		$t++;
  	}
	$$fprhash{'t'} = $t; $$fprhash{'r'} = $r;
	if ( !defined($skipPrintResults) || $skipPrintResults==0 ) {
	printf "  Total proteins:  %d targets and %d decoys (FDR = %.2f%%)\n", $t-$r,$r, ($r/($t+0.01))*100;
	printf LOGFILE "  Total proteins:  %d targets and %d decoys (FDR = %.2f%%)\n", $t-$r,$r, ($r/($t+0.01))*100;
	}

	$$fprhash{'t1'} = $t1; $$fprhash{'r1'} = $r1;
	$$fprhash{'t2'} = $t2; $$fprhash{'r2'} = $r2;
	$$fprhash{'t3'} = $t3; $$fprhash{'r3'} = $r3;
}

sub create_sumpeptidehash{
	shift @_;
	my ($sumhash, $peptidehash) = @_;
	
	for my $group (keys %$sumhash){
		my %pephash = %{$$sumhash{$group}{'peptide'}};
		if (scalar(keys %$peptidehash) == 0){
			%$peptidehash = %pephash;
		} else {
			while (my ($peptide, $hash) = each %pephash){
				if (!defined($$peptidehash{$peptide})){
					%{$$peptidehash{$peptide}} = %$hash;
				} else {
					while (my ($outfile, $outhash) = each %{$$hash{'outfiles'}}){
						%{$$peptidehash{$peptide}{'outfiles'}{$outfile}} = %$outhash;
					}
				}
			}
		}
	}
}

sub peptide_Qvalue
{
	shift @_;
	my ($peptidehash,$featurehash) = @_;
	my (%smallhash,$lda); 

	$lda=0;
	while (my ($intpep, $intpephash) = each %$peptidehash)
	{
		$smallhash{$intpep}{random}=$$peptidehash{$intpep}{random};
		$smallhash{$intpep}{XCorr}=0; 
		while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}})
		{
			my $score;
			if (defined($$featurehash{$outfile}{lda_p})) { $score=1-$$featurehash{$outfile}{lda_p}; $lda=1; }
			if ( $smallhash{$intpep}{XCorr}<$score ) { $smallhash{$intpep}{XCorr}=$score; }
		}
	}

	if ($lda) { _XCorr_qValue(\%smallhash); }

	while (my ($intpep, $intpephash) = each %$peptidehash)
	{
		if (defined($smallhash{$intpep}{qvalue}))
		{ $$peptidehash{$intpep}{Qvalue}=$smallhash{$intpep}{qvalue}*100; }
		else { die "No peptide q-value for $intpep (peptide_Qvalue)\n"; }
	}
}

sub _XCorr_qValue
{
	my ($smallhash)=@_;

	#sort by XCorr (ascending): from low to high
        my @array;
        for my $outhash (sort {$$a{XCorr} <=> $$b{XCorr}} values %$smallhash)
        {
                push(@array, $outhash);
        }

        #count decoy numbers: from high score to low
        my @dCount; $dCount[scalar(@array)]=0;
        for (my $i=scalar(@array)-1; $i>=0; $i--)
        {
                my $tmp=$array[$i];
                $dCount[$i]=$dCount[$i+1];
                if ( $$tmp{random} == 1 ) { $dCount[$i]++; }
        }

	#open(OUT,">>_XCorr_qValue");
        #compute q values
        for (my $i=0; $i<scalar(@array); $i++ )
        {
                my $tmp=$array[$i];
                #$$tmp{qvalue}=2*$dCount[$i]/(scalar(@array)-$i);
                if ( (scalar(@array)-$i-$dCount[$i]) )
		{
	                $$tmp{qvalue}=$dCount[$i]/(scalar(@array)-$i-$dCount[$i]);
		}
		else { $$tmp{qvalue}=1; }
		#print OUT "$$tmp{XCorr},$$tmp{random},$$tmp{qvalue}\n"; 
        }
	#print OUT "\n";
	#close OUT;

	# correct/streamline scan p-value:
	# 1) sort search score: low => high
	# 2) q-value: high => low
	my $t=$array[0];
	my $highestQ=$$t{qvalue};
	for (my $i=1; $i<scalar(@array); $i++ )
	{
		my $tmp=$array[$i];
		if ( $$tmp{qvalue}>$highestQ) 
		{ 
			if ($$tmp{random} == 0) # only for targets; leave decoys untouched
			{ $$tmp{qvalue}=$highestQ;}  # streamline q value
		}
		else { $highestQ=$$tmp{qvalue}; }
	}
}

sub buildFeatureHash
{
	shift @_;
        my ($peptidehash,$featurehash)=@_;

        while (my ($intpep, $intpephash) = each %$peptidehash)
        {
                $$intpephash{'random'}=1;
                while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}})
                {
                        $$featurehash{$outfile}{'XCorr'}=$$outhash{'XCorr'};
                        $$featurehash{$outfile}{'dCn'}=$$outhash{'dCn'};
                        $$featurehash{$outfile}{'ppm'}=$$outhash{'ppm'};
                        #$$featurehash{$outfile}{'ppm_passed'}=$$outhash{'ppm'};
                        $$featurehash{$outfile}{'charge'}=$$outhash{'charge'};
                        $$featurehash{$outfile}{'mis'}=$$outhash{'mis'};
                        $$featurehash{$outfile}{'mod'}=$$outhash{'mod'};
                        $$featurehash{$outfile}{'tryptic'}=$$outhash{'tryptic'};

                        my $nomod=$intpep;
                        $nomod =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;

                        $$featurehash{$outfile}{'peptide'}=$intpep;
                        $$featurehash{$outfile}{'peptideLength'}=length($nomod);

                        if ($$outhash{'protein'} =~ /Random__/ || $$outhash{'protein'} =~ /Decoy__/)
                        {
                                $$featurehash{$outfile}{'random'} = 1;
                        }
                        else
                        {
                                $$featurehash{$outfile}{'random'} = 0; $$intpephash{'random'}=0;
                        }

                }
        }
}

sub printFeatureHash
{
	shift @_;
	my ($output,$featurehash)=@_;

	open(FTR,">$output");
        print FTR "outfile\txcorr\tdCn\tpeptideLength\ttryptic\tppm\tmod\tmc\tcharge1\tcharge2\tcharge3\tcharge4\ttargetDecoy\n";
        foreach my $outfile (keys %{$featurehash})
        {
                 print FTR "$outfile\t$$featurehash{$outfile}{XCorr}\t$$featurehash{$outfile}{dCn}\t$$featurehash{$outfile}{peptideLength}\t$$featurehash{$outfile}{tryptic}\t$$featurehash{$outfile}{ppm}\t$$featurehash{$outfile}{mod}\t$$featurehash{$outfile}{mis}\t";
                 for (my $i=1; $i<=4; $i++)
                 {
                         if ( $$featurehash{$outfile}{charge} == $i ) { print FTR "1\t"; }
                         else { print FTR "0\t"; }
                 }
                 if ($$featurehash{$outfile}{random}==1) { print FTR "decoy"; }
                 else { print FTR "target"; }
                 print FTR "\n";
         }
         close FTR;

}

sub attchLDAresult2FeatureHash
{
	shift @_;
	my ($inputFile,$featurehash)=@_;

	my $inputDirectory = dirname($inputFile);
	#my $input=$inputDirectory."\/".(split(/\./, $inputFile))[0]."_LDA.txt";
	my $input=basename($inputFile); $input =~ s/\.[^.]+$//; # print "$input\n";
	my $input=$inputDirectory."\/".$input."_LDA.txt"; #print "$input\n";
	#my $input.="_LDA.txt";

	open(IN,$input) || die "Cannot open LDA result file $input!!!\n";
	my $line=<IN>;
	while (<IN>)
	{
		chomp;
		my @t=split(/\t/,$_);
		if (defined($$featurehash{$t[0]}))
		{
			$$featurehash{$t[0]}{lda_p}=$t[13];
		}
		#else {die "Not ";}
	}
	close IN;
}

sub rLDA_analysis
{
	shift @_;
	my ($inputFile,$search_engine,$R_script)=@_;

	my $nTotLines = `wc -l $inputFile`;
	chomp($nTotLines);
	$nTotLines = (split(/ /, $nTotLines))[0] - 1;
	my $inputDirectory = dirname($inputFile);
	my $inputFileBaseName = basename($inputFile);
	if ($inputDirectory eq "\.") 
	{
        	$inputDirectory = cwd();
	}

	## Running LDA in R
	my $R = Statistics::R -> new();
	$R -> set('inputFile', $inputFileBaseName);
	$R -> set('inputDirectory', $inputDirectory);
	$R -> set('nTotLines', $nTotLines);
	$R -> set('search_engine', $search_engine);
	$R -> run_from_file($R_script);
=head
	$R -> run(q`setwd(inputDirectory);
	#suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
        library(MASS);
        library(qvalue);
	#search_engine=\$search_engine;
        ## Data loading
        #table5rows = read.table(inputFile, sep = "\t", header = T, nrows = 5);
        table5rows = read.table(inputFile, sep = "\t", header = T, nrows = 30);
        classes = sapply(table5rows, class);
        data = read.table(inputFile, sep = "\t", header = T, nrows = nTotLines, colClasses = classes);
        ## Initialization
        numSamples = dim(data)[1];
        pval = rep(1, numSamples);
        qval = rep(1, numSamples);
        fdr = rep(1, numSamples);
        ## Independent LDA model for each charge state
        chargeStates = c(1, 2, 3, 4);
        for (charge in chargeStates) {
                colName = paste("charge", charge, sep = "");
                subInd = which(data[colnames(data) == colName] == 1);
                if (length(subInd) > 30) {
                        subData = data[subInd, ];
                        ## Generation of a data matrix as in Du et al.
                        ## Three variables: log(XCorr)/log(peptideLength), sqrt(dCn) and ppm
                        if (!is.na(subData$ppm))
			{
                        if ( search_engine == 'jump' )
			{
				if (sd(subData$dCn)>0)
				{
	                        	X = as.data.frame(cbind(subData$xcorr, sqrt(subData$dCn), abs(subData$ppm))); # maybe for Jscore?
				}
				else
				{
					X = as.data.frame(cbind(subData$xcorr,abs(subData$ppm))); # dJn is constant
				}
                        }
			else
			{
	                        X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn), abs(subData$ppm)));
                        	#X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn)));
			}
			} # end of if (sd(subData$ppm)>0)
			else # ppm not considered
			{
                        if ( search_engine == 'jump' )
			{
				if (sd(subData$dCn)>0)
				{
	                        	X = as.data.frame(cbind(subData$xcorr, sqrt(subData$dCn))); # for Jscore
				}
				else
				{
					X = as.data.frame(cbind(subData$xcorr)); # dJn is constant
				}
                        }
			else
			{
	                        #X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn), abs(subData$ppm)));
                        	X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn)));
			}
			} # end of else (sd(subData$ppm)>0)
                        group = subData$targetDecoy;
                        ## LDA modeling
                        if (length(table(group))>=2 & table(group)[1]>1 & table(group)[2]>1) {
				ldaModel = lda(group ~ ., X);
				## Calculation of LDA scores of the samples
				scores = as.numeric(predict(ldaModel, X)$x);
				## Let the scores of "target"s be always greater than "decoy"s
				if (mean(scores[group == "target"]) < mean(scores[group == "decoy"])) {
					scores = -scores;
				}
				## Null hypothesis distribution by fitting "decoy" scores to a normal distribution
				H0 = fitdistr(scores[group == "decoy"], "normal");
				p = 1 - pnorm(scores, mean = H0$estimate[1], sd = H0$estimate[2]);
				pval[subInd] = p;
				## Calculate q-values and FDR (Benjamini-Hochberg)
				qval[subInd] = qvalue(p)$qvalue;
				fdr[subInd] = p.adjust(p, method = "BH");
			} else {
				pval[subInd] = 1;
				qval[subInd] = 1;
				fdr[subInd] = 1;
			}
                }
		else
		{
			pval[subInd] = 1;
			qval[subInd] = 1;
			fdr[subInd] = 1;
		}
        }
        write.table(cbind(pval, qval, fdr), file = "LDA_Result.txt", sep = "\t");`);
=cut
	$R -> stop();
	print "LDA modeling in R has been finished\n";

	## Obtain p-values, q-values and FDR values
	my $LDAfile = $inputDirectory."\/"."LDA_Result.txt";
	open (LDAFILE, $LDAfile) or die "Cannot open $LDAfile\n";
	my (@pvalue, @qvalue, @fdr) = ();
	my $k = 0;
	while (<LDAFILE>) {
        	chomp;
	        if ($k > 0) {
        	        my @elements = split(/\t/, $_);
	                push (@pvalue, $elements[1]);
        	        push (@qvalue, $elements[2]);
                	push (@fdr, $elements[3]);
	        }
        	print "Gathering the results ($k / $nTotLines lines)\r";
	        $k++;
	}
	close (LDAFILE);
	system ("rm $LDAfile");
	print "\n";

	## Output file generation
	my $outputDirectory = $inputDirectory;
	#my $outputFile = $inputDirectory."\/".(split(/\./, $inputFile))[0]."_LDA.txt";
	#my $outputFile = (split(/\./, $inputFile))[0]."_LDA.txt";
	my $outputFile = basename($inputFile);
	$outputFile =~ s/\.txt/_LDA\.txt/;
	$outputFile = dirname($inputFile) . "/" . $outputFile;
	open (OUTFILE, ">", $outputFile) or die "Cannot generate $outputFile\n";
	open (INFILE, "<", $inputFile) or die "Cannot open $inputFile\n";
	my $nLine = -1;
	while (<INFILE>) {
        	## Writing to an output file
	        chomp;
        	if ($nLine == -1) {
	                my $header = $_."\tPvalue\tQvalue\tFDR\n";
        	        print OUTFILE $header;
	        } else {
        	        my $line = $_."\t$pvalue[$nLine]\t$qvalue[$nLine]\t$fdr[$nLine]\n";
	                print OUTFILE $line;
        	}
	        $nLine++;
	}
	close (INFILE);
	close (OUTFILE);

	return($outputFile);
}

sub find_XCorrdCn_cutoffs
{
  	open (LOGFILE, ">>$log_file");
	#fsync(\*LOGFILE) or die "fsync: $!"; sync();
	select((select(LOGFILE), $|=1)[0]);

	shift @_;
	my $pepfpr;
	my $peptide_fpr;
	my ($peptidehash, $paramhash, $fpr, $featurehash, $del, $printResults) = @_;
	my ($overall_target, $overall_decoy) = (0,0);

	if (!defined($$paramhash{min_decoy_num_for_XCorr_filter}) && !defined($$paramhash{min_outfile_num_for_XCorr_filter})) { die "Missing min_decoy_num_for_XCorr_filter/min_outfile_num_for_XCorr_filter!!! Please update params file!!!\n"; }

	my %combohash = %{$$paramhash{'12combinations'}};

	#print "peptidehash in find_XCorrdCn_cutoffs:",scalar(keys %$peptidehash),"\n";

	my %bighash; my $dCnGroupNum=100;
	while (my ($intpep, $intpephash) = each %$peptidehash)
	{
		$$intpephash{'random'}=1;
		while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}})
		{
			#next if ($$outhash{status}<1);
			# 12 combinations filtering
			my $t;
			if ($$outhash{'tryptic'} == 0 ) { $t='NT'; }
			elsif ($$outhash{'tryptic'} == 1 ) { $t='PT'; }
			elsif ($$outhash{'tryptic'} == 2 ) { $t='FT'; }
			if ( !defined($combohash{"$t$$outhash{'charge'}"}) ) { die "Not define 12 combinations item: $t$$outhash{'charge'}!!!\n"; }
			if ( $combohash{"$t$$outhash{'charge'}"}==0 ) { delete ($$intpephash{'outfiles'}{$outfile}); next; }

			# $del option for deletion
			if ($del == 1)
			{
				if (!defined($$intpephash{'outfiles'}{$outfile}{'filtered'})) { delete($$intpephash{'outfiles'}{$outfile}); next;}
			}

			# priority score calculation: the larger, the less confident
			my $ps='';

#=head
			my $nomod=$intpep;
			if ( defined($$paramhash{group_by_pepLength}) and $$paramhash{group_by_pepLength} )
			{
				$nomod =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;;
				my $pepL=length($nomod);
				if ( $pepL<10 ) { $ps.="0$pepL"; }
				else { $ps.="$pepL"; }
			}
#=cut
			#if ( length($intpep) < $$paramhash{'long_peptide_length'} ) { $ps.='1'; }
			#else { $ps.='0'; }

			if ( defined($$paramhash{group_by_tryp}) and $$paramhash{group_by_tryp} )
			{ $ps.=2-$$outhash{'tryptic'};}
			if ( defined($$paramhash{group_by_mod}) and $$paramhash{group_by_mod} )
			{ $ps.=$$outhash{'mod'}; }
			if ( defined($$paramhash{group_by_mis}) and $$paramhash{group_by_mis} )
			{ $ps.=$$outhash{'mis'}; }
			if ( defined($$paramhash{group_by_charge}) and $$paramhash{group_by_charge} )
			{
				if ( $$outhash{'charge'}==0 ) { $ps.=9; } 
				elsif ( $$outhash{'charge'}==1 ) { $ps.=8; }
				else { $ps.=$$outhash{'charge'}; }
			}
=head
			if ( length($intpep) < $$paramhash{'long_peptide_length'} ) { $ps.='1'; }
			else { $ps.='0'; }
=cut
			if ( defined($$paramhash{group_by_dCn}) and $$paramhash{group_by_dCn} )
			{
				my $dCn=$$outhash{'dCn'}; if ($dCn>=1) { $dCn=0.99; }
				my $dCnIndex=$dCnGroupNum-int($dCn*$dCnGroupNum)-1;
				if ($dCnIndex<10) { $ps.="0$dCnIndex"; }
				else { $ps.=$dCnIndex; }
			}
=head
			if (defined($$outhash{'tagLength'}))
			{
				if ($$outhash{'tagLength'}>=4) { $ps.=0; }
				else { $ps.=9-$$outhash{'tagLength'}; }
			}
=cut
			# bighash construction
			if (defined($$paramhash{'FDR_filtering_method'}) and $$paramhash{'FDR_filtering_method'} eq 'LDA') { $bighash{$outfile}{'XCorr'}=1-$$featurehash{$outfile}{lda_p}; }
			else  { $bighash{$outfile}{'XCorr'}=$$outhash{'XCorr'};}
			#$bighash{$outfile}{'XCorr'}=1-$$featurehash{$outfile}{lda_p};
			$bighash{$outfile}{'dCn'}=$$outhash{'dCn'};
			$bighash{$outfile}{'run'}=$$outhash{'run'};
			$bighash{$outfile}{'ps'}=$ps;
			$bighash{$outfile}{'peptide'}=$intpep;
			if ($$outhash{'protein'} =~ /Random__/ || $$outhash{'protein'} =~ /Decoy__/)
			{
				$bighash{$outfile}{'random'} = 1;
			}
			else
			{
				$bighash{$outfile}{'random'} = 0; $$intpephash{'random'}=0;
			}
		}
	}

	#print "big hash in find_XCorrdCn_cutoffs:",scalar(keys %bighash),"\n";

	my $decoyNum=0; my %smallhash; #my ($overall_target,$overall_decoy)=(0,0);
	my $psmK=0;
	my %TDhash; my $outfileCount=0;
	open(OUT,">>\.idsumTmp\/Grouping_for_XCorr_filtering.txt");
	printf OUT "\nMinimum Peptide FPR (might take minutes): %.2f\n",$fpr;
	print OUT "\ncount\tps\toutfile\tpeptide\tXCorr\tdCn\trandom\tqvalue\n";

if (defined($$paramhash{'FDR_filtering_method'}) and $$paramhash{'FDR_filtering_method'} eq 'LDA')
{
	my %smallhash=%bighash;
	_XCorr_qValue(\%smallhash);
	my $psm_accepted=1;
	for my $out ( sort {$smallhash{$b}{XCorr} <=> $smallhash{$a}{XCorr}} keys %smallhash )
	{
=head
		if ( !defined($$peptidehash{$smallhash{$out}{peptide}}{Qvalue}) || 
		$$peptidehash{$smallhash{$out}{peptide}}{Qvalue}>$smallhash{$out}{qvalue} )
		{
			$$peptidehash{$smallhash{$out}{peptide}}{Qvalue}=$smallhash{$out}{qvalue};
		}
=cut
		print OUT "$outfileCount\t$smallhash{$out}{ps}\t$out\t$smallhash{$out}{peptide}\t$smallhash{$out}{XCorr}\t$smallhash{$out}{dCn}\t$smallhash{$out}{random}\t$smallhash{$out}{qvalue}\t";

		my $run=$smallhash{$out}{run};
		if (!defined($TDhash{$run}))  { $TDhash{$run}{rm}{D}=$TDhash{$run}{rm}{T}=$TDhash{$run}{ps}{D}=$TDhash{$run}{ps}{T}=0; }

		if (  $psm_accepted==1 and $smallhash{$out}{qvalue}<=$fpr/100 )
		{
			print OUT "accepted\n";
			if ( $smallhash{$out}{random}==1 ) { $overall_decoy++; $TDhash{$run}{ps}{D}++;  } else { $overall_target++; $TDhash{$run}{ps}{T}++; }
		}
		else
		{
			$psm_accepted=0;
			print OUT "filtered\n";
			if ( $smallhash{$out}{random}==1 ) { $TDhash{$run}{rm}{D}++; } else  { $TDhash{$run}{rm}{T}++; }
			delete ($$peptidehash{$smallhash{$out}{peptide}}{outfiles}{$out});
		}
	}
}
else
{
	for my $outfile ( sort { $bighash{$b}{ps} <=> $bighash{$a}{ps} } keys %bighash )
	{
		foreach my $k (keys %{$bighash{$outfile}}) { $smallhash{$outfile}{$k}=$bighash{$outfile}{$k}; }
		#print OUT "$outfileCount\t$smallhash{$outfile}{ps}\t$outfile\t$smallhash{$outfile}{XCorr}\t$smallhash{$outfile}{random}\t$smallhash{$outfile}{qvalue}\n";

		if ( $bighash{$outfile}{'random'} == 1 ) { $decoyNum++; }
		$psmK++;# counting outfile number for bin
		$outfileCount++; # counting total outfile number

		#if ( $decoyNum>=$$paramhash{'min_decoy_num_for_XCorr_filter'} ) # controlled by min decoy number
		#if ( $psmK>=$$paramhash{'min_outfile_num_for_XCorr_filter'} ) # controlled by min outfile number
		if ( $psmK>=$$paramhash{'min_outfile_num_for_XCorr_filter'} and scalar(keys %bighash)-$outfileCount>=$$paramhash{'min_outfile_num_for_XCorr_filter'} or scalar(keys %bighash)==$outfileCount ) # controlled by min outfile number
		{
			print OUT "\n";
			_XCorr_qValue(\%smallhash);#print "";
			# use qvalue to clean peptidehash
			for my $out ( sort {$smallhash{$b}{qvalue} <=> $smallhash{$a}{qvalue}} keys %smallhash )
			#my $psm_accepted=1;
			#for my $out ( sort {$smallhash{$b}{XCorr} <=> $smallhash{$a}{XCorr}} keys %smallhash )
			{
				print OUT "$outfileCount\t$smallhash{$out}{ps}\t$out\t$smallhash{$out}{peptide}\t$smallhash{$out}{XCorr}\t$smallhash{$out}{dCn}\t$smallhash{$out}{random}\t$smallhash{$out}{qvalue}\t";

				my $run=$smallhash{$out}{run};
				if (!defined($TDhash{$run}))
				{ $TDhash{$run}{rm}{D}=$TDhash{$run}{rm}{T}=$TDhash{$run}{ps}{D}=$TDhash{$run}{ps}{T}=0; }

				#if ( $psm_accepted==1 and $smallhash{$out}{qvalue}<=$fpr/100 )
				if ( $smallhash{$out}{qvalue}<=$fpr/100 )
				{
					print OUT "accepted\n";
					if ( $smallhash{$out}{random}==1 ) { $overall_decoy++; $TDhash{$run}{ps}{D}++;  }
					else { $overall_target++; $TDhash{$run}{ps}{T}++; }
				}
				else
				{
					#$psm_accepted=0;
					print OUT "filtered\n";
					if ( $smallhash{$out}{random}==1 ) { $TDhash{$run}{rm}{D}++; }
					else  { $TDhash{$run}{rm}{T}++; }
					delete ($$peptidehash{$smallhash{$out}{peptide}}{outfiles}{$out});
				}
			}

			undef(%smallhash); $decoyNum=0; $psmK=0;
		}
	}
}
	close OUT;

	if (defined($$paramhash{printRun_infor}) and $$paramhash{printRun_infor}==1)
	{
        foreach my $run (keys %TDhash)
        {
                print "For run $run:\n  $TDhash{$run}{ps}{T} targets and $TDhash{$run}{ps}{D} decoys passed\n";
                print "  $TDhash{$run}{rm}{T} targets and $TDhash{$run}{rm}{D} decoys deleted\n";
                print LOGFILE "For run $run:\n  $TDhash{$run}{ps}{T} targets and $TDhash{$run}{ps}{D} decoys passed\n";
                print LOGFILE "  $TDhash{$run}{rm}{T} targets and $TDhash{$run}{rm}{D} decoys deleted\n";
        }
	}

	#print "     Deleted $num outfiles that shared the same scans\n";
	if ($printResults) {
	printf "Accepting %d outfiles: %d targets and $overall_decoy decoys (FDR = %.2f%%)\n", $overall_target+$overall_decoy,$overall_target,100*$overall_decoy/($overall_target+0.01);
	printf LOGFILE "Accepting %d outfiles: %d targets and $overall_decoy decoys (FDR = %.2f%%)\n", $overall_target+$overall_decoy,$overall_target,100*$overall_decoy/($overall_target+0.01);
	}

	# if there are multiple outfiles for one scan (due to un-predicted charge state), pick the best one
	# establish scanhash
	my %scanhash;
	while (my ($peptide, $pephash) = each %$peptidehash){
	    	while (my ($outfile, $outhash) = each %{$$pephash{'outfiles'}}){
      			my ($run, $scan) = ($$outhash{'run'}, $$outhash{'scan'});
			next if (!defined($outfile)||!defined($run)||!defined($scan));
		      	$scanhash{$run}{$scan}{$outfile}{'peptide'} = $peptide;
		      	$scanhash{$run}{$scan}{$outfile}{'charge'} = $$outhash{'charge'};
		      	$scanhash{$run}{$scan}{$outfile}{'XCorr'} = $$outhash{'XCorr'};
    		}
  	}
	#if there are multiple outfiles (e.g., with different charges) for a scan, pick the best one
  	my $num = 0;
	for my $runhash (values %scanhash){# %runhash{$scan}{$outfile}
  		for my $hash (values %$runhash){# %hash{$outfile}
		    	if (scalar(keys %$hash) > 1)# multiple outfiles for one scan: with different charges
			{
				my $bestXCorr = 0;
				my $bestoutfile;
			      	for my $outfile (keys %$hash)
				{
					my $orig_XCorr = $$hash{$outfile}{'XCorr'};
					if ($$hash{$outfile}{'charge'} == 3)
					{
						if ($bestXCorr<$orig_XCorr-1)
						{
							$bestoutfile = $outfile;
							$bestXCorr = $orig_XCorr-1;
						}
					} 
					else 
					{
						if ($bestXCorr <= $orig_XCorr)
						{
							$bestoutfile = $outfile;
							$bestXCorr = $orig_XCorr;
						}
					}	
      				}
				for my $outfile (keys %$hash){
					next if ($outfile eq $bestoutfile);
					my $peptide = $$hash{$outfile}{'peptide'};
					$num++;
					delete $$peptidehash{$peptide}{'outfiles'}{$outfile};
				}
    			}#end of if
  		}
	}

	# calculate the peptide FPR 
	# %temp_pep{$peptide}=$proteinName
	my %temp_pep;
        while (my ($intpep, $intpephash) = each %$peptidehash){
                while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}}){
			$temp_pep{$intpep}=$$outhash{'protein'};
		}
	}
	my $Decoy_num=0;
	
	foreach my $peptide (keys %temp_pep)
	{
		next if (!$temp_pep{$peptide});
		if($temp_pep{$peptide} =~ /Decoy/)
		{
			$Decoy_num++;
		}
	}
	
	$pepfpr = $Decoy_num/((scalar keys %temp_pep)-$Decoy_num+0.01)*100;# fpr = decoy#/target#

	my $target_peptide_num = (scalar keys %temp_pep)-$Decoy_num;
	my $adjust_pepfpr = $Decoy_num/(0.01+$target_peptide_num)*100;# fpr = decoy#/target#

	#if no outfiles left, delete that petide
	while (my ($intpep, $intpephash) = each %$peptidehash){
		#while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}}){ #ADDED next 5 lines DMD 10/16/2009
		#	 if ($$outhash{'protein'} =~ /Decoy__/ || $$outhash{'protein'} =~ /Decoy/){
		#			delete $$intpephash{'outfiles'}{$outfile};
		#		}
		#}
		if (scalar(keys %{$$intpephash{'outfiles'}}) == 0){
			delete ($$peptidehash{$intpep});
		}
		
	}

	# unique scan: checi if there are outfiles with  2 or more charges (multi unique scans) for one peptide
	my $uniqScanD=0; my $uniqScanT=0;
	while (my ($intpep, $intpephash) = each %$peptidehash)
	{
		my %chargeHash;
		while (my ($outfile, $outhash) = each %{$$intpephash{'outfiles'}})
		{
			$chargeHash{$$outhash{charge}}='';
		}
		$$intpephash{'uniq_scan_num'}=scalar(keys %chargeHash);
		if ( $$intpephash{random}==1 ) { $uniqScanD+=$$intpephash{'uniq_scan_num'}; } 
		else { $uniqScanT+=$$intpephash{'uniq_scan_num'}; }
	}

	if ($printResults) {
	printf "  Unique precursor ions: %d targets and %d decoys (FDR = %.2f%%)\n", $uniqScanT,$uniqScanD,$uniqScanD*100/($uniqScanT+0.01);
	printf LOGFILE "  Unique precursor ions: %d targets and %d decoys (FDR = %.2f%%)\n", $uniqScanT,$uniqScanD,$uniqScanD*100/($uniqScanT+0.01);
	}

	close LOGFILE;

 	return ($pepfpr,$Decoy_num,$target_peptide_num,$adjust_pepfpr);
}

sub ppm_hard_cut
{
	my ($runhash,$ppmCut)=@_;

	print "Not enough good scans: use $ppmCut ppm to filter scans without mass calibration\n";
	while (my ($outfile, $outhash) = each %$runhash)
	{
		if (defined($$outhash{ppm}))
		{
			if (abs($$outhash{ppm})<=$ppmCut) { $$outhash{status}=1; }
			else { $$outhash{status}=-1; }
		}
		else {die "Not defined ppm: $outfile!!!\n";}
	}
}

sub acumass
{
	open (LOGFILE, ">>$log_file");
	my $deltaC=1.0033548378;
	#fsync(\*LOGFILE) or die "fsync: $!";	sync();
	#select((select(LOGFILE), $|=1)[0]);

	shift @_;
	my ($run, $runhash, $paramhash) = @_;
	my $masstype = $$paramhash{'mass_consideration'}{'value'};
	my $bf = $$paramhash{'scan_events'};
	#print "\nRunning Mass Accuracy Filtering for $run\n";
        #print LOGFILE "\nRunning Mass Accuracy Filtering for $run\n";

	if (!defined($$paramhash{'XCorr'}))
	{
		if ( $$paramhash{'search_engine'} eq 'jump' ) { $$paramhash{'XCorr'}=20; }
		else { $$paramhash{'XCorr'}=2.5; }
	}

	# Create a hash of only good outfiles using thresholds set in parameter file
	my ($XCorr, $dCn, $tryptic, $charge) = ($$paramhash{'XCorr'}, $$paramhash{'dCn'}, $$paramhash{'trypticity'}, $$paramhash{'charge'});
	$XCorr = .1 if ($$paramhash{'norm_XCorr'});
	my %goodhash;
	my $largeDiff_k=0;
	my $pass_k=0;
	my $rest_k=0;
	#open(OUT,">acumass_debug_outfiles");
	#print OUT "out\txcorr\tdcn\ttryptic\tcharge\tpass\n";
	while (my ($outfile, $outhash) = each %$runhash){
		my ($mass, $expmass, $charge) = ($$outhash{'MH'}, $$outhash{'expMH'}, $$outhash{'charge'});
		# correct mass for jump (how many units of H?)
		if ( $$paramhash{'search_engine'} eq 'jump' )
		{
			my $dM=$mass - $expmass;
			#my $unit=int($dM);
			my $unit=floor($dM+0.5);
			#if (abs($unit) >= 1) { print "$outfile, $dM=$mass - $expmass,$unit\n"; }
			$mass=$mass-$unit*$deltaC;
		}

		# mass diff
		my $diff = $mass - $expmass;
		my ($status,$Mtype);
		if ( abs($mass - $expmass)<0.4) { $status = $Mtype=0; }
		else {  $status=$Mtype = -1;}
		
		#if |mass - expmass|<0.4 => diff = mass - expmass; status = mtype=0; 
		#else diff = mass - expmass; status=mtype = -1;
		#my ($diff, $status, $Mtype) = get_realdiff($mass, $expmass, $masstype);
		#print "$outfile;$mass;$expmass;$charge;$diff;$status;$Mtype\n";
		#my $ppm = ($diff/$mass/$charge)*(10**6);
		my $ppm = ($diff/$expmass)*(10**6); # modified on 10/10/14
    $$outhash{ppm} = $ppm; $$outhash{status} = $status; $$outhash{Mtype} = $Mtype;
		#next if ($status == -1);
		if ($status == -1) { $largeDiff_k++; next; }
		#print OUT "$outfile\t$$outhash{'XCorr'}\t$$outhash{'dCn'}\t$$outhash{'tryptic'}\t$$outhash{'charge'}\t";
		#if ($$outhash{'XCorr'} >= $XCorr && $$outhash{'dCn'} >= $dCn && $$outhash{'tryptic'} == $tryptic && $$outhash{'charge'} == $charge)
		if ($$outhash{'XCorr'} >= $XCorr && $$outhash{'dCn'} >= $dCn && $$outhash{'tryptic'} == $tryptic)
		{
			$pass_k++;#print OUT "1\n";

####### changed by xusheng ###########
#
#			next if ($$outhash{'protein'} =~ /Decoy__/ || $$outhash{'protein'} =~ /Decoy/);
			next if ($$outhash{'protein'} =~ /Random__/ || $$outhash{'protein'} =~ /Decoy__/);			
			%{$goodhash{$outfile}} = %$outhash; $$outhash{status} = 1;
		}
#		else { print OUT "0\n"; }
		$rest_k++;
	}
#	print "largeDiff_k: $largeDiff_k; pass_k: $pass_k; rest_k: $rest_k\n";
#	close OUT;
	#exit;

  # Do preliminary cyclic 2SD filter of good outfiles hash
  my ($sum, $sum2, $n) = (0,0,scalar(keys %goodhash));
  while (my ($outfile, $hash) = each %goodhash){
    my $ppm = $$hash{'ppm'};
    $sum += $ppm; $sum2 += $ppm**2;
  }

  # if not enoughg good scans
  if (($n < 1) || ($n <= 3*$bf))
  {
    if (!defined($$paramhash{static_cutoff_without_mass_calib}))
    {
	while (my ($outfile, $outhash) = each %$runhash){
		$$outhash{status} = 1;
	}
	print "Skipping accurate mass filtering because of not enough good scans.\n";
	print LOGFILE "Skipping accurate mass filtering because of not enough good scans.\n";
	return;
    }
    else
    {
      ppm_hard_cut($runhash,$$paramhash{static_cutoff_without_mass_calib});
      return;
    }
  }
  my $premean = $sum/$n;
	
  my $prestdev = $utils->calc_stdev($sum, $sum2, $n);
	print "Accepting $n high-quality outfiles as internal standards (search score>$$paramhash{'XCorr'} and fully tryptic)\n";
	print LOGFILE "Accepting $n high-quality outfiles as internal standards (search score>$$paramhash{'XCorr'} and fully tryptic)\n";
                                                                                                                                                             
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
        $sum += $ppm; $sum2 += $ppm**2; $n++;
      }
    }
    $premean = $sum/$n;
    $prestdev = $utils->calc_stdev($sum, $sum2, $n);
    @vgarray = @array;#outfiles sorted by scan number
  }
  my $ppmshift = $sum/$n;
  my $stdev = $utils->calc_stdev($sum, $sum2, $n);
  printf "Computing ppm based on $n internal standards: %.2f +/- %.2f\n",$ppmshift,$stdev;
  printf LOGFILE "Computing ppm based on $n internal standards: %.2f +/- %.2f\n",$ppmshift,$stdev;

  # Calculate ppmshift using moving average +/- specified number of scan events
  #print Dumper(\@vgarray);
  for (my $i=$bf; $i < $n-$bf; $i++)
  {
    #establish %smallhash: to store 'good' (non-outlier) scans in one sliding window
    my %smallhash;
    my $outfile = $vgarray[$i];
    next if (defined($goodhash{$outfile}{'outlier'}));
    %{$smallhash{$outfile}} = %{$goodhash{$outfile}};
    $smallhash{$outfile}{'outfile'} = $outfile;
    for (my $j=$i-$bf; $j <= $i+$bf; $j++)
    {
      next if (defined($goodhash{$vgarray[$j]}{'outlier'}));
      %{$smallhash{$vgarray[$j]}} = %{$goodhash{$vgarray[$j]}};
      $smallhash{$vgarray[$j]}{'outfile'} = $outfile;
    }

    #Pick up outliers accoding to Dixon's Q test; rm outliers from %smallhash; set status in %runhash and %goodhash
    my %outliers;
    remove_outliers(\%outliers, \%smallhash, "ppm");
    for my $outlier (keys %outliers)
    {
      $$runhash{$outlier}{'status'} = -1;
      $goodhash{$outlier}{'outlier'} = 1;
    }

    #calculate mean ppm in this sliding window
    my ($ppmsum, $num) = ($goodhash{$outfile}{'ppm'}, 1);
    for my $goodfile (keys %smallhash)
    {
      $num++; $ppmsum += $smallhash{$goodfile}{'ppm'};
    }
    my $mean = $ppmsum/$num;

    #correctedppm & ppmshift
    $goodhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'ppm'} - $mean;
    $goodhash{$outfile}{'ppmshift'} = $mean;

    $$runhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'correctedppm'};
    $$runhash{$outfile}{'ppmshift'} = $mean;
    $$runhash{$outfile}{'status'} = 1;
  }

  # Correct ppm shift for front and rear end of data
  #print Dumper(\%goodhash);exit;
  my $shift = $goodhash{$vgarray[$bf]}{'ppmshift'};
  for (my $i=0; $i<$bf; $i++){
    my $outfile = $vgarray[$i];
    next if (defined($goodhash{$outfile}{'outlier'}));
		#print Dumper(\%goodhash);
    $goodhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'ppm'} - $shift;
    $goodhash{$outfile}{'ppmshift'} = $shift;
		$$runhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'correctedppm'};
    $$runhash{$outfile}{'ppmshift'} = $shift;
  }
  $shift = $goodhash{$vgarray[$n-$bf-1]}{'ppmshift'};
	if ($n>$bf){
  	for (my $i=$n-$bf; $i<$n; $i++){
    	my $outfile = $vgarray[$i];
    	next if (defined($goodhash{$outfile}{'outlier'}));
    	$goodhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'ppm'} - $shift;
    	$goodhash{$outfile}{'ppmshift'} = $shift;
			$$runhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'correctedppm'};
    	$$runhash{$outfile}{'ppmshift'} = $shift;
  	}
	}

  # Calculate sd to use as real data filtering cutoff
  my ($aftersum, $aftersum2, $aftern) = (0,0,0);
  while (my ($outfile, $hash) = each %goodhash)
  {
    next if (defined($goodhash{$outfile}{'outlier'}));
    $aftersum += $$hash{'correctedppm'}; $aftersum2 += $$hash{'correctedppm'}**2; $aftern++;
  }
  printf "Final ppm standard deviation: %.2f (after moving average correction)\n", $utils->calc_stdev($aftersum, $aftersum2, $aftern);
  printf LOGFILE "Final ppm standard deviation: %.2f (after moving average correction)\n", $utils->calc_stdev($aftersum, $aftersum2, $aftern);

	#calculate ppm cutoff using corrected SD
	my $finalSD=$utils->calc_stdev($aftersum, $aftersum2, $aftern);
	my $sdcutoff = $$paramhash{'sd'} * $utils->calc_stdev($aftersum, $aftersum2, $aftern);
	if ($$paramhash{'sd_or_static'} eq "static"){
		$sdcutoff = $$paramhash{'static_cutoff'};
	}
	printf "Mass accuracy cutoff ($$paramhash{sd} SD) = \+\/\- %.2f ppm\n",$sdcutoff;
	printf LOGFILE "Mass accuracy cutoff ($$paramhash{sd} SD) = \+\/\- %.2f ppm\n",$sdcutoff;
	
	# Make small array to help with real data filtering
	my @filter; # @filter=(good_scan_NO)
	my %filterhash;# %filterhash{$good_scan_NO}=ppm_value
	my $last_goodscan = 0;
	for my $outfile (sort {$goodhash{$a}{'scan'}<=>$goodhash{$b}{'scan'}} keys %goodhash)
	{
		next if (defined($goodhash{$outfile}{'outlier'}));
		$last_goodscan = $goodhash{$outfile}{'scan'};
		$filterhash{$goodhash{$outfile}{'scan'}} = $goodhash{$outfile}{'ppmshift'};
		push (@filter, $goodhash{$outfile}{'scan'});
	}

	# Filter the data that did not pass the good cutoff
	my $prev_scan = 0; my $forward_scan = shift @filter;
	my ($prev_shift, $forward_shift) = (0, $filterhash{$forward_scan});
	my %TDhash;
	for my $outfile (sort {$$runhash{$a}{'scan'}<=>$$runhash{$b}{'scan'}} keys %$runhash)
	{
		#next if (abs($$runhash{$outfile}{'status'}) == 1); Filter all including good ones
		my $scan = $$runhash{$outfile}{'scan'}; my $ppm = $$runhash{$outfile}{'ppm'};
		my $ppmshift = 0;

		#calculate ppmshift using 'best' scans (selected from previous processes)
		if ($scan > $last_goodscan)
		{
			$ppmshift = $filterhash{$last_goodscan};
		} 
		else 
		{
			while ($scan > $forward_scan)
			{
				$prev_scan = $forward_scan; $prev_shift = $forward_shift;
				$forward_scan = shift @filter; $forward_shift = $filterhash{$forward_scan};
			}
			$ppmshift = ($prev_shift+$forward_shift)/2;
			$ppmshift = $forward_shift if ($prev_scan == 0);
		}
		if (!defined($ppm) || !defined($ppmshift) || !defined($sdcutoff))
		{
			print Dumper($$runhash{$outfile});exit;
                        print LOGFILE Dumper($$runhash{$outfile});exit;
		}

		# loop to calculate T & D for 1-10 SD / ppm
		#for (my $i=2; $i<=20; $i+=2)
		for (my $i=1; $i<=10; $i++)
		{
			my $tmpCutoff;
			if ($$paramhash{sd_or_static} eq 'sd') { $tmpCutoff=$finalSD*$i;}
			else {$tmpCutoff=$i;}

			if (!defined($TDhash{$i})) { $TDhash{$i}{rm}{D}=$TDhash{$i}{rm}{T}=$TDhash{$i}{ps}{D}=$TDhash{$i}{ps}{T}=0; }

			if (abs($ppm-$ppmshift)>$tmpCutoff)
			{
				if ($$runhash{$outfile}{protein} =~ /Decoy/) { $TDhash{$i}{rm}{D}++; }
				else { $TDhash{$i}{rm}{T}++; }
			}
			else
			{
				if ($$runhash{$outfile}{protein} =~ /Decoy/) { $TDhash{$i}{ps}{D}++; }
				else { $TDhash{$i}{ps}{T}++; }
			}
		}
		

		#filter all scans by the cutoff
		if (abs($ppm-$ppmshift)>$sdcutoff)
		#if (abs($ppm)>$sdcutoff) # for debug
		{
			$$runhash{$outfile}{'status'} = -1;
			$$runhash{$outfile}{'correctedppm'} = $ppm-$ppmshift;
	   		$$runhash{$outfile}{'ppmshift'} = $ppmshift;
		} 
		else 
		{
			$$runhash{$outfile}{'status'} = 1;
			$$runhash{$outfile}{'correctedppm'} = $ppm-$ppmshift;
   			$$runhash{$outfile}{'ppmshift'} = $ppmshift;
		}
		$prev_scan = $scan; $prev_shift = $ppmshift;
	}
	#print T & D counts for 1-10 SD / ppm
	#for (my $i=2; $i<=20; $i=$i+2)
	for (my $i=1; $i<=10; $i++)
	{
		if ($$paramhash{sd_or_static} eq 'sd') {print "  $i SD";}
		else {print "  $i ppm";}
		#$TDhash{$i}{rm}{D}+=0.01;
		printf ": $TDhash{$i}{rm}{T} targets and $TDhash{$i}{rm}{D} decoys deleted (FDR for deleted outfiles = %.2f%%)\n", 200*$TDhash{$i}{rm}{D}/($TDhash{$i}{rm}{D}+$TDhash{$i}{rm}{T}+0.01);
#=head
		if ($$paramhash{sd_or_static} eq 'sd') {print LOGFILE "  $i SD";}
		else {print LOGFILE "  $i ppm";}
		#$TDhash{$i}{rm}{D}+=0.01;
		printf LOGFILE ": $TDhash{$i}{rm}{T} targets and $TDhash{$i}{rm}{D} decoys deleted (FDR for deleted outfiles = %.2f%%)\n", 200*$TDhash{$i}{rm}{D}/($TDhash{$i}{rm}{D}+$TDhash{$i}{rm}{T}+0.01);
#=cut
	}

	#??
	for my $outfile (sort {$$runhash{$a}{'scan'}<=>$$runhash{$b}{'scan'}} keys %$runhash)
	{
		next if (abs($$runhash{$outfile}{'status'}) != 1);
		if ($$runhash{$outfile}{'Mtype'} == -1){ $MX++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 0){ $M0++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 1){ $M1++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 2){ $M2++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 3){ $M3++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 4){ $M4++;
		}
	}
	close LOGFILE;
}

sub staticacumass{
	shift @_;
	my ($run, $runhash, $paramhash, $staticppm) = @_;
	my $masstype = $$paramhash{'mass_consideration'}{'value'};
	my $bf = $$paramhash{'scan_events'};
	#print "       Filtering $run\n";

	# Create a hash of only good outfiles using thresholds set in parameter file
	my ($XCorr, $dCn, $tryptic, $charge) = ($$paramhash{'XCorr'}, $$paramhash{'dCn'}, $$paramhash{'trypticity'}, $$paramhash{'charge'});
	my %goodhash;
	while (my ($outfile, $outhash) = each %$runhash){
		my ($mass, $expmass, $charge) = ($$outhash{'MH'}, $$outhash{'expMH'}, $$outhash{'charge'});
		my ($diff, $status, $Mtype) = get_realdiff($mass, $expmass, $masstype);
		my $ppm = ($diff/$mass/$charge)*(10**6);
    $$outhash{ppm} = $ppm; $$outhash{status} = $status; $$outhash{Mtype} = $Mtype;
		next if ($status == -1);
		if ($$outhash{'XCorr'} >= $XCorr && $$outhash{'dCn'} >= $dCn && $$outhash{'tryptic'} == $tryptic && $$outhash{'charge'} == $charge){


####### changed by xusheng ###########
#
#			next if ($$outhash{'protein'} =~ /Decoy__/ || $$outhash{'protein'} =~ /Decoy/);
			 next if ($$outhash{'protein'} =~ /Random__/ || $$outhash{'protein'} =~ /Decoy__/);
			%{$goodhash{$outfile}} = %$outhash; $$outhash{status} = 1;
		}
	}
	

	# Do preliminary cyclic 2SD filter of good outfiles hash
  my ($sum, $sum2, $n) = (0,0,scalar(keys %goodhash));
  while (my ($outfile, $hash) = each %goodhash){
    my $ppm = $$hash{'ppm'};
    $sum += $ppm; $sum2 += $ppm**2;
  }
  my $premean= $sum/$n;
  my $prestdev = $utils->calc_stdev($sum, $sum2, $n);
	#print "       Automatically accepted MS/MS = $n\n";
	if ($n < $bf){
		while (my ($outfile, $outhash) = each %$runhash){
			$$outhash{status} = 1;
		}
		print "Skipping accurate mass filtering because of not enough good scans.\n";
		return;
	}
                                                                                                                                                             
  my @vgarray;
  for (my $cycle=1; $cycle<=2; $cycle++){
    my @array;
    ($sum, $sum2, $n) = (0,0,0);
    for my $outfile (sort {$goodhash{$a}{'scan'}<=>$goodhash{$b}{'scan'}} keys %goodhash){
      my $ppm = $goodhash{$outfile}{'ppm'};
      if (($ppm < $premean - 2*$prestdev) || ($ppm > $premean + 2*$prestdev)){
				$$runhash{$outfile}{'status'} = -1;
        delete ($goodhash{$outfile});
      } else {
        push (@array, $outfile);
        $sum += $ppm; $sum2 += $ppm**2; $n++;
      }
    }
    $premean = $sum/$n;
    $prestdev = $utils->calc_stdev($sum, $sum2, $n);
    @vgarray = @array;
  }
  my $ppmshift = $sum/$n;
  my $stdev = $utils->calc_stdev($sum, $sum2, $n);
  #print "         Average ppm shift: $ppmshift +/- $stdev --> $n measurements\n";

	# Calculate ppmshift using moving average +/- specified number of scan events
	for (my $i=$bf; $i < $n-$bf; $i++){
    my %smallhash;
    my $outfile = $vgarray[$i];
    next if (defined($goodhash{$outfile}{'outlier'}));
    %{$smallhash{$outfile}} = %{$goodhash{$outfile}};
    $smallhash{$outfile}{'outfile'} = $outfile;
    for (my $j=$i-$bf; $j <= $i+$bf; $j++){
      next if (defined($goodhash{$vgarray[$j]}{'outlier'}));
      %{$smallhash{$vgarray[$j]}} = %{$goodhash{$vgarray[$j]}};
      $smallhash{$vgarray[$j]}{'outfile'} = $outfile;
    }
    my %outliers;
    remove_outliers(\%outliers, \%smallhash, "ppm");
    for my $outlier (keys %outliers){
			$$runhash{$outlier}{'status'} = -1;
      $goodhash{$outlier}{'outlier'} = 1;
    }
    my ($ppmsum, $num) = ($goodhash{$outfile}{'ppm'}, 1);
    for my $goodfile (keys %smallhash){
      $num++; $ppmsum += $smallhash{$goodfile}{'ppm'};
    }
    my $mean = $ppmsum/$num;
    $goodhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'ppm'} - $mean;
    $goodhash{$outfile}{'ppmshift'} = $mean;
		$$runhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'correctedppm'};
    $$runhash{$outfile}{'ppmshift'} = $mean;
    $$runhash{$outfile}{'status'} = 1;
  }

	# Correct ppm shift for front and rear end of data
	my $shift = $goodhash{$vgarray[$bf]}{'ppmshift'};
  for (my $i=0; $i<$bf; $i++){
    my $outfile = $vgarray[$i];
    next if (defined($goodhash{$outfile}{'outlier'}));
    $goodhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'ppm'} - $shift;
    $goodhash{$outfile}{'ppmshift'} = $shift;
		$$runhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'correctedppm'};
    $$runhash{$outfile}{'ppmshift'} = $shift;
  }
  $shift = $goodhash{$vgarray[$n-$bf-1]}{'ppmshift'};
	if ($n>$bf){
  	for (my $i=$n-$bf; $i<$n; $i++){
    	my $outfile = $vgarray[$i];
    	next if (defined($goodhash{$outfile}{'outlier'}));
    	$goodhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'ppm'} - $shift;
    	$goodhash{$outfile}{'ppmshift'} = $shift;
			$$runhash{$outfile}{'correctedppm'} = $goodhash{$outfile}{'correctedppm'};
    	$$runhash{$outfile}{'ppmshift'} = $shift;
  	}
	}

	# Calculate sd to use as real data filtering cutoff
  my ($aftersum, $aftersum2, $aftern) = (0,0,0);
  while (my ($outfile, $hash) = each %goodhash){
    next if (defined($goodhash{$outfile}{'outlier'}));
    $aftersum += $$hash{'correctedppm'}; $aftersum2 += $$hash{'correctedppm'}**2; $aftern++;
  }
  #printf "         Calculated ppm standard deviation: %.2f\n", $utils->calc_stdev($aftersum, $aftersum2, $aftern);
	my $sdcutoff = $$paramhash{'sd'} * $utils->calc_stdev($aftersum, $aftersum2, $aftern);
	$sdcutoff = $staticppm;
	#print "         Mass accuracy cutoff = \+\/\- $sdcutoff ppm\n";
	
	# Make small array to help with real data filtering
	my @filter; my %filterhash;
	my $last_goodscan = 0;
	for my $outfile (sort {$goodhash{$a}{'scan'}<=>$goodhash{$b}{'scan'}} keys %goodhash){
		next if (defined($goodhash{$outfile}{'outlier'}));
		$last_goodscan = $goodhash{$outfile}{'scan'};
		$filterhash{$goodhash{$outfile}{'scan'}} = $goodhash{$outfile}{'ppmshift'};
		push (@filter, $goodhash{$outfile}{'scan'});
	}
	# Filter the data that did not pass the good cutoffs
	my $prev_scan = 0; my $forward_scan = shift @filter;
	my ($prev_shift, $forward_shift) = (0, $filterhash{$forward_scan});
	for my $outfile (sort {$$runhash{$a}{'scan'}<=>$$runhash{$b}{'scan'}} keys %$runhash){
		next if (abs($$runhash{$outfile}{'status'}) == 1);
		my $scan = $$runhash{$outfile}{'scan'}; my $ppm = $$runhash{$outfile}{'ppm'};
		my $ppmshift = 0;
		if ($scan > $last_goodscan){
			$ppmshift = $filterhash{$last_goodscan};
		} else {
			while ($scan > $forward_scan){
				$prev_scan = $forward_scan; $prev_shift = $forward_shift;
				$forward_scan = shift @filter; $forward_shift = $filterhash{$forward_scan};
			}
			$ppmshift = ($prev_shift+$forward_shift)/2;
			$ppmshift = $forward_shift if ($prev_scan == 0);
		}
		if (abs($ppm-$ppmshift)>$sdcutoff){
			$$runhash{$outfile}{'status'} = -1;
		} else {
			$$runhash{$outfile}{'status'} = 1;
			$$runhash{$outfile}{'correctedppm'} = $ppm-$ppmshift;
   		$$runhash{$outfile}{'ppmshift'} = $ppmshift;
		}
		$prev_scan = $scan; $prev_shift = $ppmshift;
	}
	
	for my $outfile (sort {$$runhash{$a}{'scan'}<=>$$runhash{$b}{'scan'}} keys %$runhash){
		next if (abs($$runhash{$outfile}{'status'}) != 1);
		if ($$runhash{$outfile}{'Mtype'} == -1){ $MX++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 0){ $M0++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 1){ $M1++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 2){ $M2++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 3){ $M3++;
		} elsif ($$runhash{$outfile}{'Mtype'} == 4){ $M4++;
		}
	}
}


sub get_realdiff
{
  my ($mass, $expmass, $masstype) = @_;
  my $mtype = 0;                                                                                                                                                           
  my $realdiff = $mass-$expmass; my $diff = $realdiff;
  my $status = 0;
	if ($masstype==1)
	{
	  	if ($realdiff > -0.4 && $realdiff <= 0.4)
		{   
			$diff = $realdiff;  
		} 
		else 
		{ 
			$status = -1; $mtype = -1;
		}
	} 
	elsif ($masstype==2){
  		if ($realdiff > -0.4 && $realdiff <= 0.4){   $diff = $realdiff;  
  		} elsif ($realdiff > 0.6 && $realdiff <= 1.4){    $diff = $realdiff-$H; $mtype = 1;
		} else { $status = -1; $mtype = -1; }
	} elsif ($masstype==3){
  		if ($realdiff > -0.4 && $realdiff <= 0.4){   $diff = $realdiff;  
  		} elsif ($realdiff > 0.6 && $realdiff <= 1.4){    $diff = $realdiff-$H; $mtype = 1;
  		} elsif ($realdiff > 1.6 && $realdiff <= 2.4){    $diff = $realdiff-2*$H; $mtype = 2;
		} else { $status = -1; $mtype = -1;}
	} elsif ($masstype==4){
  		if ($realdiff > -0.4 && $realdiff <= 0.4){   $diff = $realdiff;  
  		} elsif ($realdiff > 0.6 && $realdiff <= 1.4){    $diff = $realdiff-$H; $mtype = 1;
  		} elsif ($realdiff > 1.6 && $realdiff <= 2.4){    $diff = $realdiff-2*$H; $mtype = 2;
  		} elsif ($realdiff > 2.6 && $realdiff <= 3.4){    $diff = $realdiff-3*$H; $mtype = 3;
		} else { $status = -1; $MX++; $mtype = -1}
	} elsif ($masstype==5){
  		if ($realdiff > -0.4 && $realdiff <= 0.4){   $diff = $realdiff;  
  		} elsif ($realdiff > 0.6 && $realdiff <= 1.4){    $diff = $realdiff-$H; $mtype = 1;
  		} elsif ($realdiff > 1.6 && $realdiff <= 2.4){    $diff = $realdiff-2*$H; $mtype = 2;
  		} elsif ($realdiff > 2.6 && $realdiff <= 3.4){    $diff = $realdiff-3*$H; $mtype = 3;
  		} elsif ($realdiff > 3.6 && $realdiff <= 4.4){    $diff = $realdiff-4*$H; $mtype = 4;
		} else { $status = -1; $mtype = -1;}
	} elsif ($masstype==6){
  		if ($realdiff > -1.4 && $realdiff <= -0.6){   $diff = $realdiff+$H; $mtype = 6;  
  		} elsif ($realdiff > -0.4 && $realdiff <= 0.4){   $diff = $realdiff;  
  		} elsif ($realdiff > 0.6 && $realdiff <= 1.4){    $diff = $realdiff-$H; $mtype = 1;
  		} elsif ($realdiff > 1.6 && $realdiff <= 2.4){    $diff = $realdiff-2*$H; $mtype = 2;
		} else { $status = -1; $mtype = -1;}

	} elsif ($masstype==7){
  		if ($realdiff > -2.4 && $realdiff <= -1.6){   $diff = $realdiff+2*$H; $mtype = 7;  
  		} elsif ($realdiff > -1.4 && $realdiff <= -0.6){   $diff = $realdiff+$H; $mtype = 6;  
  		} elsif ($realdiff > -0.4 && $realdiff <= 0.4){   $diff = $realdiff;  
  		} elsif ($realdiff > 0.6 && $realdiff <= 1.4){    $diff = $realdiff-$H; $mtype = 1;
  		} elsif ($realdiff > 1.6 && $realdiff <= 2.4){    $diff = $realdiff-2*$H; $mtype = 2;
		} else { $status = -1; $mtype = -1;}
	}
  return ($diff, $status, $mtype);
}

sub remove_outliers
{
  my ($outlier, $hash, $type) = @_;

  #store outhash address in @array; out hashes are sorted by PPM value  
  my @array;
  for my $outhash (sort {$$a{$type} <=> $$b{$type}} values %$hash){
    push(@array, $outhash);
  }

  my ($round, $removed) = (1, 1);
  while ($removed >= 1)
  {
    $removed = 0;
    my $n = scalar(@array); next if ($n < 3);
    my $q = 0;
    if ($n < 8){ $q = $dixon{$n}{$n}{'95'};} else {$q = $dixon{20}{$n}{'95'}; }
    my ($X1, $XN) = ($array[0], $array[$n-1]);
    my ($low, $high) = calc_critical_value($n, \@array, $type);#calculate Q=gap/range for the two ends
    if ($low > $q){ $removed++; shift @array; delete($$hash{$$X1{'outfile'}}); $$outlier{$$X1{'outfile'}} = 1;}
    if ($high > $q){ $removed++;pop @array; delete($$hash{$$XN{'outfile'}}); $$outlier{$$XN{'outfile'}} = 1;}
  }
}

sub calc_critical_value
{
  my ($num, $array, $type) = @_;
  my ($low, $high) = (0,0);

#for (my $i=0;$i<$num;$i++) {print "$$array[$i]{$type} ";} print "\n";

  if ($num>=8){
    my ($X1, $X3, $XN, $XNmin2) = ($$array[0], $$array[2], $$array[$num-1], $$array[$num-3]);
		if ($$XNmin2{$type}-$$X1{$type}==0 || $$XN{$type}-$$X3{$type} == 0){
			$low = 0; $high = 0;
		} else {
    	$low = ($$X3{$type}-$$X1{$type})/($$XNmin2{$type}-$$X1{$type});
    	$high = ($$XN{$type}-$$XNmin2{$type})/($$XN{$type}-$$X3{$type});
		}
  } elsif($num>=3 && $num<=7){
    my ($X1, $X2, $XN, $XNmin1) = ($$array[0], $$array[1], $$array[$num-1], $$array[$num-2]);
		if ($$XN{$type}-$$X1{$type} == 0){
			$low = 0; $high = 0;
		} else {
    	$low = ($$X2{$type}-$$X1{$type})/($$XN{$type}-$$X1{$type});
    	$high = ($$XN{$type}-$$XNmin1{$type})/($$XN{$type}-$$X1{$type});
		}
  }
                                                                                                                                                             
  return ($low, $high);
}

sub calc_emPAI{
	shift @_;

	my ($paramhash, $proteinhash) = @_;

	for my $protein (keys %$proteinhash){
		my $ftpeps = $utils->get_peptides_emPAI($$proteinhash{$protein}{'sequence'}, $$paramhash{'min_peptide_length'}, 
																			$$paramhash{'enzyme'}, $$paramhash{'max_peptide_length'});
		my $good = 0;
		for my $peptide (@$ftpeps){
			$peptide =~ s/[\-A-Z]+\.([A-Z]+)\.[A-Z\-]+\Z/$1/;
			my @temp = split("", $peptide);
			my $pho = $utils->get_Pho(\@temp);
			$good++ if (($pho <= $$paramhash{'max_peptide_hydro'}) && ($pho >= $$paramhash{'min_peptide_hydro'}));
		}
		$$proteinhash{$protein}{'emPAI_ftpeps'} = $good;
	}
}


sub parse_JUMPparams
{
	shift @_;
        my ($seqparamhash, $fraction) = @_;


        open (LOGFILE, ">>$log_file");                                                                                         
	#select((select(LOGFILE), $|=1)[0]);
        my $paramfile = "$fraction\/jump.params";
	unless (-e $paramfile) { $paramfile = 'jump.params'; }

        while (!(-e $paramfile)){
                print "\n$paramfile does not exists!!!";
                print LOGFILE "\n$paramfile does not exists!!!";
                $paramfile = idsum2::CL->ask("\nPlease enter the correct full path to the jump.params parameter file.", "");
        }

        open (IN, "<$paramfile");

        while (<IN>)
        {
                next if (/^#/ || /^\s+\Z/);
                chomp;
                s/[;#].*//;
                s/\s+\Z//;              #Q: remove empty lines
                my ($key, $value) = split(" = ", $_);
                next if (!defined($value));

                if ($key =~ /^dynamic_([A-Z]+)$/)       #dynamic mod
                {
                        $$seqparamhash{$fraction}{"dynamicmods"}{$1} = $value;
                        next;
                }
                if ($key eq 'add_Cterm_peptide')        #cterm mod
                {
                        $$seqparamhash{$fraction}{"staticmods"}{'cterm'} = $value;
                        next;
                }
                if ($key eq 'add_Nterm_peptide')        #nterm mod
                {
                        $$seqparamhash{$fraction}{"staticmods"}{'nterm'} = $value;
                        next;
                }
                #add_A_Alanine
                if ($key =~ /^add_([A-Z])_(\w+)$/)      #static mod
                {
                        $$seqparamhash{$fraction}{"staticmods"}{$1} = $value;
                        next;
                }

        }

        close IN;

}

sub pepXMLhashesTransfer
{
  shift @_;
  my ($paramhash,$run,$outfile,$paraHash,$outInforHash,$peptidehash,$runhash,$delhash,$blank,$cpr,$cpn,$pep2reads) = @_;
  #my ($paramhash,$peptidehash,$runhash,$run,$outfile,$delhash,$blank,$cpr,$cpn,$pep2reads) = @_;

  my ($Hydrogen_mass, $HydrogenMinus_mass, $Cterm_mass)=(1.007825032,1.00335,19.017806);
  
  my %proteins;
  my $fullpath = "$run\/$outfile";

  # check if RNAseq data is vailable
  my $noRNAseq;
  if (!defined($cpr) || !defined($cpn) || !defined($pep2reads)) { $noRNAseq=1; } else { $noRNAseq=0; }
 
  #open (IN, $fullpath);
  #u170k_siRNA_02ug.14423.1.3.spout
  $outfile =~ m/\.(\d+)\.(\d+)\.(\d)\.spout/;
  my ($scan,$charge) = ("$1.$2",$3);
  my $ppi=$2; $ppi=substr($ppi,0,1);

  # check if filter ppi
  if (defined($$paramhash{ppi_included}))
  {
        my @t=split(/_/,$$paramhash{ppi_included});
        my $incl=0;
        for (my $i=0; $i<=$#t; $i++) { if ( $t[$i] == $ppi ) {$incl=1;last;} }
        unless ($incl) { $$delhash{ppi_filtered}++; return }
  }

  #my ($sorcerer) =  grep(/SageNResearch/ || /SORCERER/, <IN>);
  #if (defined($sorcerer)){ $$paramhash{'sorcererouts'} = 1; }
  my $outkey=$outfile; $outkey =~ s/(\.spout)//;  #print "$outkey\n";

  #get PSM lines
  #seek(IN, 0,0);
  #my @hits=grep(/^\s+\d+\.\s+\d+\s+\//,<IN>);
  if (!defined($$outInforHash{$outkey}{hitShown}) || scalar($$outInforHash{$outkey}{hitShown})==0) { $$blank++; return; }

  #get red protain names
  #seek(IN,0,0);
  #my $red_count=0;
  #my @redProNames=grep(/^\s+\d+\s+[0-9A-Za-z#]/,<IN>);
  #close IN;

  #decide which PSM and protein name to use
  my ($mxlv,$mxtd,$mxG,$mxRpk,$dCn,$mxProName,$mxRank,$mxReadSpt)=(10,'','',0,0,'',0,0,0);
  my @tmpProteins; my $oldPep='';
  for (my $i=0; $i<$$outInforHash{$outkey}{hitShown}; $i++)
  {
    #remove redundant space
    #$hits[$i] =~ s/\s*\/\s*/\//g; $hits[$i] =~ s/(\s)\s+/$1/g; $hits[$i] =~ s/^\s+//; $hits[$i] =~ s/\s+\Z//;
    #my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $hits[$i]);
    #$dCn=$c;
    #my ($peptide, $red)  = ($g, 0);
    #if(defined($h)) { $red = $g; $red =~ s/\+//; $peptide = $h; }

    # pass values from %outInforHash
    my $k=$i+1;# corespond to %outInforHash
    my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein) = ($k,'','',$$outInforHash{$outkey}{'calc_neutral_pep_mass'}[$k]+$Hydrogen_mass,0,$$outInforHash{$outkey}{'WeightedEvalue'}[$k],'',"$$outInforHash{$outkey}{'num_matched_ions'}[$k]\/$$outInforHash{$outkey}{'tot_num_ions'}[$k]",$$outInforHash{$outkey}{'protein'}[$k]);
    if ($i==0) { $dCn=0; } else { $dCn=$$outInforHash{$outkey}{'deltacn'}[$i]; }

	#print "$outkey,$i,$dCn\n";
    #if ($i==0) { $dCn=0; } else { $dCn=$$outInforHash{$outkey}{'deltacn'}[$i]; }
    my ($peptide, $red)  = ("$$outInforHash{$outkey}{'peptide_prev_aa'}[$k]\.$$outInforHash{$outkey}{'fullPeptide'}[$k]\.$$outInforHash{$outkey}{'peptide_next_aa'}[$k]",$$outInforHash{$outkey}{'num_tot_proteins'}[$k]-1);

    #last if ( $i>=1 and $dCn>0 ); # old critera

    # criteria for jumping out the loop
    last if ( $i>=1 and $dCn>0.05 ); 
    my $pepSim=1;
    if ( $oldPep ne '' ) { $pepSim=peptide_similarity($oldPep,$peptide); }
	#if ( $pepSim<0 ) {die "uneuqla length for peptide sequences with small dCn!!!\n $outfile,$oldPep,$peptide\n";}
    last if ( $pepSim<0.9 ); 
    last if ( $i>=10 );

    # for 1st protein name
    my $Decoy=0;
    if ($noRNAseq) 
    { 
      if ($i==0) 
      { 
        $mxProName=$protein; 
        if ( $protein =~ /Decoy/ ) { $Decoy=1; }
      } 
    }
    else {_update_max_levels(\$mxlv,\$mxtd,\$mxG,\$mxRpk,\$mxProName,\$mxRank,$protein,$cpr,$cpn,$i);}

    if ( $protein =~ /\#\#Decoy/ ) { $protein =~ s/\#\#//; }

    # for redundant protein names
    my $j;
    for ( $j=1; $j<=$red; $j++)
    {
	#if (!defined($redProNames[$red_count])) {die "undefined red proten names lines:\n $outfile,$i,$j,$red,$red_count\n";}
      #chomp($redProNames[$red_count]); my $tmpP=$redProNames[$red_count]; $red_count++;
      ##$tmpP =~ s/^\s+\d+\s+//;  if ( $tmpP =~ /\s/ ) { $tmpP =~ s/^(.*?)\s+/$1/; }   $tmpProteins[$i][$j]=$tmpP;# old
      #if($tmpP=~/^\s+\d*\s*([\#a-zA-Z0-9\.\_\-\|\:]+)\s*/) { $tmpProteins[$i][$j]=$1;} else {die "Ero:red pro\n";}

      my $tmpP=$$outInforHash{$outkey}{'alternative_protein'}[$j][$k]; $tmpProteins[$i][$j]=$tmpP;
#  print "$outfile,$i,$k,$j,$tmpP,$tmpProteins[$i][$j]\n";

      if ($noRNAseq) 
      {
	if ($Decoy and $tmpP !~ /Decoy/) {   $mxProName=$tmpP;}
      }
      else {_update_max_levels(\$mxlv,\$mxtd,\$mxG,\$mxRpk,\$mxProName,\$mxRank,$tmpP,$cpr,$cpn,$i);}
    }
    $tmpProteins[$i][0]=$j-1;

    $oldPep=$peptide;

    if ($noRNAseq) {}else {
    my $readSpt=get_read_spt($peptide,$protein,$paramhash,$pep2reads);
#if ( $outfile eq 'u170k_siRNA_02ug.62439.62439.3.out' ) { print "$outfile,$i,$mxReadSpt,$readSpt,$peptide,$protein\n"; }
    if ($mxlv==4 and $mxReadSpt<$readSpt )
    { $mxReadSpt=$readSpt; $mxRank=$i; $mxProName=$protein; }}
  }
  if ($noRNAseq) { $mxlv=1; $mxRank=0;}
  if ( $$outInforHash{$outkey}{hitShown}==1 ) { $dCn=1; }
  my $final_rank=$mxRank;
#  for (my $j=1; $j<=$tmpProteins[$final_rank][0]; $j++) {print "$outfile,$final_rank,$j,$tmpProteins[$final_rank][$j]\n";}

  #continue to use old code, but make an interface
  #open (IN, $fullpath);
  #seek(IN, 0,0);
  #mass and database
=head
  my ($mass, $database) = grep(/\+\smass\s\=/ || /proteins\s\=/, <IN>);
  my $temp=$hits[$final_rank];#interface
  $temp =~ s/\s*\/\s*/\//g; $temp =~ s/(\s)\s+/$1/g; $temp =~ s/^\s+//; $temp =~ s/\s+\Z//;
  $mass =~ s/mass\s+\=\s+([\d\.]+)\s//; $mass = $1;
  $database =~ s/\s+//g; my @temparray = split (/\,/, $database); $database = pop(@temparray);
  if($database=~/hdr/)	{$database=~s/\.hdr$//;	}
=cut
  my ($mass, $database) = ($$outInforHash{$outkey}{'calc_neutral_pep_mass'}[1]+$Hydrogen_mass+$$outInforHash{$outkey}{'massdiff'}[1],$$paraHash{'DB'});
  if($database=~/mdx/) {$database=~s/\.mdx$//;}

  # PSM line
  my $k=$final_rank+1;
  my $tag=''; if (defined($$outInforHash{$outkey}{'TagSeq'}[$k])) {$tag=$$outInforHash{$outkey}{'TagSeq'}[$k];}
  #my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $temp);
  my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e) = ($k,'','',$$outInforHash{$outkey}{'calc_neutral_pep_mass'}[$k]+$Hydrogen_mass,0,$$outInforHash{$outkey}{'WeightedEvalue'}[$k],'',"$$outInforHash{$outkey}{'num_matched_ions'}[$k]\/$$outInforHash{$outkey}{'tot_num_ions'}[$k]");
  my $protein=$mxProName;
  my $precursor_peak_intensity_percentage=$$outInforHash{$outkey}{'precursor_peak_intensity_percentage'}; 
  #print "$outkey $precursor_peak_intensity_percentage\n";

  #if (defined($sorcerer)){    $$paramhash{'sorcererouts'} = 1;  }
  
  if ($dCn < $$paramhash{'min_dCn'} || $XCorr < $$paramhash{'min_XCorr'}){
    $$delhash{'DX'}++;
    if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'DX'}++; } else { $$delhash{'target'}{'DX'}++; }
    return ($database);
  }
  #my $peptide = $g; my $red = 0; 

#my $ions = $e;
my $ions = "$$outInforHash{$outkey}{'num_matched_ions'}[$k]\/$$outInforHash{$outkey}{'tot_num_ions'}[$k]";
#if (!defined($protein)) {die "$outfile,$protein,$mxRank,$noRNAseq,$tmpProteins[0][1]\n";}
#  $protein =~ s/\,.*//;

  if (defined($$paramhash{'filter_contaminants'})){
    if ($$paramhash{'filter_contaminants'} == 1){
#      if (defined($$paramhash{'badproteins'}{$protein}) || $protein =~ /keratin/ || $protein=~ /KERATIN/){
     if (defined($$paramhash{'badproteins'}{$protein}) || $protein =~ /keratin/ || $protein=~ /CON_/){ 
       $$delhash{'contaminants'}++;
       if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'contaminants'}++; } else { $$delhash{'target'}{'contaminants'}++; }
        return ($database);
      }
    }
  }
  $proteins{$protein} = 1;
  my $target_decoy_mixed = 0;

  my ($peptide, $red)  = ("$$outInforHash{$outkey}{'peptide_prev_aa'}[$k]\.$$outInforHash{$outkey}{'fullPeptide'}[$k]\.$$outInforHash{$outkey}{'peptide_next_aa'}[$k]",$$outInforHash{$outkey}{'num_tot_proteins'}[$k]-1);
  undef($k);
  if($red>0)
  {
    #if(defined($h) && $protein !~ /Random__/ && $protein !~ /Decoy__/  && $protein !~/CON_/){
    #$red = $g; $red =~ s/\+//; $peptide = $h;
    for (my $j=1; $j<scalar(@{$tmpProteins[$final_rank]}); $j++) 
    {
      my $pro=$tmpProteins[$final_rank][$j];#print "$outfile,$final_rank,$j,$tmpProteins[$final_rank][$j]\n";

	#Remove the decoy in the redundancy??
        #mark targets in %proteins
        #if ($pro !~ /Random__/ && $pro !~ /Decoy__/ && $pro !~ /CON_/){
        if ($pro !~ /CON_/ && $pro !~ /Decoy/)
        {
                if (defined($$paramhash{'filter_contaminants'})){
                        if ($$paramhash{'filter_contaminants'} == 1){
                                if (!defined($$paramhash{'badproteins'}{$pro}) && $pro !~ /keratin/ && $pro !~ /KERATIN/ && $pro !~ /CON_/){
                                        $proteins{$pro} = 1;
                                }
                        }
                        else {
                                $proteins{$pro} = 1;
				#print "$outfile,$pro\n";
                        }
                }
        }
        if($pro =~ /Decoy/)
        {
                $target_decoy_mixed = 1;
        }
	last if ($j>=20);
    }

  }

  if($protein !~ /Decoy/ && $target_decoy_mixed == 1)
  {
	$$delhash{'sharedecoy'}++;
  }
  my $testpeptide = $peptide;
  my $mods = ($testpeptide =~ s/([\@\%\&\^\~\$\#\*\^\~]+)/$1/g) || 0;
  $testpeptide =~ s/.\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\../$1/;
  my $mis = ($testpeptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;
  if (defined($$paramhash{'max_peptide_mod'})){
    if ($mods > $$paramhash{'max_peptide_mod'}){
      #$$max_mod_del++;
      
      $$delhash{'maxmod'}++;
      if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'maxmod'}++; } else { $$delhash{'target'}{'maxmod'}++; }
      return ($database);
    }
  }

  if (defined($$paramhash{'max_peptide_mis'})){
    if ($mis > $$paramhash{'max_peptide_mis'}){
      #$$max_mis_del++;
      #print $mis."mismismis\t".$$paramhash{'max_peptide_mis'}."hashmishashmis\t".$max_mis_del."\tccccccccccc\n";
      $$delhash{'maxmis'}++;
      if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'maxmis'}++; } else { $$delhash{'target'}{'maxmis'}++; }
      return ($database);
    }
  }
  my $tryptic = $utils->istryptic($peptide);#returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
  $tryptic -= 1 if ($tryptic>1);
  my $intpep = $peptide;
  $intpep =~ s/[A-Z\-]\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\.[A-Z\-]/$1/;
  if ($$paramhash{'mix_label'} ne 0){  #???
    my $labels = "\[$$paramhash{'mix_label'}\]";
    if ($intpep =~ /$labels[\*\#\@\%\&\^\~\$\^\~]/ && ($intpep =~ /$labels[A-Z]/ || $intpep =~ /$labels\Z/)){
        $$delhash{'mixlabel'}++;
        if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'mixlabel'}++; } else { $$delhash{'target'}{'mixlabel'}++; }
				#print "$outfile $peptide $XCorr $dCn\n";
        return ($database);
    }
  }
  if (defined($$paramhash{'peptide_mod_removal'})){
    if ($$paramhash{'peptide_mod_removal'} ne 0){
      my @labels = split("", $$paramhash{'peptide_mod_removal'});
      for my $aa (@labels){
        my $label = "\[$aa\]";
        if ($intpep =~ /$label[\*\#\@\%\&\^\~\$\^\~]/){
          $$delhash{'modremoval'}++;
          if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'modremoval'}++; } else { $$delhash{'target'}{'modremoval'}++; }
          return ($database);
        }
      }
    }
  }
  if (defined($$paramhash{'peptide_aa_removal'})){
    if ($$paramhash{'peptide_aa_removal'} ne 0){
      my @labels = split("", $$paramhash{'peptide_aa_removal'});
      for my $aa (@labels){
        my $label = "\[$aa\]";
        if ($intpep =~ /$label/){
          $$delhash{'aaremoval'}++;
          if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'aaremoval'}++; } else { $$delhash{'target'}{'aaremoval'}++; }
          return ($database);
        }
      }
    }
  }
  my $nomod = $intpep; $nomod =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;
  if (length($nomod) < $$paramhash{'min_peptide_length'}){
    $$delhash{'length'}++;
    if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'length'}++; } else { $$delhash{'target'}{'length'}++; }
    return ($database);
  }
  my @array = split("", $nomod);
  my $calcmw = $utils->get_MW(\@array); #returns molecular weight of given sequence
  my $num = 0;
	if (defined($$paramhash{$run}{staticmods}{'nterm'})){
    $calcmw += $$paramhash{$run}{staticmods}{'nterm'};
  } elsif (defined($$paramhash{$run}{staticmods}{'cterm'})){
    $calcmw += $$paramhash{$run}{staticmods}{'cterm'};
  }
  for my $aa (keys %{$$paramhash{$run}{staticmods}}){
    $num += ($intpep =~ s/$aa/$aa/g); #print "$outfile: $peptide, $aa, $num\n";
    $calcmw += $num * $$paramhash{$run}{staticmods}->{$aa};
    $num = 0;
  }
  for my $aa (keys %{$$paramhash{$run}{dynamicmods}}){
    my $test = "\[$aa\]";
    $num += ($intpep =~ s/($test[\@\%\&\^\~\$\#\*\^\~])/$1/g);
    $calcmw += $num * $$paramhash{$run}{dynamicmods}->{$aa};
    $num = 0;
  }

#	if ($expmass != $calcmw) { die "$outfile: $expmass, $calcmw, $peptide\n"; }	
  $expmass = $calcmw;#$expmass predefined; why bothers recalculation? why diff?
  if ($$paramhash{'norm_XCorr'}){
    my $length = length($nomod);
    if ($charge == 2){
      if ($length >= 15){
        $XCorr = log($XCorr)/log(15*2);
      } else {
        $XCorr = log($XCorr)/log($length*2);
      }
    } elsif ($charge == 3){
      if ($length >= 25){
        $XCorr = log($XCorr)/log(25*4);
      } else {
        $XCorr = log($XCorr)/log($length*4);
      }
    }
  }

#print "XCorr:",$XCorr,"\n";
#print "Sp:",$Sp,"\n";
#print "rank:",$rank,"\n";
#print "dCn:",$dCn,"\n";
#print "protein:",$protein,"\n";
#print "MH:",$mass,"\n";

###### added by yanji
  # get scan number from outfile name
#  my @outfile_element = split/\./, $outfile;
#  pop(@outfile_element);
#  pop(@outfile_element);
#  my $scan_number = pop(@outfile_element);

  # get original hash, origmsms_hash
#  my $hash_dir = dirname($run)."\/".".hashes";
#  my $hash_ref = retrieve("$hash_dir/origmsms_hash");
#  my %origmsms_hash = %$hash_ref;

  # get rention time, rt, and intensity from the hash table
#  $scan_number = int($scan_number || 0);
#  my $rt = $origmsms_hash{$scan_number}{'rt'};
#  my $intensity = $origmsms_hash{$scan_number}{'intensity'};
###############added by xusheng on 06/14/2012 ###########
#  my $prec_int = $origmsms_hash{$scan_number}{'prec_int'};

  # save rt and intensity
#  $$runhash{$run}{$outfile}{'rt'} = $rt;
#  $$runhash{$run}{$outfile}{'intensity'} = $intensity;
#  $$runhash{$run}{$outfile}{'prec_int'} = $prec_int;

	
	if (defined($$peptidehash{$intpep})){ ##ADDED DMD NOVEMBER 25, 2008 to remove peptide redundancy (decoys) for peptides with real and decoy matches
		my $decoy = 0;
		my $target = 0;
		for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){
		#	print $intpep,"\t",$pros,"\n";
########## changed by xusheng #####################
#Remove all start test ^ sign from the statement
#			$decoy++ if ($pros =~ /Decoy/ || $pros =~ /Decoy__/);
#
			$target ++ if ($pros !~ /Decoy__/ && $pros !~ /Random__/);
			$decoy++ if ($pros =~ /Decoy__/ || $pros =~ /Random__/);
		}

		if ($protein !~ /Decoy__/ && $protein !~ /Random__/ && $decoy > 0){#current out is target; but previous have decoys
			for my $delout (keys %{$$peptidehash{$intpep}{'outfiles'}}){#deltete all previous outfiles of this peptide
				my $run = $$peptidehash{$intpep}{'outfiles'}{$delout}{'run'};
				delete $$runhash{$run}{$delout};
			}
			delete $$peptidehash{$intpep};#delete previous item of this peptide (for renew?)
		} elsif (($protein =~ /Decoy__/ || $protein =~ /Random__/) && $decoy == 0) {#current decoy; deltete current out
			delete $$runhash{$run}{$outfile};
			return;
		}

                for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){#no error??? delete $$peptidehash{$intpep};???
                        if($target>0 and ($pros =~ /Decoy__/ || $pros =~ /Random__/))
                        {
                                $$delhash{'sharedecoy'}++;
				delete $$peptidehash{$intpep}{'proteins'}{$pros};#rm decoy items;keep target; why not delete outfiles?
                        }
		}

	}



 	$$runhash{$run}{$outfile}{'scan'} = $scan;
        $$runhash{$run}{$outfile}{'run'} = $run;
       $$runhash{$run}{$outfile}{'peptide'} = $peptide;
        $$runhash{$run}{$outfile}{'tryptic'} = $tryptic;
        $$runhash{$run}{$outfile}{'mod'} = $mods;
        $$runhash{$run}{$outfile}{'mis'} = $mis;
        $$runhash{$run}{$outfile}{'red'} = $red;
        $$runhash{$run}{$outfile}{'charge'} = $charge;
        $$runhash{$run}{$outfile}{'XCorr'} = $XCorr;
       $$runhash{$run}{$outfile}{'Sp'} = $Sp;#??
       $$runhash{$run}{$outfile}{'rank'} = $rank;
       $$runhash{$run}{$outfile}{'dCn'} = $dCn;
       $$runhash{$run}{$outfile}{'protein'} = $protein;
       $$runhash{$run}{$outfile}{'ions'} = $ions;
       $$runhash{$run}{$outfile}{'MH'} = $mass;
       $$runhash{$run}{$outfile}{'expMH'} = $expmass;
       $$runhash{$run}{$outfile}{'path'} = $fullpath;
       $$runhash{$run}{$outfile}{'intpep'} = $intpep;
       $$runhash{$run}{$outfile}{'protein_priorty_level'} = int($mxlv);
       $$runhash{$run}{$outfile}{'rpk'} = $mxRpk;
       $$runhash{$run}{$outfile}{'psmRank'} = $final_rank;
       $$runhash{$run}{$outfile}{'tagLength'} = length($tag);

	if ($noRNAseq) {} else {
	if ( $protein =~ /Decoy/ ) 
	{ 
		my $origPep=orig_pep_for_decoy($paramhash,$peptide);#print "$outfile,$peptide,$origPep\n";
		my @aa=split(//,$origPep); pop(@aa); shift(@aa);
		$origPep=join('',@aa);
		if (defined($$pep2reads{$origPep})) {$$runhash{$run}{$outfile}{'readSpt'} = $$pep2reads{$origPep};}
		else { $$runhash{$run}{$outfile}{'readSpt'} = 0; }
	}
	else 
	{
		if (defined($$pep2reads{$nomod})) {$$runhash{$run}{$outfile}{'readSpt'} = $$pep2reads{$nomod};}
		else { $$runhash{$run}{$outfile}{'readSpt'} = 0; }
	}}


  $$peptidehash{$intpep}{'orig_peptide'} = $peptide;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'scan'} = $scan;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'run'} = $run;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'peptide'} = $peptide;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'tryptic'} = $tryptic;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'mod'} = $mods;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'mis'} = $mis;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'red'} = $red;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'charge'} = $charge;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'XCorr'} = $XCorr;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'Sp'} = $Sp;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'rank'} = $rank;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'dCn'} = $dCn;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'protein'} = $protein;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'ions'} = $ions;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'MH'} = $mass;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'expMH'} = $expmass;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'path'} = $fullpath;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'precursor_peak_intensity_percentage'} = $precursor_peak_intensity_percentage;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'tagLength'} = length($tag);
  #$$peptidehash{$intpep}{'outfiles'}{$outfile}{'PTM_sites'} = $PTM_sites;
  #$$peptidehash{$intpep}{'outfiles'}{$outfile}{'PTM_score'} = $PTM_score;
  $$peptidehash{$intpep}{'tryptic'} = $tryptic;
  $$peptidehash{$intpep}{'mis'} = $mis;
  $$peptidehash{$intpep}{'mod'} = $mods;
  $$peptidehash{$intpep}{'readSpt'} = $$runhash{$run}{$outfile}{'readSpt'};

    # PTM_score
    #my ($PTM_sites,$PTM_score);
    if ( defined($$outInforHash{$outkey}{PTM_sites}) and 
	defined($$outInforHash{$outkey}{PTM_score}) )
    { 
	$$peptidehash{$intpep}{'outfiles'}{$outfile}{'PTM_sites'}=$$outInforHash{$outkey}{PTM_sites};
	$$peptidehash{$intpep}{'outfiles'}{$outfile}{'PTM_score'}=$$outInforHash{$outkey}{PTM_score};
    }
###### added by yanji, add rention time, rt, and intensity
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'rt'} = $rt;
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'intensity'} = $intensity;
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'prec_int'} = $prec_int;
###### end of addition
  my $target = 0;
  for my $key (keys %proteins){
    if($key !~ /Decoy/)
    {
	$target++;
     }
	$key =~ s/\#\#//;
	if ($key =~ /\#\#/) { die "$key\n\n"; }
    $$peptidehash{$intpep}{'proteins'}{$key}++;
  }

########## count the scan number rather than protein number #####################
  if($target>0)
  {
	my $flag = 0;
       for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){#rm decoy protein names
           if($pros =~ /Decoy/)
           {
		$flag=1;
                #  $$delhash{'sharedecoy'}++;
                  delete $$peptidehash{$intpep}{'proteins'}{$pros};
            }
      }
	if($flag>0)
	{
      		$$delhash{'sharedecoy'}++;
	}
   }
  
  $$peptidehash{$intpep}{'expMH'} = $expmass;
  #%{$$runhash{$run}{$outfile}} = %{$$peptidehash{$intpep}{'outfiles'}{$outfile}};
  #$$runhash{$run}{$outfile}{'intpep'} = $intpep;

  return ($database);
}

sub pepXMLhashesTransfer_Sequest
{
  shift @_;
  my ($paramhash,$run,$outfile,$paraHash,$outInforHash,$peptidehash,$runhash,$delhash,$blank,$cpr,$cpn,$pep2reads) = @_;
  #my ($paramhash,$peptidehash,$runhash,$run,$outfile,$delhash,$blank,$cpr,$cpn,$pep2reads) = @_;

  my ($Hydrogen_mass, $HydrogenMinus_mass, $Cterm_mass)=(1.007825032,1.00335,19.017806);
  
  my %proteins;
  my $fullpath = "$run\/$outfile";

  # check if RNAseq data is vailable
  my $noRNAseq;
  if (!defined($cpr) || !defined($cpn) || !defined($pep2reads)) { $noRNAseq=1; } else { $noRNAseq=0; }
 
  #open (IN, $fullpath);
  #u170k_siRNA_02ug.14423.1.3.spout
  $outfile =~ m/\.(\d+)\.(\d+)\.(\d)\.out/;
  my ($scan,$charge) = ("$1\.$2",$3);
  my $ppi=$2; $ppi=substr($ppi,0,1);

  # check if filter ppi
  if (defined($$paramhash{ppi_included}))
  {
        my @t=split(/_/,$$paramhash{ppi_included});
        my $incl=0;
        for (my $i=0; $i<=$#t; $i++) { if ( $t[$i] == $ppi ) {$incl=1;last;} }
        unless ($incl) { $$delhash{ppi_filtered}++; return }
  }

  #my ($sorcerer) =  grep(/SageNResearch/ || /SORCERER/, <IN>);
  #if (defined($sorcerer)){ $$paramhash{'sorcererouts'} = 1; }
  my $outkey=$outfile; $outkey =~ s/(\.out)//;  #print "$outkey\n";

  #get PSM lines
  #seek(IN, 0,0);
  #my @hits=grep(/^\s+\d+\.\s+\d+\s+\//,<IN>);
  if (!defined($$outInforHash{$outkey}{hitShown}) || scalar($$outInforHash{$outkey}{hitShown})==0) { $$blank++; return; }

  #get red protain names
  #seek(IN,0,0);
  #my $red_count=0;
  #my @redProNames=grep(/^\s+\d+\s+[0-9A-Za-z#]/,<IN>);
  #close IN;

  #decide which PSM and protein name to use
  my ($mxlv,$mxtd,$mxG,$mxRpk,$dCn,$mxProName,$mxRank,$mxReadSpt)=(10,'','',0,0,'',0,0,0);
  my @tmpProteins; my $oldPep='';
  for (my $i=0; $i<$$outInforHash{$outkey}{hitShown}; $i++)
  {
    #remove redundant space
    #$hits[$i] =~ s/\s*\/\s*/\//g; $hits[$i] =~ s/(\s)\s+/$1/g; $hits[$i] =~ s/^\s+//; $hits[$i] =~ s/\s+\Z//;
    #my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $hits[$i]);
    #$dCn=$c;
    #my ($peptide, $red)  = ($g, 0);
    #if(defined($h)) { $red = $g; $red =~ s/\+//; $peptide = $h; }

    # pass values from %outInforHash
    my $k=$i+1;# corespond to %outInforHash
    my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein) = ($k,'','',$$outInforHash{$outkey}{'calc_neutral_pep_mass'}[$k]+$Hydrogen_mass,0,$$outInforHash{$outkey}{'xcorr'}[$k],'',"$$outInforHash{$outkey}{'num_matched_ions'}[$k]\/$$outInforHash{$outkey}{'tot_num_ions'}[$k]",$$outInforHash{$outkey}{'protein'}[$k]);
    if ($i==0) { $dCn=0; } else { $dCn=$$outInforHash{$outkey}{'deltacn'}[$i]; }
    my ($peptide, $red)  = ("$$outInforHash{$outkey}{'peptide_prev_aa'}[$k]\.$$outInforHash{$outkey}{'fullPeptide'}[$k]\.$$outInforHash{$outkey}{'peptide_next_aa'}[$k]",$$outInforHash{$outkey}{'num_tot_proteins'}[$k]-1);

    #last if ( $i>=1 and $dCn>0 ); # old critera

    # criteria for jumping out the loop
    #if (!defined($dCn)) { die "$outkey, $i\n"; }
    last if ( $i>=1 and $dCn>0.05 ); 
    my $pepSim=1;
    if ( $oldPep ne '' ) { $pepSim=peptide_similarity($oldPep,$peptide); }
	#if ( $pepSim<0 ) {die "uneuqla length for peptide sequences with small dCn!!!\n $outfile,$oldPep,$peptide\n";}
    last if ( $pepSim<0.9 ); 
    #last if ( $i>=10 );

    # for 1st protein name
    if ($noRNAseq) { if ($i==0) { $mxProName=$protein; } }
    else {_update_max_levels(\$mxlv,\$mxtd,\$mxG,\$mxRpk,\$mxProName,\$mxRank,$protein,$cpr,$cpn,$i);}

    # for redundant protein names
    my $j;
    for ( $j=1; $j<=$red; $j++)
    {
	#if (!defined($redProNames[$red_count])) {die "undefined red proten names lines:\n $outfile,$i,$j,$red,$red_count\n";}
      #chomp($redProNames[$red_count]); my $tmpP=$redProNames[$red_count]; $red_count++;
      ##$tmpP =~ s/^\s+\d+\s+//;  if ( $tmpP =~ /\s/ ) { $tmpP =~ s/^(.*?)\s+/$1/; }   $tmpProteins[$i][$j]=$tmpP;# old
      #if($tmpP=~/^\s+\d*\s*([\#a-zA-Z0-9\.\_\-\|\:]+)\s*/) { $tmpProteins[$i][$j]=$1;} else {die "Ero:red pro\n";}

      my $tmpP=$$outInforHash{$outkey}{'alternative_protein'}[$j][$k]; $tmpProteins[$i][$j]=$tmpP;
#  print "$outfile,$i,$k,$j,$tmpP,$tmpProteins[$i][$j]\n";

      if ($noRNAseq) {}
      else {_update_max_levels(\$mxlv,\$mxtd,\$mxG,\$mxRpk,\$mxProName,\$mxRank,$tmpP,$cpr,$cpn,$i);}
    }
    $tmpProteins[$i][0]=$j-1;

    $oldPep=$peptide;

    if ($noRNAseq) {}else {
    my $readSpt=get_read_spt($peptide,$protein,$paramhash,$pep2reads);
#if ( $outfile eq 'u170k_siRNA_02ug.62439.62439.3.out' ) { print "$outfile,$i,$mxReadSpt,$readSpt,$peptide,$protein\n"; }
    if ($mxlv==4 and $mxReadSpt<$readSpt )
    { $mxReadSpt=$readSpt; $mxRank=$i; $mxProName=$protein; }}
  }
  if ($noRNAseq) { $mxlv=1; $mxRank=0;}
  if ( $$outInforHash{$outkey}{hitShown}==1 ) { $dCn=1; }
  my $final_rank=$mxRank;
#  for (my $j=1; $j<=$tmpProteins[$final_rank][0]; $j++) {print "$outfile,$final_rank,$j,$tmpProteins[$final_rank][$j]\n";}

  #continue to use old code, but make an interface
  #open (IN, $fullpath);
  #seek(IN, 0,0);
  #mass and database
=head
  my ($mass, $database) = grep(/\+\smass\s\=/ || /proteins\s\=/, <IN>);
  my $temp=$hits[$final_rank];#interface
  $temp =~ s/\s*\/\s*/\//g; $temp =~ s/(\s)\s+/$1/g; $temp =~ s/^\s+//; $temp =~ s/\s+\Z//;
  $mass =~ s/mass\s+\=\s+([\d\.]+)\s//; $mass = $1;
  $database =~ s/\s+//g; my @temparray = split (/\,/, $database); $database = pop(@temparray);
  if($database=~/hdr/)	{$database=~s/\.hdr$//;	}
=cut
  my ($mass, $database) = ($$outInforHash{$outkey}{'calc_neutral_pep_mass'}[1]+$Hydrogen_mass+$$outInforHash{$outkey}{'massdiff'}[1],$$paraHash{'DB'});
  if($database=~/mdx/) {$database=~s/\.mdx$//;}

  # PSM line
  my $k=$final_rank+1;
  #my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $temp);
  my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e) = ($k,'','',$$outInforHash{$outkey}{'calc_neutral_pep_mass'}[$k]+$Hydrogen_mass,0,$$outInforHash{$outkey}{'xcorr'}[$k],'',"$$outInforHash{$outkey}{'num_matched_ions'}[$k]\/$$outInforHash{$outkey}{'tot_num_ions'}[$k]");
  my $protein=$mxProName;
  #print "\r$protein";

  #if (defined($sorcerer)){    $$paramhash{'sorcererouts'} = 1;  }
  
  if ($dCn < $$paramhash{'min_dCn'} || $XCorr < $$paramhash{'min_XCorr'}){
    $$delhash{'DX'}++;
    if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'DX'}++; } else { $$delhash{'target'}{'DX'}++; }
    return ($database);
  }
  #my $peptide = $g; my $red = 0; 

#my $ions = $e;
my $ions = "$$outInforHash{$outkey}{'num_matched_ions'}[$k]\/$$outInforHash{$outkey}{'tot_num_ions'}[$k]";
#if (!defined($protein)) {die "$outfile,$protein,$mxRank,$noRNAseq,$tmpProteins[0][1]\n";}
#  $protein =~ s/\,.*//;

  if (defined($$paramhash{'filter_contaminants'})){
    if ($$paramhash{'filter_contaminants'} == 1){
#      if (defined($$paramhash{'badproteins'}{$protein}) || $protein =~ /keratin/ || $protein=~ /KERATIN/){
     if (defined($$paramhash{'badproteins'}{$protein}) || $protein =~ /keratin/ || $protein=~ /CON_/){ 
       $$delhash{'contaminants'}++;
       if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'contaminants'}++; } else { $$delhash{'target'}{'contaminants'}++; }
        return ($database);
      }
    }
  }
  $proteins{$protein} = 1;
  my $target_decoy_mixed = 0;

  my ($peptide, $red)  = ("$$outInforHash{$outkey}{'peptide_prev_aa'}[$k]\.$$outInforHash{$outkey}{'fullPeptide'}[$k]\.$$outInforHash{$outkey}{'peptide_next_aa'}[$k]",$$outInforHash{$outkey}{'num_tot_proteins'}[$k]-1);
  undef($k);
  if($red>0)
  {
    #if(defined($h) && $protein !~ /Random__/ && $protein !~ /Decoy__/  && $protein !~/CON_/){
    #$red = $g; $red =~ s/\+//; $peptide = $h;
    for (my $j=1; $j<scalar(@{$tmpProteins[$final_rank]}); $j++) 
    {
      my $pro=$tmpProteins[$final_rank][$j];#print "$outfile,$final_rank,$j,$tmpProteins[$final_rank][$j]\n";

	#Remove the decoy in the redundancy??
        #mark targets in %proteins
        #if ($pro !~ /Random__/ && $pro !~ /Decoy__/ && $pro !~ /CON_/){
        if ($pro !~ /CON_/ && $pro !~ /Decoy/)
        {
                if (defined($$paramhash{'filter_contaminants'})){
                        if ($$paramhash{'filter_contaminants'} == 1){
                                if (!defined($$paramhash{'badproteins'}{$pro}) && $pro !~ /keratin/ && $pro !~ /KERATIN/ && $pro !~ /CON_/){
                                        $proteins{$pro} = 1;
                                }
                        }
                        else {
                                $proteins{$pro} = 1;
				#print "$outfile,$pro\n";
                        }
                }
        }
        if($pro =~ /Decoy/)
        {
                $target_decoy_mixed = 1;
        }
	last if ($j>=20);
    }

  }

  if($protein !~ /Decoy/ && $target_decoy_mixed == 1)
  {
	$$delhash{'sharedecoy'}++;
  }
  my $testpeptide = $peptide;
  my $mods = ($testpeptide =~ s/([\@\%\&\^\~\$\#\*\^\~]+)/$1/g) || 0;
  $testpeptide =~ s/.\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\../$1/;
  my $mis = ($testpeptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;
  if (defined($$paramhash{'max_peptide_mod'})){
    if ($mods > $$paramhash{'max_peptide_mod'}){
      #$$max_mod_del++;
      
      $$delhash{'maxmod'}++;
      if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'maxmod'}++; } else { $$delhash{'target'}{'maxmod'}++; }
      return ($database);
    }
  }

  if (defined($$paramhash{'max_peptide_mis'})){
    if ($mis > $$paramhash{'max_peptide_mis'}){
      #$$max_mis_del++;
      #print $mis."mismismis\t".$$paramhash{'max_peptide_mis'}."hashmishashmis\t".$max_mis_del."\tccccccccccc\n";
      $$delhash{'maxmis'}++;
      if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'maxmis'}++; } else { $$delhash{'target'}{'maxmis'}++; }
      return ($database);
    }
  }
  my $tryptic = $utils->istryptic($peptide);#returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
  $tryptic -= 1 if ($tryptic>1);
  my $intpep = $peptide;
  $intpep =~ s/[A-Z\-]\.([A-Z\@\%\&\^\~\$\#\*\^\~]+)\.[A-Z\-]/$1/;
  if ($$paramhash{'mix_label'} ne 0){  #???
    my $labels = "\[$$paramhash{'mix_label'}\]";
    if ($intpep =~ /$labels[\*\#\@\%\&\^\~\$\^\~]/ && ($intpep =~ /$labels[A-Z]/ || $intpep =~ /$labels\Z/)){
        $$delhash{'mixlabel'}++;
        if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'mixlabel'}++; } else { $$delhash{'target'}{'mixlabel'}++; }
				#print "$outfile $peptide $XCorr $dCn\n";
        return ($database);
    }
  }
  if (defined($$paramhash{'peptide_mod_removal'})){
    if ($$paramhash{'peptide_mod_removal'} ne 0){
      my @labels = split("", $$paramhash{'peptide_mod_removal'});
      for my $aa (@labels){
        my $label = "\[$aa\]";
        if ($intpep =~ /$label[\*\#\@\%\&\^\~\$\^\~]/){
          $$delhash{'modremoval'}++;
          if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'modremoval'}++; } else { $$delhash{'target'}{'modremoval'}++; }
          return ($database);
        }
      }
    }
  }
  if (defined($$paramhash{'peptide_aa_removal'})){
    if ($$paramhash{'peptide_aa_removal'} ne 0){
      my @labels = split("", $$paramhash{'peptide_aa_removal'});
      for my $aa (@labels){
        my $label = "\[$aa\]";
        if ($intpep =~ /$label/){
          $$delhash{'aaremoval'}++;
          if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'aaremoval'}++; } else { $$delhash{'target'}{'aaremoval'}++; }
          return ($database);
        }
      }
    }
  }
  my $nomod = $intpep; $nomod =~ s/[\@\%\&\^\~\$\#\*\^\~]//g;
  if (length($nomod) < $$paramhash{'min_peptide_length'}){
    $$delhash{'length'}++;
    if ($protein =~ /Decoy/) { $$delhash{'decoy'}{'length'}++; } else { $$delhash{'target'}{'length'}++; }
    return ($database);
  }
  my @array = split("", $nomod);
  my $calcmw = $utils->get_MW(\@array); #returns molecular weight of given sequence
  my $num = 0;
	if (defined($$paramhash{$run}{staticmods}{'nterm'})){
    $calcmw += $$paramhash{$run}{staticmods}{'nterm'};
  } elsif (defined($$paramhash{$run}{staticmods}{'cterm'})){
    $calcmw += $$paramhash{$run}{staticmods}{'cterm'};
  }
  for my $aa (keys %{$$paramhash{$run}{staticmods}}){
    $num += ($intpep =~ s/$aa/$aa/g); #print "$outfile: $peptide, $aa, $num\n";
    $calcmw += $num * $$paramhash{$run}{staticmods}->{$aa};
    $num = 0;
  }
  for my $aa (keys %{$$paramhash{$run}{dynamicmods}}){
    my $test = "\[$aa\]";
    $num += ($intpep =~ s/($test[\@\%\&\^\~\$\#\*\^\~])/$1/g);
    $calcmw += $num * $$paramhash{$run}{dynamicmods}->{$aa};
    $num = 0;
  }

#	if ($expmass != $calcmw) { die "$outfile: $expmass, $calcmw, $peptide\n"; }	
  $expmass = $calcmw;#$expmass predefined; why bothers recalculation? why diff?
  if ($$paramhash{'norm_XCorr'}){
    my $length = length($nomod);
    if ($charge == 2){
      if ($length >= 15){
        $XCorr = log($XCorr)/log(15*2);
      } else {
        $XCorr = log($XCorr)/log($length*2);
      }
    } elsif ($charge == 3){
      if ($length >= 25){
        $XCorr = log($XCorr)/log(25*4);
      } else {
        $XCorr = log($XCorr)/log($length*4);
      }
    }
  }

#print "XCorr:",$XCorr,"\n";
#print "Sp:",$Sp,"\n";
#print "rank:",$rank,"\n";
#print "dCn:",$dCn,"\n";
#print "protein:",$protein,"\n";
#print "MH:",$mass,"\n";

###### added by yanji
  # get scan number from outfile name
#  my @outfile_element = split/\./, $outfile;
#  pop(@outfile_element);
#  pop(@outfile_element);
#  my $scan_number = pop(@outfile_element);

  # get original hash, origmsms_hash
#  my $hash_dir = dirname($run)."\/".".hashes";
#  my $hash_ref = retrieve("$hash_dir/origmsms_hash");
#  my %origmsms_hash = %$hash_ref;

  # get rention time, rt, and intensity from the hash table
#  $scan_number = int($scan_number || 0);
#  my $rt = $origmsms_hash{$scan_number}{'rt'};
#  my $intensity = $origmsms_hash{$scan_number}{'intensity'};
###############added by xusheng on 06/14/2012 ###########
#  my $prec_int = $origmsms_hash{$scan_number}{'prec_int'};

  # save rt and intensity
#  $$runhash{$run}{$outfile}{'rt'} = $rt;
#  $$runhash{$run}{$outfile}{'intensity'} = $intensity;
#  $$runhash{$run}{$outfile}{'prec_int'} = $prec_int;

	
	if (defined($$peptidehash{$intpep})){ ##ADDED DMD NOVEMBER 25, 2008 to remove peptide redundancy (decoys) for peptides with real and decoy matches
		my $decoy = 0;
		my $target = 0;
		for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){
		#	print $intpep,"\t",$pros,"\n";
########## changed by xusheng #####################
#Remove all start test ^ sign from the statement
#			$decoy++ if ($pros =~ /Decoy/ || $pros =~ /Decoy__/);
#
			$target ++ if ($pros !~ /Decoy__/ && $pros !~ /Random__/);
			$decoy++ if ($pros =~ /Decoy__/ || $pros =~ /Random__/);
		}

		if ($protein !~ /Decoy__/ && $protein !~ /Random__/ && $decoy > 0){#current out is target; but previous have decoys
			for my $delout (keys %{$$peptidehash{$intpep}{'outfiles'}}){#deltete all previous outfiles of this peptide
				my $run = $$peptidehash{$intpep}{'outfiles'}{$delout}{'run'};
				delete $$runhash{$run}{$delout};
			}
			delete $$peptidehash{$intpep};#delete previous item of this peptide (for renew?)
		} elsif (($protein =~ /Decoy__/ || $protein =~ /Random__/) && $decoy == 0) {#current decoy; deltete current out
			delete $$runhash{$run}{$outfile};
			return;
		}

                for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){#no error??? delete $$peptidehash{$intpep};???
                        if($target>0 and ($pros =~ /Decoy__/ || $pros =~ /Random__/))
                        {
                                $$delhash{'sharedecoy'}++;
				delete $$peptidehash{$intpep}{'proteins'}{$pros};#rm decoy items;keep target; why not delete outfiles?
                        }
		}

	}





 	$$runhash{$run}{$outfile}{'scan'} = $scan;
        $$runhash{$run}{$outfile}{'run'} = $run;
       $$runhash{$run}{$outfile}{'peptide'} = $peptide;
        $$runhash{$run}{$outfile}{'tryptic'} = $tryptic;
        $$runhash{$run}{$outfile}{'mod'} = $mods;
        $$runhash{$run}{$outfile}{'mis'} = $mis;
        $$runhash{$run}{$outfile}{'red'} = $red;
        $$runhash{$run}{$outfile}{'charge'} = $charge;
        $$runhash{$run}{$outfile}{'XCorr'} = $XCorr;
       $$runhash{$run}{$outfile}{'Sp'} = $Sp;#??
       $$runhash{$run}{$outfile}{'rank'} = $rank;
       $$runhash{$run}{$outfile}{'dCn'} = $dCn;
       $$runhash{$run}{$outfile}{'protein'} = $protein;
       $$runhash{$run}{$outfile}{'ions'} = $ions;
       $$runhash{$run}{$outfile}{'MH'} = $mass;
       $$runhash{$run}{$outfile}{'expMH'} = $expmass;
       $$runhash{$run}{$outfile}{'path'} = $fullpath;
       $$runhash{$run}{$outfile}{'intpep'} = $intpep;
       $$runhash{$run}{$outfile}{'protein_priorty_level'} = int($mxlv);
       $$runhash{$run}{$outfile}{'rpk'} = $mxRpk;
       $$runhash{$run}{$outfile}{'psmRank'} = $final_rank;

	if ($noRNAseq) {} else {
	if ( $protein =~ /Decoy/ ) 
	{ 
		my $origPep=orig_pep_for_decoy($paramhash,$peptide);#print "$outfile,$peptide,$origPep\n";
		my @aa=split(//,$origPep); pop(@aa); shift(@aa);
		$origPep=join('',@aa);
		if (defined($$pep2reads{$origPep})) {$$runhash{$run}{$outfile}{'readSpt'} = $$pep2reads{$origPep};}
		else { $$runhash{$run}{$outfile}{'readSpt'} = 0; }
	}
	else 
	{
		if (defined($$pep2reads{$nomod})) {$$runhash{$run}{$outfile}{'readSpt'} = $$pep2reads{$nomod};}
		else { $$runhash{$run}{$outfile}{'readSpt'} = 0; }
	}}


  $$peptidehash{$intpep}{'orig_peptide'} = $peptide;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'scan'} = $scan;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'run'} = $run;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'peptide'} = $peptide;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'tryptic'} = $tryptic;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'mod'} = $mods;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'mis'} = $mis;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'red'} = $red;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'charge'} = $charge;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'XCorr'} = $XCorr;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'Sp'} = $Sp;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'rank'} = $rank;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'dCn'} = $dCn;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'protein'} = $protein;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'ions'} = $ions;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'MH'} = $mass;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'expMH'} = $expmass;
  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'path'} = $fullpath;
  $$peptidehash{$intpep}{'tryptic'} = $tryptic;
  $$peptidehash{$intpep}{'mis'} = $mis;
  $$peptidehash{$intpep}{'mod'} = $mods;
  $$peptidehash{$intpep}{'readSpt'} = $$runhash{$run}{$outfile}{'readSpt'};

###### added by yanji, add rention time, rt, and intensity
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'rt'} = $rt;
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'intensity'} = $intensity;
#  $$peptidehash{$intpep}{'outfiles'}{$outfile}{'prec_int'} = $prec_int;
###### end of addition
  my $target = 0;
  for my $key (keys %proteins){
    if($key !~ /Decoy/)
    {
	$target++;
     }
	$key =~ s/\#\#//;
    $$peptidehash{$intpep}{'proteins'}{$key}++;
  }

########## count the scan number rather than protein number #####################
  if($target>0)
  {
	my $flag = 0;
       for my $pros (keys %{$$peptidehash{$intpep}{'proteins'}}){#rm decoy protein names
           if($pros =~ /Decoy/)
           {
		$flag=1;
                #  $$delhash{'sharedecoy'}++;
                  delete $$peptidehash{$intpep}{'proteins'}{$pros};
            }
      }
	if($flag>0)
	{
      		$$delhash{'sharedecoy'}++;
	}
   }
  
  $$peptidehash{$intpep}{'expMH'} = $expmass;
  #%{$$runhash{$run}{$outfile}} = %{$$peptidehash{$intpep}{'outfiles'}{$outfile}};
  #$$runhash{$run}{$outfile}{'intpep'} = $intpep;

  return ($database);
}

sub printAcceptedPSMpepXML
{
	shift @_;
	my ($outInforHash,$paraHash,$peptidehash,$xml)=@_;

	# store accepted outfiles in a hash
	my %acpOut;
	while ( my ($pep, $hash) = each %$peptidehash )
	{
		foreach my $out (keys %{$$hash{'outfiles'}})
		{
			#print "  Parsing $out\r";
			$out =~ s/\.out$//;
			$out =~ s/\.spout$//;
			$acpOut{$out}='';
		}
	}

	# filter %outInforHash using %acpOut
	foreach my $out (keys %$outInforHash)
	{
		if (!defined($acpOut{$out}))
		{
			delete $$outInforHash{$out};
		}
	}

	# print pepXML using pepXML_parser module
	my $pepParser=pepXML_parser->new();
	$pepParser->printPepXML($paraHash,$outInforHash,$xml,1);
}

