#!/usr/bin/perl -wT
package pepXML_parser;
use strict;
use vars qw($Hydrogen_mass $HydrogenMinus_mass $Cterm_mass);

BEGIN
{
	($Hydrogen_mass, $HydrogenMinus_mass, $Cterm_mass)=(1.007825032,1.00335,19.017806);
}

sub new
{
	my($class) = @_;
	my $self = {};
	bless ($self,$class);
	$self->getNowDate();

	return $self;
}

#---------------------------------------------------------------------------------------------------------
1;

sub printPSMinfor
{
	shift @_;
	my ($hash,$output)=@_;

	open(OUT,">$output");

	print OUT "outfile\tscore\tdeltaScore\tprecursor_neutral_mass\tcalc_neutral_pep_mass\tpeptide\tprotein\tTD\n";
	foreach my $out (keys %$hash)
	{
		my ($score,$deltaScore);
		if (defined($$hash{$out}{WeightedEvalue}[1])) { $score=$$hash{$out}{WeightedEvalue}[1]; $deltaScore=$$hash{$out}{deltacn}[1]; }
		elsif (defined($$hash{$out}{xcorr}[1])) { $score=$$hash{$out}{xcorr}[1]; $deltaScore=$$hash{$out}{deltacn}[1]; }
		else {next;}

		my $precursor_neutral_mass=$$hash{$out}{precursor_neutral_mass};
		my $calc_neutral_pep_mass=$$hash{$out}{calc_neutral_pep_mass}[1];
		my $peptide=$$hash{$out}{peptide}[1];
		my $protein=$$hash{$out}{protein}[1];
		my $TD; if ($protein =~ m/Decoy/) {$TD='decoy'; $protein =~ s/\#\#//;} else {$TD='target';}

		print OUT "$out\t$score\t$deltaScore\t$precursor_neutral_mass\t$calc_neutral_pep_mass\t$peptide\t$protein\t$TD\n";
	}

	close OUT;
}

sub initialize_outfile
{
	my ($oldhash,$newhash,$out)=@_;

	my $hash1=$$oldhash{$out};
	$$newhash{$out}{'scan'}=$$hash1{'scan'};
	$$newhash{$out}{'precursor_neutral_mass'}=$$hash1{'precursor_neutral_mass'};
	$$newhash{$out}{'precursor_peak_intensity_order'}=$$hash1{'precursor_peak_intensity_order'};
	$$newhash{$out}{'charge'}=$$hash1{'charge'};
	$$newhash{$out}{'hitShown'}=0;
}

sub replace_seq_ID
{
	shift @_;
	my ($hash,$idhash)=@_;

	for ( my $i=1; $i<=$$hash{'hitShown'}; $i++ )
	{
		my $pro=$$hash{'protein'}[$i];
		if ( $pro =~ /Decoy/ ) { $pro =~ s/\#\#Decoy__//; }
		if (defined($$idhash{$pro})) {$$hash{'protein'}[$i] =~ s/$pro/$$idhash{$pro}/; }

		for (my $j=1; $j<=$$hash{'alternative_protein'}[0][$i]; $j++)
		{
			$pro=$$hash{'alternative_protein'}[$j][$i];
			if ( $pro =~ /Decoy/ ) { $pro =~ s/\#\#Decoy__//; }
			if (defined($$idhash{$pro})) {$$hash{'alternative_protein'}[$j][$i] =~ s/$pro/$$idhash{$pro}/; }
		}
	}
}

sub replace_one_seq_ID
{
	my ($protein,$idhash)=@_;
	my $pro=$protein;
	if ( $pro =~ /Decoy/ ) { $pro =~ s/\#\#Decoy__//; }
}

sub add_filtered_outfile
{
	shift @_;
	my ($paraHash,$oldhash,$newhash,$out,$pep2read,$threshold)=@_;

	initialize_outfile($oldhash,$newhash,$out);

	my $hash1=$$oldhash{$out};
	my ($hash3)=($$newhash{$out});

	#if (!defined($$hash1{'hitShown'})) {die "$out,$$hash1{'precursor_neutral_mass'}\n";}
	for ( my $i=1; $i<=$$hash1{'hitShown'}; $i++ )
	{
		my $pep=$$hash1{'peptide'}[$i];
		my $pro=$$hash1{'protein'}[$i];

		if ( $pro =~ /Decoy/ ) { $pep=orig_pep_for_decoy($pep,'R2'); }
		if ( defined($$pep2read{$pep}) && $$pep2read{$pep}>=$threshold )
		{
			add_psm($paraHash,$hash1,$hash3,$i);
		}
	}
	
	# recalculate dCn
	calculate_dCn($hash3);
}

sub orig_pep_for_decoy
{
	my ($peptide,$decoy_strategy)=@_;

	if ( $decoy_strategy eq 'R2' )
        {
                $peptide =~ s/[\@\#\*\^\~\$\.]//g;
                my @aa=split(//,$peptide);
                my $lst=pop(@aa);
                unshift(@aa,$lst);

                $peptide=join('',@aa);
                $peptide=reverse($peptide);
                return $peptide;
        }
 	else { die "undefined decoy strategy: $decoy_strategy!!!\n"; }
}

sub add_outfile
{
	shift @_;
	my ($paraHash,$oldhash,$newhash,$out)=@_;

	my $hash1=$$oldhash{$out};
	$$newhash{$out}{'scan'}=$$hash1{'scan'};
	$$newhash{$out}{'precursor_neutral_mass'}=$$hash1{'precursor_neutral_mass'};
	$$newhash{$out}{'precursor_peak_intensity_order'}=$$hash1{'precursor_peak_intensity_order'};
	$$newhash{$out}{'charge'}=$$hash1{'charge'};
	$$newhash{$out}{'hitShown'}=0;
	my ($hash3)=($$newhash{$out});

	for ( my $i=1; $i<=$$hash1{'hitShown'}; $i++ ) { add_psm($paraHash,$hash1,$hash3,$i); }
	
	# recalculate dCn
	calculate_dCn($hash3);
}

sub merge_outfiles
{
	shift @_;
	my ($paraHash,$outhash1,$outhash2,$outhash3,$out)=@_;

	my ($hash1,$hash2)=($$outhash1{$out},$$outhash2{$out});

	$$outhash3{$out}{'scan'}=$$hash1{'scan'};
	$$outhash3{$out}{'precursor_neutral_mass'}=$$hash1{'precursor_neutral_mass'};
	$$outhash3{$out}{'precursor_peak_intensity_order'}=$$hash1{'precursor_peak_intensity_order'};
	$$outhash3{$out}{'charge'}=$$hash1{'charge'};
	$$outhash3{$out}{'hitShown'}=0;
	my ($hash3)=($$outhash3{$out});

	# an array to mark whether a psm in hash2 is used
	my @array;
	for ( my $i=1; $i<=$$hash2{'hitShown'}; $i++ ) { $array[$i]=0; }

	# merge hash1 and hash2
	my ($hitnum1,$hitnum2)=(1,1);

	my $scoreType;
	if ( defined($$hash1{'xcorr'}[$hitnum1]) or defined($$hash1{'xcorr'}[$hitnum2]) )
	{ $scoreType='xcorr'; }
	else { $scoreType='WeightedEvalue'; }

	while ( $hitnum1<=$$hash1{'hitShown'} and $hitnum2<=$$hash2{'hitShown'} )
	{
		if ( $array[$hitnum2] ) { $hitnum2++; next; }  # check if this PSM in hash2 is used
		if ( $$hash1{'fullPeptide'}[$hitnum1] eq $$hash2{'fullPeptide'}[$hitnum2] )
		{
			if ( $$hash1{$scoreType}[$hitnum1] >= $$hash2{$scoreType}[$hitnum2] ) # equal scores or peptide in hash1 has a higher score
			{
				#die "Same peptide, not equal score for the same scan: $out,$hitnum1,$hitnum2\n$$hash1{'fullPeptide'}[$hitnum1],$$hash1{$scoreType}[$hitnum1],$$hash2{'fullPeptide'}[$hitnum2],$$hash2{$scoreType}[$hitnum2]\n";
				#print "$out,$hitnum1,$hitnum2: $$hash1{'fullPeptide'}[$hitnum1],$$hash1{$scoreType}[$hitnum1],$$hash2{'fullPeptide'}[$hitnum2],$$hash2{$scoreType}[$hitnum2]\n";
				#merge_psm($hash1,$hash2,$hash3,$hitnum1,$hitnum2); $hitnum1++; $hitnum2++;
				add_psm($paraHash,$hash1,$hash3,$hitnum1);
				add_redundant_protein($hash2,$hash3,$hitnum2);
			}
			else	# peptide in hash2 has a higher score
			{
				add_psm($paraHash,$hash2,$hash3,$hitnum2);
				add_redundant_protein($hash1,$hash3,$hitnum1);
			}
			$hitnum1++; $array[$hitnum2]=1; $hitnum2++;
		}
		else
		{
			if ( $$hash1{$scoreType}[$hitnum1] == $$hash2{$scoreType}[$hitnum2] ) 
			{
				add_psm($paraHash,$hash1,$hash3,$hitnum1);
				# search potential I/L ambiguity peptides in hash2
				my $i=$hitnum2+1;
				while ( $i<=$$hash2{'hitShown'} && $$hash2{$scoreType}[$i]==$$hash1{$scoreType}[$hitnum1] )
				{
					if ( $$hash2{'fullPeptide'}[$i] eq $$hash1{'fullPeptide'}[$hitnum1])
					{
						add_redundant_protein($hash2,$hash3,$i);  $array[$i]=1;
						last;
					}
					$i++;
				}
				$hitnum1++;
			}
			elsif ( $$hash1{$scoreType}[$hitnum1] > $$hash2{$scoreType}[$hitnum2] )
			{
				add_psm($paraHash,$hash1,$hash3,$hitnum1); $hitnum1++;
			}
			elsif ( $$hash1{$scoreType}[$hitnum1] < $$hash2{$scoreType}[$hitnum2] )
			{
				add_psm($paraHash,$hash2,$hash3,$hitnum2); $array[$hitnum2]=1;$hitnum2++;
			}
			
		}
	}

	# add the rest of hash1
	for (my $i=$hitnum1; $i<=$$hash1{'hitShown'}; $i++)
	{
		add_psm($paraHash,$hash1,$hash3,$i); 
	}

	# add the rest of hash2
	for (my $i=$hitnum2; $i<=$$hash2{'hitShown'}; $i++)
	{
		add_psm($paraHash,$hash2,$hash3,$i);  $array[$i]=1;
	}

	# recalculate dCn
	calculate_dCn($hash3);

	# delete items
	delete($$outhash1{$out});
	delete($$outhash2{$out});
}

sub calculate_dCn
{
	my ($hash)=@_;

	my $hitnum=$$hash{'hitShown'};

	my $scoreType;
	if ( defined($$hash{'xcorr'}[1]) ) { $scoreType='xcorr'; }
	else { $scoreType='WeightedEvalue'; }

	for (my $i=1; $i<$$hash{'hitShown'}; $i++)
	{
		#$$hash{'deltacn'}[$i]=($$hash{'xcorr'}[$i]-$$hash{'xcorr'}[$i+1])/$$hash{'xcorr'}[$i]; # compare adjacent psm
		$$hash{'deltacn'}[$i]=($$hash{$scoreType}[1]-$$hash{$scoreType}[$i+1])/$$hash{$scoreType}[1]; # sequest method: always compare with the top hit
	}
	$$hash{'deltacn'}[$$hash{'hitShown'}]=1;# last psm
}

sub add_redundant_protein
{
	my ($oldhash,$newhash,$hitnum)=@_;

	my $newhit=$$newhash{'hitShown'};

	# store the existed protein names
	my %existedPro;
	for ( my $i=1; $i<=$$newhash{'alternative_protein'}[0][$newhit]; $i++ ) 
	{ 
		$existedPro{$$newhash{'alternative_protein'}[$i][$newhit]}=''; 
	}

	# main protein in hash2 psm
	if (!defined($existedPro{$$oldhash{'protein'}[$hitnum]}))
	{
		$$newhash{'alternative_protein'}[0][$newhit]++;
		$$newhash{'alternative_protein'}[$$newhash{'alternative_protein'}[0][$newhit]][$newhit]=$$oldhash{'protein'}[$hitnum];
	}

	# alternative proteins
	for ( my $i=1; $i<=$$oldhash{'alternative_protein'}[0][$hitnum]; $i++ )
	{
		if (!defined($existedPro{$$oldhash{'alternative_protein'}[$i][$hitnum]}))
		{
			$$newhash{'alternative_protein'}[0][$newhit]++;
			$$newhash{'alternative_protein'}[$$newhash{'alternative_protein'}[0][$newhit]][$newhit]=$$oldhash{'alternative_protein'}[$i][$hitnum];
		}
	}
	#$$newhash{'alternative_protein'}[0][$newhit]+=$$oldhash{'alternative_protein'}[0][$hitnum];
}

sub add_psm
{
	my ($paraHash,$oldhash,$newhash,$hitnum)=@_;

	$$newhash{'hitShown'}++; my $newhit=$$newhash{'hitShown'};

	# required fields
	foreach my $k (keys %{$$paraHash{hitRequiredParams}})
	{ $$newhash{$k}[$newhit]=$$oldhash{$k}[$hitnum]; }

	# optional fields
	foreach my $k (keys %{$$paraHash{hitOptionalParams}})
	{ $$newhash{$k}[$newhit]=$$oldhash{$k}[$hitnum]; }

=head
	# required and optional fields
	foreach my $k (keys %$oldhash)
	{
		if ( defined($$oldhash{$k}[$hitnum]) ) { $$newhash{$k}[$newhit]=$$oldhash{$k}[$hitnum]; }
		#next if ($k =~ /scan/ || $k =~ /precursor_neutral_mass/  || $k =~ /charge/ || $k =~ /hitShown/ || $k =~ /alternative_protein/ || $k =~ /dynamicMod/ || $k =~ /staticMod/);
		#$$newhash{$k}[$newhit]=$$oldhash{$k}[$hitnum];
	}
=cut
	# alternative_protein fields
	for ( my $i=0; $i<=$$oldhash{'alternative_protein'}[0][$hitnum]; $i++ )
	{
		$$newhash{'alternative_protein'}[$i][$newhit]=$$oldhash{'alternative_protein'}[$i][$hitnum];
	}

	# modification fields
	if (defined($$oldhash{'dynamicMod'}{'found'}[$hitnum])) {
	foreach my $k (keys %{$$oldhash{'dynamicMod'}})
	{
		if ($k =~ /found/ ) {$$newhash{'dynamicMod'}{$k}[$newhit]=$$oldhash{'dynamicMod'}{$k}[$hitnum];}
		else
		{
			for (my $j=1; $j<=$$oldhash{'dynamicMod'}{'found'}[$hitnum]; $j++) 
			{
				$$newhash{'dynamicMod'}{$k}[$j][$newhit]=$$oldhash{'dynamicMod'}{$k}[$j][$hitnum];
			}
		}
	}}
	if (defined($$oldhash{'staticMod'}{'found'}[$hitnum])) {
	foreach my $k (keys %{$$oldhash{'staticMod'}})
	{
		if ($k =~ /found/ ) {$$newhash{'staticMod'}{$k}[$newhit]=$$oldhash{'staticMod'}{$k}[$hitnum];}
		else
		{
			for (my $j=1; $j<=$$oldhash{'staticMod'}{'found'}[$hitnum]; $j++) 
			{
				$$newhash{'staticMod'}{$k}[$j][$newhit]=$$oldhash{'staticMod'}{$k}[$j][$hitnum];
			}
		}
	}}
}

sub gettime{
   my ($sec, $min, $hour, $day, $mon, $year) = localtime(time);
   my ($now_date, $now_time);
   $now_date=join("-",($year+1900,$mon+1,$day));
   $now_time=join(":",($hour,$min,$sec));
   return ($now_date, $now_time);
}

sub getNowDate
{
        my $self = shift;
        my ($now_date, $now_time) = &gettime();
        $self->{'now_date'}=$now_date;
        $self->{'now_time'}=$now_time;
}


sub enzymeStat{
my ($AA, $status)=@_;
my (%enzymeStat);
$enzymeStat{'trypsin'}{'cut'}='KR';
$enzymeStat{'trypsin'}{'nocuts'}='P';
$enzymeStat{'trypsin'}{'sense'}='C';

$enzymeStat{'tryptic'}{'cut'}='KR';
$enzymeStat{'tryptic'}{'nocuts'}='P';
$enzymeStat{'tryptic'}{'sense'}='C';

$enzymeStat{'argc'}{'cut'}='R';
$enzymeStat{'argc'}{'nocuts'}='P';
$enzymeStat{'argc'}{'sense'}='C';

$enzymeStat{'aspn'}{'cut'}='D';
#$enzymeStat{'aspn'}{'nocuts'}='';
$enzymeStat{'aspn'}{'sense'}='N';

$enzymeStat{'chymotrypsin'}{'cut'}='YWFM';
$enzymeStat{'chymotrypsin'}{'nocuts'}='P';
$enzymeStat{'chymotrypsin'}{'sense'}='C';

$enzymeStat{'clostripain'}{'cut'}='R';
$enzymeStat{'clostripain'}{'nocuts'}='-';
$enzymeStat{'clostripain'}{'sense'}='C';

$enzymeStat{'cnbr'}{'cut'}='M';
$enzymeStat{'cnbr'}{'nocuts'}='P';
$enzymeStat{'cnbr'}{'sense'}='C';

$enzymeStat{'elastase'}{'cut'}='GVLIA';
$enzymeStat{'elastase'}{'nocuts'}='P';
$enzymeStat{'elastase'}{'sense'}='C';

if (exists($enzymeStat{$AA}{$status})){ return $enzymeStat{$AA}{$status};}
else { return ''; }

}

sub printPepXML
{
        my ($self,$seqPara, $outInfor, $outXML, $topHit)=@_;
	my $run='';

        open(XML,">$outXML.pepXML");

        print XML '<?xml version="1.0" encoding="Accoding to Out2XML(TPP v2.9 GALE rev.3, Build 200611091255)"?>',"\n";
        print XML '<msms_pipeline_analysis date="',$self->{'now_date'},' ', $self->{'now_time'},'" ','summary_xml="',$outXML,'.pepXML">',"\n";
        print XML '<msms_run_summary base_name="',$run,'" raw_data_type="raw" raw_data="">',"\n";

        print XML '<sample_enzyme name="',$$seqPara{'enzyme'},'">',"\n";
        print XML '<specificity cut="',enzymeStat($$seqPara{'enzyme'},'cut');
        unless (enzymeStat($$seqPara{'enzyme'},'nocuts') eq '')
        {
                print XML '" no_cut="',enzymeStat($$seqPara{'enzyme'},'nocuts');
        }
        print XML '" sense="',enzymeStat($$seqPara{'enzyme'},'sense'),'"/>',"\n";
        print XML '</sample_enzyme>',"\n";

        print XML '<search_summary base_name="',$run,'" search_engine="',$$seqPara{'search_engine'},'" precursor_mass_type="',$$seqPara{precursor_mass_type},'" fragment_mass_type="',$$seqPara{fragment_mass_type},'" out_data_type="',$$seqPara{out_data_type},'" out_data="',$$seqPara{out_data},'" search_id="">',"\n";
        print XML '<search_database local_path="',$$seqPara{'DB'},'" type="AA"/>',"\n";
        if (defined($$seqPara{'dynamicMod'}{'count'}))
        {
                for (my $i=1; $i<=$$seqPara{'dynamicMod'}{'count'}; $i++)
                {
                        print XML '<aminoacid_modification aminoacid="',$$seqPara{'dynamicMod'}{'AA'}[$i],'" massdiff="',$$seqPara{'dynamicMod'}{'massdiff'}[$i],'" mass="',$$seqPara{'dynamicMod'}{'mass'}[$i],'" variable="Y" symbol="',$$seqPara{'dynamicMod'}{'symble'}[$i],'"/>',"\n";
                }
        }
        if (defined($$seqPara{'staticMod'}{'count'}))
        {
                for (my $i=1; $i<=$$seqPara{'staticMod'}{'count'}; $i++)
                {
                        print XML '<aminoacid_modification aminoacid="',$$seqPara{'staticMod'}{'AA'}[$i],'" massdiff="',$$seqPara{'staticMod'}{'massdiff'}[$i],'" mass="',$$seqPara{'staticMod'}{'mass'}[$i],'" variable="N"/>',"\n";
                }
        }
        if (defined($$seqPara{'NtermMod'}{'mass'}))
        {
                print XML '<terminal_modification terminus="n" massdiff="',$$seqPara{'NtermMod'}{'massdiff'},'" mass="',$$seqPara{'NtermMod'}{'mass'},'" variable="N" protein_terminus="N"/>',"\n";
        }
        if (defined($$seqPara{'CtermMod'}{'mass'}))
        {
                print XML '<terminal_modification terminus="c" massdiff="',$$seqPara{'CtermMod'}{'massdiff'},'" mass="',$$seqPara{'CtermMod'}{'mass'},'" variable="N" protein_terminus="C"/>',"\n";
        }

        #OPTIONAL PARAMETERS
        foreach my $tmpTitle (keys %{$$seqPara{'optionalParams'}})
        {#print "$tmpTitle,$$seqPara{$tmpTitle}\n";
                if (defined($$seqPara{$tmpTitle}))
                {
                        print XML "\<parameter name=\"$tmpTitle\" value=\"$$seqPara{$tmpTitle}\"\/\>\n";
                }
        }

        print XML '</search_summary>',"\n";

        #opendir(DIR,$run); my @out = grep {/\.out\Z/} readdir(DIR);
        #closedir(DIR);
        my @out = (keys %$outInfor);

        my %tmpoutfiles;
        for my $outfile (@out)
        {
                #Rat_B_100ng_Q.10000.10000.3.spout
                #print "$outfile\n";
                $outfile =~ /\.(\d+)\.(\d+)\.(\d+)/;
                $tmpoutfiles{$outfile}=$1;
        }
        @out=sort {$tmpoutfiles{$a} <=> $tmpoutfiles{$b}} keys %tmpoutfiles;

        #@out=sort(@out);
        my $count=0;
        for my $outfile (@out)
        {
                $count++;
                #$outfile =~ s/(.*?)(\.out)/$1$2/;  print "$outfile\n";
                print XML "\<spectrum_query spectrum=\"$outfile\" start_scan=\"$$outInfor{$outfile}{'scan'}\" end_scan=\"$$outInfor{$outfile}{'scan'}\" precursor_neutral_mass=\"",$$outInfor{$outfile}{'precursor_neutral_mass'},"\" ";
		if (defined($$outInfor{$outfile}{'precursor_peak_intensity_percentage'})) 
		{ print XML "precursor_peak_intensity_percentage=\"$$outInfor{$outfile}{'precursor_peak_intensity_percentage'}\" "; }
		if (defined($$outInfor{$outfile}{'precursor_peak_intensity_order'})) 
		{ print XML "precursor_peak_intensity_order=\"$$outInfor{$outfile}{'precursor_peak_intensity_order'}\" "; }
		if (defined($$outInfor{$outfile}{'precursor_matched_peptides'})) 
		{ print XML "precursor_matched_peptides=\"$$outInfor{$outfile}{'precursor_matched_peptides'}\" "; }
		print XML "assumed_charge=\"$$outInfor{$outfile}{'charge'}\" index=\"$count\"\>\n";
                print XML '<search_result>',"\n";

                #print "$$outInfor{$outfile}{'hitShown'}\n";
                for (my $i=1; $i<=$$outInfor{$outfile}{'hitShown'} and $i<=$topHit; $i++ )
                {

                        print XML "\<search_hit hit_rank=\"",$i, '" peptide="',$$outInfor{$outfile}{'peptide'}[$i],'" peptide_prev_aa="',$$outInfor{$outfile}{'peptide_prev_aa'}[$i],'" peptide_next_aa="',$$outInfor{$outfile}{'peptide_next_aa'}[$i],
                        '" protein="',$$outInfor{$outfile}{'protein'}[$i],'" num_tot_proteins="',$$outInfor{$outfile}{'num_tot_proteins'}[$i],'" num_matched_ions="',$$outInfor{$outfile}{'num_matched_ions'}[$i],'" tot_num_ions="',
                        $$outInfor{$outfile}{'tot_num_ions'}[$i],'" calc_neutral_pep_mass="',$$outInfor{$outfile}{'calc_neutral_pep_mass'}[$i],'" massdiff="',$$outInfor{$outfile}{'massdiff'}[$i],
                        '" num_tol_term="',$$outInfor{$outfile}{'num_tol_term'}[$i],'" num_missed_cleavages="',$$outInfor{$outfile}{'num_missed_cleavages'}[$i],'" is_rejected="',0,"\"\>\n";

                        #alternative_protein
                        for (my $j=1; $j<=$$outInfor{$outfile}{'alternative_protein'}[0][$i]; $j++)
                        {
                                unless (defined($$outInfor{$outfile}{'alternative_protein'}[$j][$i])) {print "$outfile,red=$$outInfor{$outfile}{'red'}[$i],i=$i,j=$j\n";next;}
                                print XML '<alternative_protein protein="',$$outInfor{$outfile}{'alternative_protein'}[$j][$i],'"/>',"\n";
                        }

                        #modification
                        if (defined($$seqPara{'modification'}{'Nterm'}{'mass'}) or defined($$seqPara{'modification'}{'Cterm'}{'mass'}) or
                        (defined($$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]) and $$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]) or
                        (defined($$outInfor{$outfile}{'staticMod'}{'found'}[$i]) and $$outInfor{$outfile}{'staticMod'}{'found'}[$i]) )
                        {
                                print XML '<modification_info';
                                if (defined($$seqPara{'modification'}{'Nterm'}{'mass'}))
                                {
                                        print XML  ' mod_nterm_mass="',$$seqPara{'modification'}{'Nterm'}{'mass'},"\"";
                                }
                                elsif (defined($$seqPara{'modification'}{'Cterm'}{'mass'}))
                                {
                                        print XML  ' mod_cterm_mass="',$$seqPara{'modification'}{'Cterm'}{'mass'},"\"";
                                }
                                print XML "\>\n";

                                if (defined($$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]))
                                {#print "\nmark,$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]";
                                        for (my $j=1; $j<=$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]; $j++)
                                        {
                                                print XML '<mod_aminoacid_mass position="',$$outInfor{$outfile}{'dynamicMod'}{'position'}[$j][$i],'" mass="',$$outInfor{$outfile}{'dynamicMod'}{'mass'}[$j][$i],'"/>',"\n";
                                        }
                                }

                                if (defined($$outInfor{$outfile}{'staticMod'}{'found'}[$i]))
                                {#print "\nmark,$$outInfor{$outfile}{'staticMod'}{'found'}[$i]";
                                        for (my $j=1; $j<=$$outInfor{$outfile}{'staticMod'}{'found'}[$i]; $j++)
                                        {
                                                print XML '<mod_aminoacid_mass position="',$$outInfor{$outfile}{'staticMod'}{'position'}[$j][$i],'" mass="',$$outInfor{$outfile}{'staticMod'}{'mass'}[$j][$i],'"/>',"\n";
                                        }
                                }



                                print XML '</modification_info>',"\n";
                        }

                        #print XML '<search_score name="xcorr" value="',$$outInfor{$outfile}{'XCorr'}[$i],'"/>',"\n";
                        #print XML '<search_score name="deltacn" value="',$$outInfor{$outfile}{'deltCn'}[$i],'"/>',"\n";
                        #print XML '<search_score name="deltacnstar" value="',0,'"/>',"\n";
                        #print XML '<search_score name="spscore" value="',$$outInfor{$outfile}{'Sp'}[$i],'"/>',"\n";
                        #print XML '<search_score name="sprank" value="',$$outInfor{$outfile}{'spRank'}[$i],'"/>',"\n";

                        foreach my $tmpTitle (keys %{$$seqPara{'hitOptionalParams'}})
                        {
                                if (defined($$outInfor{$outfile}{$tmpTitle}[$i]))
                                {
                                        print XML "\<search_score name=\"$tmpTitle\" value=\"$$outInfor{$outfile}{$tmpTitle}[$i]\"\/\>\n";
                                }
                        }

                        print XML '</search_hit>',"\n";
                }

                print XML '</search_result>',"\n";
                print XML '</spectrum_query>',"\n";
        }

        print XML '</msms_run_summary>',"\n";
        print XML '</msms_pipeline_analysis>',"\n";

        close(XML);

}



sub pepXML2Hashes
{
	shift @_;
	my ($paraHash, $outInforHash, $xml)=@_;
	
	open (IN, $xml);  

	# read header lines of pepXML file, initilize %paraHash
	while (<IN>)
	{
		#required fields
		if (/sample_enzyme name=\"(\w+)\"\>\Z/) {$$paraHash{'enzyme'}=lc($1); next;}
		if (/search_database local_path=\"(.*?)\" /) {$$paraHash{'DB'}=$1; next;}
		if (/search_engine=\"(.*?)\" precursor_mass_type=\"(.*?)\" fragment_mass_type=\"(.*?)\" out_data_type=\"(.*?)\" out_data=\"(.*?)\"/)
		{
			$$paraHash{'search_engine'}=$1;
			$$paraHash{'precursor_mass_type'}=$2;
			$$paraHash{'fragment_mass_type'}=$3;
			$$paraHash{'out_data_type'}=$4;
			$$paraHash{'out_data'}=$5;
			next;
		}
		
		#modifications
		if (/aminoacid_modification aminoacid\=\"(\w)\" massdiff=\"([\d\.]+)\" mass=\"([\d\.]+)\" variable=\"Y\" symbol=\"(.*)\"\/\>/)
		{
			unless (defined($$paraHash{'dynamicMod'}{'count'})) {$$paraHash{'dynamicMod'}{'count'}=0;}
			
			$$paraHash{'dynamicMod'}{'count'}++;
			$$paraHash{'dynamicMod'}{'AA'}[$$paraHash{'dynamicMod'}{'count'}]=$1;
			$$paraHash{'dynamicMod'}{'massdiff'}[$$paraHash{'dynamicMod'}{'count'}]=$2;
			$$paraHash{'dynamicMod'}{'mass'}[$$paraHash{'dynamicMod'}{'count'}]=$3;
			$$paraHash{'dynamicMod'}{'symble'}[$$paraHash{'dynamicMod'}{'count'}]=$4;
			next;
		}
		if (/aminoacid_modification aminoacid\=\"(\w)\" massdiff=\"([\d\.]+)\" mass=\"([\d\.]+)\" variable=\"N\"\/\>/)
		{
			unless (defined($$paraHash{'staticMod'}{'count'})) {$$paraHash{'staticMod'}{'count'}=0;}
			
			$$paraHash{'staticMod'}{'count'}++;
			$$paraHash{'staticMod'}{'AA'}[$$paraHash{'staticMod'}{'count'}]=$1;
			$$paraHash{'staticMod'}{'massdiff'}[$$paraHash{'staticMod'}{'count'}]=$2;
			$$paraHash{'staticMod'}{'mass'}[$$paraHash{'staticMod'}{'count'}]=$3;
			next;
		}
		if (/terminal_modification terminus\=\"[Nn]\" massdiff=\"([\d\.]+)\" mass=\"([\d\.]+)\" variable=\"N\" protein_terminus=\"[Nn]\"\/\>/)
		{			
			#$$paraHash{'NtermMod'}=1;
			$$paraHash{'NtermMod'}{'massdiff'}=$1;
			$$paraHash{'NtermMod'}{'mass'}=$2;
			next;
		}
		if (/terminal_modification terminus\=\"[Cc]\" massdiff=\"([\d\.]+)\" mass=\"([\d\.]+)\" variable=\"N\" protein_terminus=\"[Cc]\"\/\>/)
		{			
			#$$paraHash{'NtermMod'}=1;
			$$paraHash{'CtermMod'}{'massdiff'}=$1;
			$$paraHash{'CtermMod'}{'mass'}=$2;
			next;
		}
		
		#optional parameters
		if (/\<parameter name\=\"(.*?)\" value\=\"(.*?)\"\/\>/) 
		{
			$$paraHash{$1}=$2; 
			$$paraHash{'optionalParams'}{$1}='';
			next;
		}
	
	
	
		last if (/^\<\/search_summary\>/);
	}
	#my $line=<IN>; 	print $$paraHash{'ion_series'};
	
	# read through pepXML, and construct %outInforHash
	my ($outfile, $altProCount);
	while (<IN>)
	{
		if (/^\<spectrum_query /)
		{
			my @tmp = / (\w+)\=\"(.*?)\"/g;
			$outfile=$tmp[1];
			for (my $i=0; $i<scalar(@tmp); $i=$i+2)
			{
				$$outInforHash{$outfile}{$tmp[$i]}=$tmp[$i+1];
			}

			# scan
			$outfile =~ m/\.(\d+)\.(\d+)\.(\d)$/;
			$$outInforHash{$outfile}{'scan'}=$1;
			$$outInforHash{$outfile}{'charge'}=$3;
			
			$$outInforHash{$outfile}{'hitShown'}=0;
			next;
		}
=head
		if (/^\<spectrum_query spectrum=\"(.*?)\" start_scan=\"(\d+)\" end_scan=\"\d+\" precursor_neutral_mass=\"([\d\.]+)\" precursor_peak_intensity_percentage=\"[\d\.\%]+\" precursor_peak_intensity_order=\"\d+\" assumed_charge=\"(\d+)\" index=\"\d+\"\>/)
		{
			$outfile=$1;
			$$outInforHash{$outfile}{'precursor_neutral_mass'}=$3;
			$$outInforHash{$outfile}{'charge'}=$4;

			# scan
			$outfile =~ m/\.(\d+)\.(\d+)\.(\d)$/;
			$$outInforHash{$outfile}{'scan'}=$1;
			
			$$outInforHash{$outfile}{'hitShown'}=0;
			next;
		}
		elsif (/^\<spectrum_query spectrum=\"(.*?)\" start_scan=\"(\d+)\" end_scan=\"\d+\" precursor_neutral_mass=\"([\d\.]+)\" assumed_charge=\"(\d+)\" index=\"\d+\"\>/)
		{
			$outfile=$1;
			$$outInforHash{$outfile}{'precursor_neutral_mass'}=$3;
			$$outInforHash{$outfile}{'charge'}=$4;

			# scan
			$outfile =~ m/\.(\d+)\.(\d+)\.(\d)$/;
			$$outInforHash{$outfile}{'scan'}=$1;
			
			$$outInforHash{$outfile}{'hitShown'}=0;
			next;
		}#else {die "Error unmatched pattern:\n$_";}
=cut		
		if (/^\<search_hit /)
		{
			$$outInforHash{$outfile}{'hitShown'}++;
			my @tmp = / (\w+)\=\"(.*?)\"/g;
			#print "$outfile: $tmp[0],$tmp[1],$tmp[2],$tmp[3]\n";
			for (my $i=0; $i<scalar(@tmp); $i=$i+2)
			{
				$$outInforHash{$outfile}{$tmp[$i]}[$$outInforHash{$outfile}{'hitShown'}]=$tmp[$i+1];
				$$paraHash{'hitRequiredParams'}{$tmp[$i]}='';
			}
			
			$altProCount=0;
			next;
		}
		
		if (/^\<alternative_protein protein=\"(.*?)\"\/\>/)
		{
			$altProCount++;
			$$outInforHash{$outfile}{'alternative_protein'}[$altProCount][$$outInforHash{$outfile}{'hitShown'}]=$1;
		
			next;
		}
		
		#modifications
		if (/^\<mod_aminoacid_mass position=\"(\d+)\" mass=\"([\d\.]+)\"\/\>/)
		{
			my $i;
			for ($i=1; $i<=$$paraHash{'dynamicMod'}{'count'}; $i++)
			{
				if ($2 == $$paraHash{'dynamicMod'}{'mass'}[$i])
				{
					unless (defined($$outInforHash{$outfile}{'dynamicMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}])) {$$outInforHash{$outfile}{'dynamicMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]=0;}
					$$outInforHash{$outfile}{'dynamicMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]++;
					$$outInforHash{$outfile}{'dynamicMod'}{'position'}[$$outInforHash{$outfile}{'dynamicMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]][$$outInforHash{$outfile}{'hitShown'}]=$1;
					$$outInforHash{$outfile}{'dynamicMod'}{'mass'}[$$outInforHash{$outfile}{'dynamicMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]][$$outInforHash{$outfile}{'hitShown'}]=$2;
					$$outInforHash{$outfile}{'dynamicMod'}{'modIndex'}[$$outInforHash{$outfile}{'dynamicMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]][$$outInforHash{$outfile}{'hitShown'}]=$i;
					last;
				}
			}
			if ($i>$$paraHash{'dynamicMod'}{'count'})
			{
				for ($i=1; $i<=$$paraHash{'staticMod'}{'count'}; $i++)
				{
					if ($2 == $$paraHash{'staticMod'}{'mass'}[$i])
					{
						unless (defined($$outInforHash{$outfile}{'staticMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}])) {$$outInforHash{$outfile}{'staticMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]=0;}
						$$outInforHash{$outfile}{'staticMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]++;
						$$outInforHash{$outfile}{'staticMod'}{'position'}[$$outInforHash{$outfile}{'staticMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]][$$outInforHash{$outfile}{'hitShown'}]=$1;
						$$outInforHash{$outfile}{'staticMod'}{'mass'}[$$outInforHash{$outfile}{'staticMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]][$$outInforHash{$outfile}{'hitShown'}]=$2;
						$$outInforHash{$outfile}{'staticMod'}{'modIndex'}[$$outInforHash{$outfile}{'staticMod'}{'found'}[$$outInforHash{$outfile}{'hitShown'}]][$$outInforHash{$outfile}{'hitShown'}]=$i;
						last;
					}
				}
			}
		
			next;
		}
		
		#optional parameters
		#if (/\<search_score name=\"(\w+)\" value=\"([\d\.]+)\"\/\>/) 
		#if (/\<search_score name=\"(\w+)\" value=\"([e\-\d\.]+)\"\/\>/) 
		if (/\<search_score name=\"(\w+)\" value=\"(.*?)\"\/\>/)  # changed on 4/12/15
		{
			$$outInforHash{$outfile}{$1}[$$outInforHash{$outfile}{'hitShown'}]=$2;  
			$$paraHash{'hitOptionalParams'}{$1}='';
			next;
		}

		# store $altProCount
		if ( /\<\/search_hit\>/ ) { $$outInforHash{$outfile}{'alternative_protein'}[0][$$outInforHash{$outfile}{'hitShown'}]=$altProCount; }
		
		last if /^\<\/msms_run_summary\>\Z/;
	}
	#print $$outInforHash{'adtmt1mg3.727.727.2'}{'deltacn'}[2];

	
	#reconstruct full peptide with dynamic modifications
	foreach $outfile (keys %{$outInforHash})
	{
		for (my $i=1; $i<=$$outInforHash{$outfile}{'hitShown'}; $i++)
		{
			if ( defined($$outInforHash{$outfile}{'dynamicMod'}{'found'}[$i]) and $$outInforHash{$outfile}{'dynamicMod'}{'found'}[$i]>=1)
			{
				#$i=1; my $j=1;
				#while ($i<=scalar() and )
				unless (defined($$outInforHash{$outfile}{'peptide'}[$i])) {die "Not defined peptide: $outfile, $i\n";}
				my @tmp=split(//,$$outInforHash{$outfile}{'peptide'}[$i]);
				#print "$outfile;@tmp;\n";
				for (my $j=1; $j<=$$outInforHash{$outfile}{'dynamicMod'}{'found'}[$i]; $j++)
				{
					#print "$$outInforHash{$outfile}{'dynamicMod'}{'position'}[$j][$i]\n";
					$tmp[$$outInforHash{$outfile}{'dynamicMod'}{'position'}[$j][$i]-1].=$$paraHash{'dynamicMod'}{'symble'}[$$outInforHash{$outfile}{'dynamicMod'}{'modIndex'}[$j][$i]];
				}

				
				$$outInforHash{$outfile}{'fullPeptide'}[$i]=join('',@tmp);#print "$outfile;@tmp;\n";

			}
			else {$$outInforHash{$outfile}{'fullPeptide'}[$i]=$$outInforHash{$outfile}{'peptide'}[$i];}
			#print "$outfile;$$outInforHash{$outfile}{'hitShown'};$$outInforHash{$outfile}{'fullPeptide'}[$i]\n";
		}
	}
	
	
	close IN;
	
}

