#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: idsum2::GenerateHtml_Tim07102014

package idsum2::GenerateHtml_Tim07102014;
#package idsum2::GenerateHtml;

use idsum2::CommonUtils;
#use CommonUtils;

sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);

  return $self;
}

### Tim Shaw Created Function
# 
#################################
sub gen_sumIDHtml {
	shift @_;
	
	my ($Sum_Hash, $Protein_Hash, $Peptide_Hash, $fprhash, $folders, $save_dir, $database, $grouphash) = @_;

	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";

	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	
	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	open (HTML, ">$save_dir/ID.html");
	#open (HTML, ">$save_dir/html/ID.html");
	#open (HTML, ">$save_dir/ID.html");
   	#print HTML "<HTML>\n";
 	#print HTML "<HEAD>\n";
	#print HTML "<TITLE>Multiple ID summary without Groups</TITLE>\n";
	open (TXT, ">$save_dir/text.txt");

	open (HTML2, ">$save_dir/ID_Redirect.html");
	
	#print_sumLegend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}	

	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	$redirect_path = "http://$server$remove/html/ID.html";
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	my ($groupnum, $subgroup_num);
	
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;		
	
	my $total_protein = scalar(keys %$Protein_Hash);
  	my $total_upeptide = 0;
  	for (sort {$$Peptide_Hash{$b}{'unique'} <=> $$Peptide_Hash{$a}{'unique'}} keys %$Peptide_Hash){
    		last if ($$Peptide_Hash{$_}{'unique'} == 0);
    		$total_upeptide++ if ($$Peptide_Hash{$_}{'unique'} == 1);
  	}
  	my $total_peptide = scalar(keys %$Peptide_Hash);
	my $total_SC = 0;
  	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
  	my $total_nomod = scalar(keys %peps_count);
	my $samples = scalar(@$folders);
	my $width = 80+($samples*10);
	
  	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=90%>\n";
  	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";
  	#print HTML "<TR>\n";
	
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC</font></TD>\n"; #changed to total peptide DMD 5/19/05
	
	#print HTML " </TR>\n";

        my $unique_fpr = 0;
        $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
        
    $uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'protein_fpr'};
    $pep_fdr = $$fprhash{'peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};

    #print HTML "</TD>\n"; 
  	#print HTML "<TR>\n";
  	
  	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";


	### these are the variables for the new html script
	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";

	my $CGILINK = "";

	$CGILINK .= "\$('td:eq(3)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[3] + '</a>');\n";
	
	
	#print HTML "<CENTER><A HREF=IDwG.html>Show Group Members</A></Center>";
	
	#print HTML "</TR></TD>\n";
  	#print HTML "<TR><TD Align=center>";
  	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  	#print HTML "</TD><TR>\n";
  	#print HTML "</TABLE>\n";

	$HTMLTABLEHEADER .= "<thead >\n";
  	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = $width%>\n";
  	
  	$HTMLTABLEHEADER .= "<tr>\n";
  	#print HTML "	<TR bgcolor=#EEEEEE>\n";
  	
  	my $index_i = 0;
  	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
  	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
	
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
  	#print HTML "		<TH width = 5%%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	print TXT "group;reference;totalpep;totalcoverage;";
	
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total SC</th>\n";
	#print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total SC</Font></B></TH>\n";
	for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_SC;";
		print TXT "$$folders[$i]_Coverage;";
		$index_i = $index_i + 1;
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">SC $$folders[$i]</th>\n";
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
  		#print HTML "		<TH width = 5%><B><Font Color = blue size=1>SC $$folders[$i]</font></Font></B></TH>\n";
	}

	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
    #print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total Peptide</Font></B></TH>\n";
    for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_TP;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">TP $$folders[$i]</th>\n";
		#print HTML "            <TH width = 5%><B><Font Color = blue size=1>TP $$folders[$i]</Font></B></TH>\n";
		print TXT "$$folders[$i]_Coverage;";
	}

	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
  	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
 	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
 	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "</tr>\n";
  	#print HTML "	</TR>\n";

	$HTMLTABLEHEADER .= "</thead>\n";
	print TXT "Description;KD;Length\n";

	### start writing the 
 	my %PrevProHash;
	my %array_abundance;
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{
		
		my $tab_group = "";
    	my $tab_reference = "";
    	my $tab_totalSC = "";
    	my $tab_SCtmt1 = "";
    	my $tab_SCtmt2 = "";
    	my $tab_totalpep = "";
    	my $tab_TPtmt1 = "";
    	my $tab_TPtmt2 = "";
    	my $tab_description = "";
    	my $tab_mass = "";
    	my $tab_abundance = "";
    
    	$HTMLTABLE .= "<tr class=\"gradeA\">\n";
		if($grouphash->{$protein}->{'order'}==1)
		{
			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};

			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			$mass = sprintf("%3.0f", $mass/1000);
		
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;

			#print HTML "<TR BGColor = #ccffcc>\n";
			if ($groupnum < 10){
				$tab_group = "000" . $groupnum;
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
				print TXT "000$groupnum;";
			}
			elsif ($groupnum < 100){
				$tab_group = "00" . $groupnum;
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
				print TXT "00$groupnum;";
			}
			elsif($groupnum < 1000) {
				$tab_group = "0" . $groupnum;
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
				print TXT "0$groupnum;";
			}
			else {
				$tab_group = $groupnum;
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
				print TXT "$groupnum;";
			}
			
			$HTMLTABLE .= "<td align=\"center\">" . $groupnum . "</td>\n";
			
			$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
			if (defined($2)){
				$tab_reference = $2;
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
				print TXT "$2;";
			} 
			else {
				$tab_reference = $protein;
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
				print TXT "$protein;";
			}

			$HTMLTABLE .= "<td align=\"center\">" . $tab_reference . "</td>\n";
			
			print TXT "$$Protein_Hash{$protein}{'total_nomod'};";
			my $Seq = $$Protein_Hash{$protein}{'sequence'};
			for my $peptide (sort keys %{$$Protein_Hash{$protein}{'peptides'}}){
				my $seq = $peptide;
				$seq =~ s/[\*\#\@]//g;
				$Seq =~ s/($seq)/\L$1/ig;
			}
			my $num = $Seq =~ s/([a-z])/$1/g;
			$$Protein_Hash{$protein}{'coverage'} = $num/length($Seq)*100;
			printf TXT "%.2f;", $num/length($Seq);

			$tab_totalSC = $$Protein_Hash{$protein}{'occurrence'};
			$HTMLTABLE .= "<td align=\"center\">" . $tab_totalSC . "</td>\n";
			#printf HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Protein_Hash{$protein}{'occurrence'}</B></center></Font></TD>\n";

			print TXT "$$Protein_Hash{$protein}{'occurrence'};";

			for (my $i=0; $i<scalar(@$folders); $i++){
				my $msdak = $$folders[$i];
				if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
					$HTMLTABLE .= "<td align=\"center\">0</td>\n";
					#$tab_SCtmt[$i] = 0;
					#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
					print TXT "0;";
				} else {
					$HTMLTABLE .= "<td align=\"center\">" . $$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'} . "</td>\n";
					#$tab_SCtmt[$i] = $$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'};
					#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</B></center></Font></TD>\n";
					print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'};";
				}
			}


			$tab_totalpep = $$Protein_Hash{$protein}{'total_nomod'};
			$HTMLTABLE .= "<td align=\"center\">" . $tab_totalpep . "</td>\n";
			#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
			print TXT "$$Protein_Hash{$protein}{'total_nomod'};";

		   for (my $i=0; $i<scalar(@$folders); $i++){
				my $msdak = $$folders[$i];
				if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
					#$tab_TPtmt[$i] = 0;
					$HTMLTABLE .= "<td align=\"center\">" . 0 . "</td>\n";		
					#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
					print TXT "0;";
					print TXT "0;";
				}
				else
				{
					my $Seq = $$Protein_Hash{$protein}{'sequence'};
					for my $peptide (sort keys %{$$Sum_Hash{$msdak}{'protein'}{$protein}{'peptides'}}){
						my $seq = $peptide;
						$seq =~ s/[\*\#\@]//g;
						$Seq =~ s/($seq)/\L$1/ig;
					}
					#$tab_TPtmt[$i] = $$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};
					$HTMLTABLE .= "<td align=\"center\">" . $$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'} . "</td>\n";
					#printf HTML "$totalcgi\n", $protein, "$save_dir/fractions/sum_$$folders[$i]", $database,$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};
					print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};";
					my $num = $Seq =~ s/([a-z])/$1/g;
					$$Sum_Hash{$msdak}{'protein'}{$protein}{'coverage'} = $num/length($Seq)*100;
					printf TXT "%.2f;", $num/length($Seq);
				}
			}

			if (!defined($description)){
				$tab_description = "&nbsp";
				#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
				print TXT ";";
			} else {
				$description =~ s/>/\&gt\;/g;
				$description =~ s/</\&lt\;/g;
				$tab_description = $description;
				#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
				print TXT "$description;";
			}
			
			$HTMLTABLE .= "<td align=\"center\">" . $tab_description . "</td>\n";	
			$tab_mass = $mass;
			$HTMLTABLE .= "<td align=\"center\">" . $mass . "</td>\n";	
			#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
			print TXT "$mass;";
			printf TXT "%d\n", length($Seq);
			$tab_abundance = $abundance;
			$HTMLTABLE .= "<td align=\"center\">" . $tab_abundance_rounded . "</td>\n";	
			#$HTMLTABLE .= "<td align=\"center\">" . $tab_abundance . "</td>\n";	
			#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
			$HTMLTABLE .= "</tr>\n";
		}
	}

	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "Multiple ID summary without Groups", "IDwG.html", "Show Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
  	close TXT;
   	print HTML $final_html;
   	close HTML;
  	
	return ($unique_fpr,$subgroup_num);
}

sub gen_IDHtml {
	shift @_;
	my ($Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $grouphash) = @_;

#foreach my $pro (keys %$grouphash) { print $pro,',', $grouphash->{$pro}->{group},';'; }

	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";
	
	
	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";
	
	open (HTML, ">$save_dir/ID.html");
  	#open (HTML, ">$save_dir/ID.html");
   	#print HTML "<HTML>\n";
   	#print HTML "<HEAD>\n";
  	#print HTML "<TITLE>ID summary w/o Groups</TITLE>\n";

	open (HTML2, ">$save_dir/ID_Redirect.html");
	
	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}	
	
	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	$redirect_path = "http://$server$remove/html/ID.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);	
	#print_Legend(*HTML);
	#print HTML "</HEAD>\n\n";

	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;	

	my $total_protein = scalar(keys %$Protein_Hash);
	my $total_upeptide = 0;
	foreach my $peptide (keys %$Peptide_Hash)
	{
		next if(!defined($$Peptide_Hash{$peptide}{'unique'}));
		last if ($$Peptide_Hash{$peptide}{'unique'} == 0);
		$total_upeptide++ if ($$Peptide_Hash{$peptide}{'unique'} == 1);
	}
	my $total_peptide = scalar(keys %$Peptide_Hash);

	my $total_SC = 0;
	
	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
	my $total_nomod = scalar(keys %peps_count);
	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=95%>\n";
	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";

	#print HTML "<TR>\n";
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC </font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

	#print HTML "<TR>\n";
	my ($overallrate, $rate1, $rate2, $rate3) = (0,0,0,0);
	my $unique_fpr = 0;
	$unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
	
	$uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'final_protein_fpr'};
    $pep_fdr = $$fprhash{'final_peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};

	#print HTML "</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";
	#print HTML "<TR>\n";
	
	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";

	### these are the variables for the new html script
	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGILINK = "";

 	#print HTML "<CENTER><A HREF=IDwG.html>Show Group Members</A></Center>";	
	
	#print HTML "</TR></TD>\n";
	#print HTML "<TR><TD Align=center>";
	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
	#print HTML "</TD><TR>\n";
	#print HTML "</TABLE>\n";

	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = 95%>\n";
	#print HTML "	<TR bgcolor=#EEEEEE>\n";
	
	$HTMLTABLEHEADER .= "<thead >\n";
	$HTMLTABLEHEADER .= "<tr>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
	#print HTML "		<TH width = 5%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">AX</th>\n";
		#print HTML "		<TH width = 2%><B><Font Color = blue>AX</Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Spectral Count</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Spectral Count</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Total Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Unique Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Unique Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Shared Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Shared Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
	#print HTML "		<TH width = 8%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">obsbl</th>\n";
		#print HTML "		<TH width = 2%><B><Font Color = blue>N<font size=1>obsbl</font></Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	$HTMLTABLEHEADER .= "</tr>\n";
	#print HTML "	</TR>\n";
	
	$HTMLTABLEHEADER .= "</thead>\n";

	$CGILINK .= "\$('td:eq(3)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[3] + '</a>');\n";
	$CGILINK .= "\$('td:eq(4)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[4] + '</a>')\n";
	$CGILINK .= "\$('td:eq(5)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[5] + '</a>')\n";	
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{
		if($grouphash->{$protein}->{'order'}==1)
		{
			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			my $SC = $$Protein_Hash{$protein}{'occurrence'};
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
			if (defined($$fprhash{'params'}))
			{
				if ($$fprhash{'params'}{'abundance_index'} == 1)
				{
					my $top = 0;
					if ($$fprhash{'params'}{'PAI_top'} eq "TP")
					{
						$top = $$Protein_Hash{$protein}{'total_nomod'};
					} else {
						$top = $SC;
					}
					my $PAI = $top/$$Protein_Hash{$protein}{'emPAI_ftpeps'};
				}
			}
		
			$mass = sprintf("%3.0f", $mass/1000);
		
			$HTMLTABLE .= "<tr class=\"gradeA\">\n";
			#print HTML "<TR BGColor = #ccffcc>\n";

		#added if loop DMD 5/19/05
			if ($groupnum < 10)
			{
				$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
			} elsif ($groupnum < 100)
			{
				$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
			} elsif($groupnum < 1000) 
			{
				$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
			} else 
			{
				$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
			}
			$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
			if (defined($2))
			{
				$HTMLTABLE .= "<td align=\"center\">" . $2 . "</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
			} else 
			{
				$HTMLTABLE .= "<td align=\"center\">" . $protein . "</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
			}
		
			$HTMLTABLE .= "<td align=\"center\">$SC</td>\n";
			#printf HTML "<TD><Font Color = blue><Center>%d</Center></Font></TD>\n", $SC;
			

			$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'total_nomod'} . "</td>\n";

			#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
			if (!$$Protein_Hash{$protein}{'unique'})
			{
				$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'unique_nomod'} . "</td>\n";	
				#print HTML "<TD><Font Size = 4 Color = Red><Center>$$Protein_Hash{$protein}{'unique_nomod'}</Center></Font></TD>\n";
			} else 
			{
				$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'unique_nomod'} . "</td>\n";
				#printf HTML "$uniquecgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'unique_nomod'};
			}
			if (!$$Protein_Hash{$protein}{'shared'})
			{
				$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'shared_nomod'} . "</td>\n";
				#print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'shared_nomod'}</Center></Font></TD>\n";
			} else 
			{
				$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'shared_nomod'} . "</td>\n";
				#printf HTML "$sharedcgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'shared_nomod'};
			}
			if (!defined($description))
			{
				$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
				#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
			} else 
			{
				$description =~ s/>/\&gt\;/g;
				$description =~ s/</\&lt\;/g;
				$HTMLTABLE .= "<td align=\"center\">$description</td>\n";
				#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
			}
			$HTMLTABLE .= "<td align=\"center\">" . $mass . "</td>\n";
			#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
			$HTMLTABLE .= "<td align=\"center\">" . $tab_abundance_rounded . "</td>\n";
			#$HTMLTABLE .= "<td align=\"center\">" . $abundance . "</td>\n";
			#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
			$HTMLTABLE .= "</tr>\n";
		}
		
	
	}
	
		
	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "ID summary w/o Groups", "IDwG.html", "Show Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
		
	return ($unique_fpr,$subgroup_num);	
}

sub gen_IDHtml_old {
	shift @_;
	my ($Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $grouphash) = @_;

	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	my $http = "HREF=http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";
	
  	open (HTML, ">$save_dir/ID.html");
   	print HTML "<HTML>\n";
   	print HTML "<HEAD>\n";
  	print HTML "<TITLE>ID summary w/o Groups</TITLE>\n";

	
	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}	
	
	print_Legend(*HTML);
	print HTML "</HEAD>\n\n";

	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;	

	my $total_protein = scalar(keys %$Protein_Hash);
	my $total_upeptide = 0;
	foreach my $peptide (keys %$Peptide_Hash)
	{
		next if(!defined($$Peptide_Hash{$peptide}{'unique'}));
		last if ($$Peptide_Hash{$peptide}{'unique'} == 0);
		$total_upeptide++ if ($$Peptide_Hash{$peptide}{'unique'} == 1);
	}
	my $total_peptide = scalar(keys %$Peptide_Hash);

	my $total_SC = 0;
	
	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
	my $total_nomod = scalar(keys %peps_count);
	print HTML "<BODY>\n";
	print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=95%>\n";
	print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";

	print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	print HTML "SC = $total_SC </font></TD>\n"; #changed to total peptide DMD 5/19/05
	print HTML " </TR>\n";

	print HTML "<TR>\n";
	my ($overallrate, $rate1, $rate2, $rate3) = (0,0,0,0);
	my $unique_fpr = 0;
	$unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
    printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};

	print HTML "</font></TD>\n"; #changed to total peptide DMD 5/19/05
	print HTML " </TR>\n";
	print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	print HTML " </TR>\n";
	print HTML "<TR><TD align=center>\n";

 	print HTML "<CENTER><A HREF=IDwG.html>Show Group Members</A></Center>";	
	
	print HTML "</TR></TD>\n";
	print HTML "<TR><TD Align=center>";
	print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
	print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
	print HTML "</TD><TR>\n";
	print HTML "</TABLE>\n";

	print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = 95%>\n";
	print HTML "	<TR bgcolor=#EEEEEE>\n";
	print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
	print HTML "		<TH width = 5%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		print HTML "		<TH width = 2%><B><Font Color = blue>AX</Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Spectral Count</Font></B></TH>\n";
	print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Total Peptide</Font></B></TH>\n";
	print HTML "		<TH width = 7%><B><Font Color = blue size=1>Unique Peptide</Font></B></TH>\n";
	print HTML "		<TH width = 7%><B><Font Color = blue size=1>Shared Peptide</Font></B></TH>\n";
	print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
	print HTML "		<TH width = 8%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
	print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		print HTML "		<TH width = 2%><B><Font Color = blue>N<font size=1>obsbl</font></Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	print HTML "	</TR>\n";
	
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{
		if($grouphash->{$protein}->{'order'}==1)
		{
			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			my $SC = $$Protein_Hash{$protein}{'occurrence'};
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
			if (defined($$fprhash{'params'}))
			{
				if ($$fprhash{'params'}{'abundance_index'} == 1)
				{
					my $top = 0;
					if ($$fprhash{'params'}{'PAI_top'} eq "TP")
					{
						$top = $$Protein_Hash{$protein}{'total_nomod'};
					} else {
						$top = $SC;
					}
					my $PAI = $top/$$Protein_Hash{$protein}{'emPAI_ftpeps'};
				}
			}
		
			$mass = sprintf("%3.0f", $mass/1000);
		
			print HTML "<TR BGColor = #ccffcc>\n";

		#added if loop DMD 5/19/05
			if ($groupnum < 10)
			{
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
			} elsif ($groupnum < 100)
			{
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
			} elsif($groupnum < 1000) 
			{
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
			} else 
			{
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
			}
			$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
			if (defined($2))
			{
				print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
			} else 
			{
				print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
			}
		#	if ($abundance != 0){
		#		printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		#	}
			printf HTML "<TD><Font Color = blue><Center>%d</Center></Font></TD>\n", $SC;


			printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
			if (!$$Protein_Hash{$protein}{'unique'})
			{
				print HTML "<TD><Font Size = 4 Color = Red><Center>$$Protein_Hash{$protein}{'unique_nomod'}</Center></Font></TD>\n";
			} else 
			{
				printf HTML "$uniquecgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'unique_nomod'};
			}
			if (!$$Protein_Hash{$protein}{'shared'})
			{
				print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'shared_nomod'}</Center></Font></TD>\n";
			} else 
			{
				printf HTML "$sharedcgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'shared_nomod'};
			}
			if (!defined($description))
			{
				print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
			} else 
			{
				$description =~ s/>/\&gt\;/g;
				$description =~ s/</\&lt\;/g;
				print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
			}
			print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
			printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $tab_abundance_rounded;
			#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		}
	}
	return ($unique_fpr,$subgroup_num);	
}

sub gen_IDwGHtml{
	shift @_;
	my ($Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $grouphash) = @_;

	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";


	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	my $http = "http://$server/cgi-bin";
	#my $http = "HREF=http://$server/cgi-bin";

	
	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	open (HTML, ">$save_dir/IDwG.html");
	#open (HTML, ">$save_dir/html/IDwG.html");
   	 #open (HTML, ">$save_dir/IDwG.html");
   	 #print HTML "<HTML>\n";
   	 #print HTML "<HEAD>\n";
   	 #print HTML "<TITLE>ID Summary w/ Groups</TITLE>\n";
	 
	open (HTML2, ">$save_dir/IDwG_Redirect.html");
	
	#print_Legend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}		
	
	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	
	$redirect_path = "http://$server$remove/IDwG.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;	

	my $total_protein = scalar(keys %$Protein_Hash);
	my $total_upeptide = 0;
	foreach my $peptide (keys %$Peptide_Hash)
	{
		next if(!defined($$Peptide_Hash{$peptide}{'unique'}));
		last if ($$Peptide_Hash{$peptide}{'unique'} == 0);
		$total_upeptide++ if ($$Peptide_Hash{$peptide}{'unique'} == 1);
	}
	my $total_peptide = scalar(keys %$Peptide_Hash);

	my $total_SC = 0;
	
	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
	my $total_nomod = scalar(keys %peps_count);
	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=95%>\n";
	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";

	#print HTML "<TR>\n";
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC </font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

	#print HTML "<TR>\n";
	my ($overallrate, $rate1, $rate2, $rate3) = (0,0,0,0);

        my $unique_fpr = 0;
        $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
     $uniqpro_fdr = $unique_fpr;
     $pro_fdr = $$fprhash{'final_protein_fpr'};
     $pep_fdr = $$fprhash{'final_peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};


	#print HTML "</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";
	#print HTML "<TR>\n";
	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";

	### these are the variables for the new html script
	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGILINK = "";

	$CGILINK .= "\$('td:eq(3)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[3] + '</a>');\n";
	$CGILINK .= "\$('td:eq(4)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[4] + '</a>')\n";
	$CGILINK .= "\$('td:eq(5)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[5] + '</a>')\n";
	#print HTML "<CENTER><A HREF=ID.html>Hide Group Members</A></Center>";

	$HTMLTABLEHEADER .= "<thead >\n";	
	$HTMLTABLEHEADER .= "<tr>\n";
	#print HTML "</TR></TD>\n";
	#print HTML "<TR><TD Align=center>";
	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
	#print HTML "</TD><TR>\n";
	#print HTML "</TABLE>\n";

	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = 95%>\n";
	#print HTML "	<TR bgcolor=#EEEEEE>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
	#print HTML "		<TH width = 5%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#print HTML "		<TH width = 2%><B><Font Color = blue>AX</Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Spectral Count</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Spectral Count</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Total Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Unique Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Unique Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Shared Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Shared Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
	#print HTML "		<TH width = 8%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
	
	$HTMLTABLEHEADER .= "</tr>\n";
	$HTMLTABLEHEADER .= "</thead>\n";
	
	if (defined($$fprhash{'params'}))
	{
		#print HTML "		<TH width = 2%><B><Font Color = blue>N<font size=1>obsbl</font></Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	#print HTML "	</TR>\n";
	
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{		
		my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
		my $description = $$Protein_Hash{$protein}{'annotation'};
		my $mass = $$Protein_Hash{$protein}{'MW'};
		my $SC = $$Protein_Hash{$protein}{'occurrence'};
		my $abundance = $$Protein_Hash{$protein}{'abundance'};
		my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
	
		$HTMLTABLE .= "<tr class=\"gradeA\">\n";
		
		if (defined($$fprhash{'params'}))
		{
			if ($$fprhash{'params'}{'abundance_index'} == 1)
			{
				my $top = 0;
				if ($$fprhash{'params'}{'PAI_top'} eq "TP")
				{
					$top = $$Protein_Hash{$protein}{'total_nomod'};
				} else {
					$top = $SC;
				}
				my $PAI = $top/$$Protein_Hash{$protein}{'emPAI_ftpeps'};
			}
		}
	
		$mass = sprintf("%3.0f", $mass/1000);
		if($grouphash->{$protein}->{'order'}==1)
		{			 
			#print HTML "<TR BGColor = #ccffcc>\n";
		}
		else
		{
			#print HTML "<TR>\n";
		}

		#added if loop DMD 5/19/05
		if ($groupnum < 10)
		{
			$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
		} elsif ($groupnum < 100)
		{
			$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
		} elsif($groupnum < 1000) 
		{
			$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
		}
		$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
		if (defined($2))
		{
			$HTMLTABLE .= "<td align=\"center\">$2</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$protein</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
		}
		$HTMLTABLE .= "<td align=\"center\">$SC</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%d</Center></Font></TD>\n", $SC;

	## added by yanji, changed path
		
		$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'total_nomod'} . "</td>\n";
		#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
		if (!$$Protein_Hash{$protein}{'unique'})
		{
			$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'unique_nomod'} . "</td>\n";
			#print HTML "<TD><Font Size = 4 Color = Red><Center>$$Protein_Hash{$protein}{'unique_nomod'}</Center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'unique_nomod'} . "</td>\n";
			#printf HTML "$uniquecgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'unique_nomod'};
		}
		if (!$$Protein_Hash{$protein}{'shared'})
		{
			$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'shared_nomod'} . "</td>\n";
			#print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'shared_nomod'}</Center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">" . $$Protein_Hash{$protein}{'shared_nomod'} . "</td>\n";
			#printf HTML "$sharedcgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'shared_nomod'};
		}
		if (!defined($description))
		{
			$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
			#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
		} else 
		{
			$description =~ s/>/\&gt\;/g;
			$description =~ s/</\&lt\;/g;
			$HTMLTABLE .= "<td align=\"center\">" . $description . "</td>\n";
			#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
		}
		$HTMLTABLE .= "<td align=\"center\">" . $mass . "</td>\n";
		#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
		$HTMLTABLE .= "<td align=\"center\">" . $tab_abundance_rounded . "</td>\n";
		#$HTMLTABLE .= "<td align=\"center\">" . $abundance . "</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		$HTMLTABLE .= "</tr>\n";
	}
	
	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "ID Summary w/ Groups", "ID.html", "Hide Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
	
}

sub gen_IDmodHtml{
	shift @_;
	my ($Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $grouphash) = @_;


	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";
	
	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";
	
	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	open (HTML, ">$save_dir/IDmod.html");
    #open (HTML, ">$save_dir/IDmod.html");
    #print HTML "<HTML>\n";
    #print HTML "<HEAD>\n";
    #print HTML "<TITLE>IDMOD summary w/o Groups</TITLE>\n";

	open (HTML2, ">$save_dir/IDmod_Redirect.html");
	
	#print_Legend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}
	
	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	
	$redirect_path = "http://$server/$current_user/$remove/IDmod.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;	

	my $total_protein = scalar(keys %$Protein_Hash);
	my $total_upeptide = 0;
	foreach my $peptide (keys %$Peptide_Hash)
	{
		next if(!defined($$Peptide_Hash{$peptide}{'unique'}));
		last if ($$Peptide_Hash{$peptide}{'unique'} == 0);
		$total_upeptide++ if ($$Peptide_Hash{$peptide}{'unique'} == 1);
	}
	my $total_peptide = scalar(keys %$Peptide_Hash);

	my $total_SC = 0;
	
	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
	my $total_nomod = scalar(keys %peps_count);
	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=95%>\n";
	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";

	#print HTML "<TR>\n";
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC </font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

	#print HTML "<TR>\n";

    my $unique_fpr = 0;
    $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
        
    $uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'final_protein_fpr'};
    $pep_fdr = $$fprhash{'final_peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};

	#print HTML "</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";
	#print HTML "<TR>\n";
	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";

	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGISPECIAL = "";

	#print HTML "<CENTER><A HREF=IDwGmod.html>Show Group Members</A></Center>";
	
	#print HTML "</TR></TD>\n";
	#print HTML "<TR><TD Align=center>";
	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
	#print HTML "</TD><TR>\n";
	#print HTML "</TABLE>\n";

	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = 95%>\n";
	#print HTML "	<TR bgcolor=#EEEEEE>\n";
	
	$HTMLTABLEHEADER .= "<thead >\n";
	$HTMLTABLEHEADER .= "<tr>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
	#print HTML "		<TH width = 5%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#print HTML "		<TH width = 2%><B><Font Color = blue>AX</Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Spectral Count</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Spectral Count</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Total Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Unique Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Unique Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Shared Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Shared Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
	#print HTML "		<TH width = 8%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#print HTML "		<TH width = 2%><B><Font Color = blue>N<font size=1>obsbl</font></Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	#print HTML "	</TR>\n";
		$HTMLTABLEHEADER .= "</tr>\n";
	#print HTML "	</TR>\n";
	
	$HTMLTABLEHEADER .= "</thead>\n";

	$CGILINK .= "\$('td:eq(3)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[3] + '</a>');\n";
	$CGILINK .= "\$('td:eq(4)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[4] + '</a>')\n";
	$CGILINK .= "\$('td:eq(5)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[5] + '</a>')\n";
	
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{
		if($grouphash->{$protein}->{'order'}==1)
		{
		my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
		my $description = $$Protein_Hash{$protein}{'annotation'};
		my $mass = $$Protein_Hash{$protein}{'MW'};
		my $SC = $$Protein_Hash{$protein}{'occurrence'};
		my $abundance = $$Protein_Hash{$protein}{'abundance'};
		my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
		if (defined($$fprhash{'params'}))
		{
			if ($$fprhash{'params'}{'abundance_index'} == 1)
			{
				my $top = 0;
				if ($$fprhash{'params'}{'PAI_top'} eq "TP")
				{
					$top = $$Protein_Hash{$protein}{'total_nomod'};
				} else {
					$top = $SC;
				}
				my $PAI = $top/$$Protein_Hash{$protein}{'emPAI_ftpeps'};
			}
		}
	
		$mass = sprintf("%3.0f", $mass/1000);
	
		if($grouphash->{$protein}->{'order'}==1)
		{			 
			#print HTML "<TR BGColor = #ccffcc>\n";
		}
		else
		{
			#print HTML "<TR>\n";
		}
		$HTMLTABLE .= "<tr class=\"gradeA\">\n";
		#added if loop DMD 5/19/05
		if ($groupnum < 10)
		{
			$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
		} elsif ($groupnum < 100)
		{
			$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
		} elsif($groupnum < 1000) 
		{
			$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
		}
		$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
		if (defined($2))
		{
			$HTMLTABLE .= "<td align=\"center\">$2</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$protein</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
		}
		$HTMLTABLE .= "<td align=\"center\">$SC</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%d</Center></Font></TD>\n", $SC;

		$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'total_nomod'}</td>\n";
		#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
		if (!$$Protein_Hash{$protein}{'unique'})
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'unique_nomod'}</td>\n";
			#print HTML "<TD><Font Size = 4 Color = Red><Center>$$Protein_Hash{$protein}{'unique_nomod'}</Center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'unique_nomod'}</td>\n";
			#printf HTML "$uniquecgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'unique_nomod'};
		}
		if (!$$Protein_Hash{$protein}{'shared'})
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'shared_nomod'}</td>\n";
			#print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'shared_nomod'}</Center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'shared_nomod'}</td>\n";
			#printf HTML "$sharedcgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'shared_nomod'};
		}
		if (!defined($description))
		{
			$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
			#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
		} else 
		{
			$description =~ s/>/\&gt\;/g;
			$description =~ s/</\&lt\;/g;
			$HTMLTABLE .= "<td align=\"center\">$description</td>\n";
			#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
		}
		$HTMLTABLE .= "<td align=\"center\">$mass</td>\n";
		#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
		$HTMLTABLE .= "<td align=\"center\">$tab_abundance_rounded</td>\n";
		#$HTMLTABLE .= "<td align=\"center\">$abundance</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		
		$HTMLTABLE .= "</tr>\n";
		}
	}
	
	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "IDMOD summary w/o Groups", "IDwGmod.html", "Show Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
	
}


	
		

sub gen_IDwGmodHtml{
	shift @_;
	my ($Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $grouphash) = @_;

	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";

	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";
	
	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	open (HTML, ">$save_dir/IDwGmod.html");
    #open (HTML, ">$save_dir/IDwGmod.html");
	
    #print HTML "<HTML>\n";
    #print HTML "<HEAD>\n";
    #print HTML "<TITLE>IDMOD Summary w/ Groups</TITLE>\n";
		
	open (HTML2, ">$save_dir/IDwGmod_Redirect.html");
	
	#print_Legend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}
	
	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	
	$redirect_path = "http://$server$remove/IDwGmod.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;	

	my $total_protein = scalar(keys %$Protein_Hash);
	my $total_upeptide = 0;
	foreach my $peptide (keys %$Peptide_Hash)
	{
		next if(!defined($$Peptide_Hash{$peptide}{'unique'}));
		last if ($$Peptide_Hash{$peptide}{'unique'} == 0);
		$total_upeptide++ if ($$Peptide_Hash{$peptide}{'unique'} == 1);
	}
	my $total_peptide = scalar(keys %$Peptide_Hash);

	my $total_SC = 0;
	
	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
	my $total_nomod = scalar(keys %peps_count);
	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=95%>\n";
	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";

	#print HTML "<TR>\n";
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC </font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

	#print HTML "<TR>\n";
        my $unique_fpr = 0;
        $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
    $uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'final_protein_fpr'};
    $pep_fdr = $$fprhash{'final_peptide_fpr'};
    
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};


	#print HTML "</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";
	#print HTML "<TR>\n";
	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";

	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGISPECIAL = "";

	#print HTML "<CENTER><A HREF=IDmod.html>Hide Group Members</A></Center>";
	
	#print HTML "</TR></TD>\n";
	#print HTML "<TR><TD Align=center>";
	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
	#print HTML "</TD><TR>\n";
	#print HTML "</TABLE>\n";

	$HTMLTABLEHEADER .= "<thead >\n";
	$HTMLTABLEHEADER .= "<tr>\n";
	
	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = 95%>\n";
	#print HTML "	<TR bgcolor=#EEEEEE>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
	#print HTML "		<TH width = 5%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#print HTML "		<TH width = 2%><B><Font Color = blue>AX</Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Spectral Count</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Spectral Count</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Total Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Unique Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Unique Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Shared Peptide</th>\n";
	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Shared Peptide</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
	#print HTML "		<TH width = 8%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
	if (defined($$fprhash{'params'}))
	{
		#print HTML "		<TH width = 2%><B><Font Color = blue>N<font size=1>obsbl</font></Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
	#print HTML "	</TR>\n";
	$HTMLTABLEHEADER .= "</tr>\n";
	#print HTML "	</TR>\n";
	
	$HTMLTABLEHEADER .= "</thead>\n";

	$CGILINK .= "\$('td:eq(3)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[3] + '</a>');\n";
	$CGILINK .= "\$('td:eq(4)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[4] + '</a>')\n";
	$CGILINK .= "\$('td:eq(5)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[5] + '</a>')\n";
	
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{		
		my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};
		my $description = $$Protein_Hash{$protein}{'annotation'};
		my $mass = $$Protein_Hash{$protein}{'MW'};
		my $SC = $$Protein_Hash{$protein}{'occurrence'};
		my $abundance = $$Protein_Hash{$protein}{'abundance'};
		my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
		if (defined($$fprhash{'params'}))
		{
			if ($$fprhash{'params'}{'abundance_index'} == 1)
			{
				my $top = 0;
				if ($$fprhash{'params'}{'PAI_top'} eq "TP")
				{
					$top = $$Protein_Hash{$protein}{'total_nomod'};
				} else {
					$top = $SC;
				}
				my $PAI = $top/$$Protein_Hash{$protein}{'emPAI_ftpeps'};
			}
		}
	
		$mass = sprintf("%3.0f", $mass/1000);
		if($grouphash->{$protein}->{'order'}==1)
		{			 
			#print HTML "<TR BGColor = #ccffcc>\n";
		}
		else
		{
			#print HTML "<TR>\n";
		}

		$HTMLTABLE .= "<tr class=\"gradeA\">\n";
		
		#added if loop DMD 5/19/05
		if ($groupnum < 10)
		{
			$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
		} elsif ($groupnum < 100)
		{
			$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
		} elsif($groupnum < 1000) 
		{
			$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
		}
		$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
		if (defined($2))
		{
			$HTMLTABLE .= "<td align=\"center\">$2</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$protein</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
		}

		$HTMLTABLE .= "<td align=\"center\">$SC</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%d</Center></Font></TD>\n", $SC;


		$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'total_nomod'}</td>\n";
		#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
		if (!$$Protein_Hash{$protein}{'unique'})
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'unique_nomod'}</td>\n";
			#print HTML "<TD><Font Size = 4 Color = Red><Center>$$Protein_Hash{$protein}{'unique_nomod'}</Center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'unique_nomod'}</td>\n";
			#printf HTML "$uniquecgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'unique_nomod'};
		}
		if (!$$Protein_Hash{$protein}{'shared'})
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'shared_nomod'}</td>\n";
			#print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'shared_nomod'}</Center></Font></TD>\n";
		} else 
		{
			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'shared_nomod'}</td>\n";
			#printf HTML "$sharedcgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'shared_nomod'};
		}
		if (!defined($description))
		{
			$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
			#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
		} else 
		{
			$description =~ s/>/\&gt\;/g;
			$description =~ s/</\&lt\;/g;
			$HTMLTABLE .= "<td align=\"center\">$description</td>\n";
			#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
		}
		$HTMLTABLE .= "<td align=\"center\">$mass</td>\n";
		#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
		$HTMLTABLE .= "<td align=\"center\">$tab_abundance_rounded</td>\n";
		#$HTMLTABLE .= "<td align=\"center\">$abundance</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		$HTMLTABLE .= "</tr>\n";
	}
	
	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "IDMOD Summary w/ Groups", "IDmod.html", "Hide Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
	
}


sub print_Legend{
	(*HTML) = @_;
	print HTML "<Script Language=JavaScript>\n";
	print HTML "<!--\n";
	print HTML "  function mainopen() {\n";
	print HTML "    descript = window.open(\"\",\"ProteinID\",\"width=300,height=450,top=200,left=700\");\n";
	print HTML "    with(descript.document){\n";
	print HTML "      write(\"<html><head>\");\n";
	print HTML "      write(\"<title>Protein Identification Legend</title>\")\n";
	print HTML "      write(\"</head><body bgcolor = #ccffcc>\");\n";
	print HTML "      write(\"<Font Size = 2>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"Acquired MS/MS spectra are searched against a database to match peptides. The assigned peptides are filtered by mass accuracy and matching scores. The summary table lists the total numbers of spectra counts, peptides, proteins, unique proteins and groups. The target-decoy strategy is used to evaluate false discovery rate (FDR) of unique proteins and peptides. Follow the links to find protein/peptide sequences.<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
        print HTML "      write(\"<Font Color = blue><B><U>Unique Protein FDR</U></B></Font><BR>\");\n";
	print HTML "      write(\"Unique Protein FDR = # decoy proteins / # unique target proteins.<BR>\");\n";
        #print HTML "      write(\"Unique Protein FDR = #decoy proteins / (#unique target proteins - #decoy proteins).<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Protein FDR</U></B></Font><BR>\");\n";
	print HTML "      write(\"Protein FDR = # decoy proteins / # target proteins.<BR>\");\n";
	#print HTML "      write(\"Protein FDR = # decoy proteins / (# target proteins - # decoy proteins).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Peptide FDR</U></B></Font><BR>\");\n";
	print HTML "      write(\"Peptide FDR = # decoy peptides / # target peptides.<BR>\");\n";
	#print HTML "      write(\"Peptide FDR = # decoy peptides / (# target peptides - # decoy peptides).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Peptides without C and M modifications</U></B></Font><BR>\");\n";
	print HTML "      write(\"Total peptide number without counting Cys and/or Met modifications.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Show/Hide Group Members</U></B></Font><BR>\");\n";
	print HTML "      write(\"Follow the links to show/hide protein group members.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Group</U></B></Font><BR>\");\n";
	print HTML "      write(\"Identified proteins that share one or more peptides are clustered into a single protein group (e.g. 0001). The group is represented by one protein with maximal spectra count (e.g. 0001.1). If other proteins in the same group are assigned with one or more uniquely identified peptides, these proteins are also counted as unique proteins (e.g. 0001.2).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Reference</U></B></Font><BR>\");\n";
	print HTML "      write(\"The accession number of identified proteins in the database.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Spectral Count</U></B></Font><BR>\");\n";
	print HTML "      write(\"The summed number of spectral counts (i.e. MS/MS scans) assigned to one identified protein.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Peptide</U></B></Font><BR>\");\n";
	print HTML "      write(\"The number of peptides assigned to one identified protein. As one peptide can be identified by multiple spectra counts, the spectral count is always equal to or larger than the total peptide number.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Unique Peptide</U></B></Font><BR>\");\n";
	print HTML "      write(\"The number of peptides uniquely assigned to one identified protein.<BR>\");\n";
	#print HTML "      write(\"<Font Color = red><B>ONLY</B></Font> with the given protein.  \");\n";

	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Shared Peptide</U></B></Font><BR>\");\n";
	print HTML "      write(\"The number of peptides shared with other proteins in a group.<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Description</U></B></Font><BR>\");\n";
	print HTML "      write(\"Protein annotation in in the database.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Mass (KD)</U></B></Font><BR>\");\n";
	print HTML "      write(\"Theoretical mass (KD) of one protein, based on the amino acid sequence in the database.<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Abundance</U></B></Font><BR>\");\n";
	print HTML "      write(\"After normalizing protein size, the SC is roughly correlated with protein abundance in the sample. Thus, abundance Index = SC x 50 (kDa) / protein size (kDa).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"</Font>\");\n";
	print HTML "      write(\"</body></html>\");\n";
	print HTML "    }\n";
	print HTML "    descript.document.close();\n";
	print HTML "  }\n";
	print HTML "-->\n";
	print HTML "</Script>\n";
}


sub print_sumLegend{
	(*HTML) = @_;
	print HTML "<Script Language=JavaScript>\n";
	print HTML "<!--\n";
	print HTML "  function mainopen() {\n";
	print HTML "    descript = window.open(\"\",\"ProteinID\",\"width=300,height=200,top=200,left=700\");\n";
	print HTML "    with(descript.document){\n";
	print HTML "      write(\"<html><head>\");\n";
	print HTML "      write(\"<title>Protein Identification Legend</title>\")\n";
	print HTML "      write(\"</head><body bgcolor = #ccffcc>\");\n";
	print HTML "      write(\"<Font Size = 2>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"Acquired MS/MS spectra are searched against a database to match peptides. The assigned peptides are filtered by mass accuracy and matching scores. The summary table lists the total numbers of spectra counts, peptides, proteins, unique proteins and groups. The target-decoy strategy is used to evaluate false discovery rate (FDR) of unique proteins and peptides. Follow the links to find protein/peptide sequences.<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
        print HTML "      write(\"<Font Color = blue><B><U>Unique Protein FDR</U></B></Font><BR>\");\n";
        print HTML "      write(\"Unique Protein FDR = #decoy proteins / # unique target proteins.<BR>\");\n";
	#print HTML "      write(\"Unique Protein FDR = #decoy proteins / (#unique target proteins - #decoy proteins).<BR>\");\n";
	
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Protein FDR</U></B></Font><BR>\");\n";
	print HTML "      write(\"Protein FDR = # decoy proteins / # target proteins.<BR>\");\n";
	#print HTML "      write(\"Protein FDR = # decoy proteins / (# target proteins - # decoy proteins).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Peptide FDR</U></B></Font><BR>\");\n";
	print HTML "      write(\"Peptide FDR = # decoy peptides / # target peptides.<BR>\");\n";
	#print HTML "      write(\"Peptide FDR = # decoy peptides / (# target peptides - # decoy peptides).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Peptides without C and M modifications</U></B></Font><BR>\");\n";
	print HTML "      write(\"Total peptide number without counting Cys and/or Met modifications.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Show/Hide Group Members</U></B></Font><BR>\");\n";
	print HTML "      write(\"Follow the links to show/hide protein group members.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Group</U></B></Font><BR>\");\n";
	print HTML "      write(\"Identified proteins that share one or more peptides are clustered into a single protein group (e.g. 0001). The group is represented by one protein with maximal spectra count (e.g. 0001.1). If other proteins in the same group are assigned with one or more uniquely identified peptides, these proteins are also counted as unique proteins (e.g. 0001.2).<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Reference</U></B></Font><BR>\");\n";
	print HTML "      write(\"The accession number of identified proteins in the database.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>SC</U></B></Font><BR>\");\n";
	print HTML "      write(\"The total number of spectral counts (i.e. MS/MS scans) assigned to one identified protein from all samples.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Total TP</U></B></Font><BR>\");\n";
	print HTML "      write(\"The total number of peptides assigned to one identified protein from all samples. As one peptide can be identified by multiple spectra counts, the total SC is always equal to or larger than the total TP.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Description</U></B></Font><BR>\");\n";
	print HTML "      write(\"Protein annotation in in the database.<BR>\");\n";
	
	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Mass (KD)</U></B></Font><BR>\");\n";
	print HTML "      write(\"Theoretical mass (KD) of one protein, based on the amino acid sequence in the database.<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"<Font Color = blue><B><U>Abundance</U></B></Font><BR>\");\n";
	print HTML "      write(\"After normalizing protein size, the SC is roughly correlated with protein abundance in the sample. Thus, abundance Index = SC x 50 (kDa) / protein size (kDa).<BR>\");\n";

	print HTML "      write(\"<BR>\");\n";
	print HTML "      write(\"</body></html>\");\n";
	print HTML "    }\n";
	print HTML "    descript.document.close();\n";
	print HTML "  }\n";
	print HTML "-->\n";
	print HTML "</Script>\n";
}



sub gen_sumIDHtml_old {
	shift @_;
	
	my ($Sum_Hash, $Protein_Hash, $Peptide_Hash, $fprhash, $folders, $save_dir, $database, $grouphash) = @_;

	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	my $http = "HREF=http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	
	open (HTML, ">$save_dir/ID.html");
   	print HTML "<HTML>\n";
 	print HTML "<HEAD>\n";
	print HTML "<TITLE>Multiple ID summary without Groups</TITLE>\n";
	open (TXT, ">$save_dir/text.txt");

	print_sumLegend(*HTML);
	print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}	

	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;		
	
	my $total_protein = scalar(keys %$Protein_Hash);
  	my $total_upeptide = 0;
  	for (sort {$$Peptide_Hash{$b}{'unique'} <=> $$Peptide_Hash{$a}{'unique'}} keys %$Peptide_Hash){
    		last if ($$Peptide_Hash{$_}{'unique'} == 0);
    		$total_upeptide++ if ($$Peptide_Hash{$_}{'unique'} == 1);
  	}
  	my $total_peptide = scalar(keys %$Peptide_Hash);
	my $total_SC = 0;
  	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
  	my $total_nomod = scalar(keys %peps_count);
	my $samples = scalar(@$folders);
	my $width = 80+($samples*10);
  	print HTML "<BODY>\n";
	print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=90%>\n";
  	print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";
  	print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	print HTML "SC = $total_SC</font></TD>\n"; #changed to total peptide DMD 5/19/05
	print HTML " </TR>\n";

        my $unique_fpr = 0;
        $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
    printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};

#    printf HTML "<TD Align=center><Font Size=2>Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};
    print HTML "</TD>\n"; 
  	print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	print HTML " </TR>\n";
	print HTML "<TR><TD align=center>\n";
	
	print HTML "<CENTER><A HREF=IDwG.html>Show Group Members</A></Center>";
	

	print HTML "</TR></TD>\n";
  	print HTML "<TR><TD Align=center>";
  	print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  	print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  	print HTML "</TD><TR>\n";
  	print HTML "</TABLE>\n";

  	print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = $width%>\n";
  	print HTML "	<TR bgcolor=#EEEEEE>\n";
  	print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";

  	print HTML "		<TH width = 5%%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
### changed on 8/9/2012 
# print HTML "          <TH width = 5%%><B><Font Color = blue>Abundance</Font></B></TH>\n";
########################
#  print HTML "		<TH width = 2%><B><Font Color = blue>Total<br><font size=1>Peptides</font></Font></B></TH>\n";
	print TXT "group;reference;totalpep;totalcoverage;";
	print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total SC</Font></B></TH>\n";
	for (my $i=0; $i<scalar(@$folders); $i++){
#		print TXT "$$folders[$i]_TP;";
		print TXT "$$folders[$i]_SC;";
		print TXT "$$folders[$i]_Coverage;";
#  		print HTML "		<TH width = 2%><B><Font Color = blue>TP<br><font size=1>$$folders[$i]</font></Font></B></TH>\n";
  		print HTML "		<TH width = 5%><B><Font Color = blue size=1>SC $$folders[$i]</font></Font></B></TH>\n";
	}

    print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total Peptide</Font></B></TH>\n";
    for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_TP;";
		print HTML "            <TH width = 5%><B><Font Color = blue size=1>TP $$folders[$i]</Font></B></TH>\n";
		print TXT "$$folders[$i]_Coverage;";
	}

  	print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
 	print HTML "		<TH width = 7%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
###### added by yanji
	print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
######
  	print HTML "	</TR>\n";

	print TXT "Description;KD;Length\n";

	
 	my %PrevProHash;
	my %array_abundance;
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{
		if($grouphash->{$protein}->{'order'}==1)
		{
			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};

			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			$mass = sprintf("%3.0f", $mass/1000);
		
	########### added by yanji
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
	######	


			print HTML "<TR BGColor = #ccffcc>\n";
		#added if loop DMD 5/19/05
			if ($groupnum < 10){
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
				print TXT "000$groupnum;";
			}
			elsif ($groupnum < 100){
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
				print TXT "00$groupnum;";
			}
			elsif($groupnum < 1000) {
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
				print TXT "0$groupnum;";
			}
			else {
				print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
				print TXT "$groupnum;";
			}
			$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
			if (defined($2)){
				print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
				print TXT "$2;";
			} 
			else {
				print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
				print TXT "$protein;";
			}

	#	printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
			print TXT "$$Protein_Hash{$protein}{'total_nomod'};";
			my $Seq = $$Protein_Hash{$protein}{'sequence'};
			for my $peptide (sort keys %{$$Protein_Hash{$protein}{'peptides'}}){
				my $seq = $peptide;
				$seq =~ s/[\*\#\@]//g;
				$Seq =~ s/($seq)/\L$1/ig;
			}
	  ##print "$protein\n$Seq\n";
			my $num = $Seq =~ s/([a-z])/$1/g;
	  ##printf "$num vs %d\n", length($Seq);exit;
			$$Protein_Hash{$protein}{'coverage'} = $num/length($Seq)*100;
			printf TXT "%.2f;", $num/length($Seq);

			printf HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Protein_Hash{$protein}{'occurrence'}</B></center></Font></TD>\n";

			print TXT "$$Protein_Hash{$protein}{'occurrence'};";

			for (my $i=0; $i<scalar(@$folders); $i++){
				my $msdak = $$folders[$i];
				if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
					print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
					print TXT "0;";
				} else {
					print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</B></center></Font></TD>\n";
					print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'};";
				}
			}

			printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
			print TXT "$$Protein_Hash{$protein}{'total_nomod'};";

		   for (my $i=0; $i<scalar(@$folders); $i++){
				my $msdak = $$folders[$i];
				if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
					print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
					print TXT "0;";
					print TXT "0;";
				}
				else
				{
					my $Seq = $$Protein_Hash{$protein}{'sequence'};
					for my $peptide (sort keys %{$$Sum_Hash{$msdak}{'protein'}{$protein}{'peptides'}}){
						my $seq = $peptide;
						$seq =~ s/[\*\#\@]//g;
						$Seq =~ s/($seq)/\L$1/ig;
					}
					
					printf HTML "$totalcgi\n", $protein, "$save_dir/fractions/sum_$$folders[$i]", $database,$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};
					print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};";
					my $num = $Seq =~ s/([a-z])/$1/g;
					$$Sum_Hash{$msdak}{'protein'}{$protein}{'coverage'} = $num/length($Seq)*100;
					printf TXT "%.2f;", $num/length($Seq);
				}
			}

			if (!defined($description)){
				print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
				print TXT ";";
			} else {
				$description =~ s/>/\&gt\;/g;
				$description =~ s/</\&lt\;/g;
				print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
				print TXT "$description;";
			}
			print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
			print TXT "$mass;";
			printf TXT "%d\n", length($Seq);
	###### added by yanji
			printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $tab_abundance_rounded;
			#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		}
	}
	
  	print HTML "</TABLE>\n";
  	print HTML "</BODY>\n";
  	print HTML "</HTML>\n";
  	close HTML;
	return ($unique_fpr,$subgroup_num);
}


sub gen_sumIDwGHtml{
	shift @_;
	
	my ($Sum_Hash, $Protein_Hash, $Peptide_Hash, $fprhash, $folders, $save_dir, $database, $grouphash) = @_;

	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";

	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";
	
	open (HTML, ">$save_dir/IDwG.html");		
	#open (HTML, ">$save_dir/IDwG.html");
	#print HTML "<HTML>\n";
	#print HTML "<HEAD>\n";
	#print HTML "<TITLE>Multiple ID Summary w/ Groups</TITLE>\n";
	open (TXT, ">$save_dir/gtext.txt");
	
	open (HTML2, ">$save_dir/IDwG_Redirect.html");
	#print_sumLegend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}

	$remove = $save_dir;

	$remove =~ s|/var/www/html||g;
	
	$redirect_path = "http://$server$remove/IDwG.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;		
	
	my $total_protein = scalar(keys %$Protein_Hash);
  	my $total_upeptide = 0;
  	for (sort {$$Peptide_Hash{$b}{'unique'} <=> $$Peptide_Hash{$a}{'unique'}} keys %$Peptide_Hash){
    		last if ($$Peptide_Hash{$_}{'unique'} == 0);
    		$total_upeptide++ if ($$Peptide_Hash{$_}{'unique'} == 1);
  	}
  	my $total_peptide = scalar(keys %$Peptide_Hash);
	my $total_SC = 0;
  	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
  	my $total_nomod = scalar(keys %peps_count);
	my $samples = scalar(@$folders);
	my $width = 80+($samples*10);
  	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=90%>\n";
  	
  	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";
  	#print HTML "<TR>\n";
	
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

        my $unique_fpr = 0;
        $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
        
    $uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'protein_fpr'};
    $pep_fdr = $$fprhash{'peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};

#    printf HTML "<TD Align=center><Font Size=2>Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};
    #print HTML "</TD>\n"; 
  	#print HTML "<TR>\n";
  	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";
	#print HTML "<CENTER><A HREF=ID.html>Hide Group Members</A></Center>";
	
	### these are the variables for the new html script
	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGILINK = "";
	
	
	$HTMLTABLEHEADER .= "<thead >\n";
	#print HTML "</TR></TD>\n";
	$HTMLTABLEHEADER .= "<tr>\n";
  	#print HTML "<TR><TD Align=center>";
  	
  	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  	#print HTML "</TD><TR>\n";
  	#print HTML "</TABLE>\n";
	
  	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = $width%>\n";
  	#print HTML "	<TR bgcolor=#EEEEEE>\n";
  	my $index_i = 0;
  	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
  	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";

	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
  	#print HTML "		<TH width = 5%%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
### changed on 8/9/2012

	$index_i = $index_i + 1; 
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
# print HTML "          <TH width = 5%%><B><Font Color = blue>Abundance</Font></B></TH>\n";
########################

	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Peptides</th>\n";
#  print HTML "		<TH width = 2%><B><Font Color = blue>Total<br><font size=1>Peptides</font></Font></B></TH>\n";
	print TXT "group;reference;totalpep;totalcoverage;";
	
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total SC</th>\n";
	#print HTML "          <TH width = 8%><B><Font Color = blue size=1>Total SC</Font></B></TH>\n";
	for (my $i=0; $i<scalar(@$folders); $i++){

		print TXT "$$folders[$i]_SC;";
		print TXT "$$folders[$i]_Coverage;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">SC $$folders[$i]</th>\n";
  		#print HTML "		<TH width = 5%><B><Font Color = blue size=1>SC $$folders[$i]</Font></B></TH>\n";
	}
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[$index_i] + '</a>');\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
        #print HTML "          <TH width = 8%><B><Font Color = blue size=1>Total Peptide</Font></B></TH>\n";
        
        for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_TP;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">TP $$folders[$i]</th>\n";
		#print HTML "            <TH width = 5%><B><Font Color = blue size=1>TP $$folders[$i]</Font></B></TH>\n";
		print TXT "$$folders[$i]_Coverage;";
	}
	
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
  	#print HTML "		<TH><B><Font Color = blue  size=1>Description</Font></B></TH>\n";
  	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
 	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
###### added by yanji
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
######
  	#print HTML "	</TR>\n";

	print TXT "Description;KD;Length\n";
	$HTMLTABLEHEADER .= "</tr>\n";
	$HTMLTABLEHEADER .= "</thead>\n";
	
 	my %PrevProHash;
	my %array_abundance;
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{

			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};

			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			$mass = sprintf("%3.0f", $mass/1000);
		
	########### added by yanji
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
	######	

			$HTMLTABLE .= "<tr class=\"gradeA\">\n";
		if($grouphash->{$protein}->{'order'}==1)
		{			 
			#print HTML "<TR BGColor = #ccffcc>\n";
		}
		else
		{
			#print HTML "<TR>\n";
		}
		#added if loop DMD 5/19/05
			if ($groupnum < 10){
				$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
				print TXT "000$groupnum;";
			}
			elsif ($groupnum < 100){
				$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
				print TXT "00$groupnum;";
			}
			elsif($groupnum < 1000) {
				$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
				print TXT "0$groupnum;";
		}
		 else {
		 	$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
			print TXT "$groupnum;";
		}
		$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
		if (defined($2)){
			$HTMLTABLE .= "<td align=\"center\">$2</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
			print TXT "$2;";
		} 
		else {
			$HTMLTABLE .= "<td align=\"center\">$protein</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
			print TXT "$protein;";
		}

		print TXT "$$Protein_Hash{$protein}{'total_nomod'};";
		my $Seq = $$Protein_Hash{$protein}{'sequence'};
		for my $peptide (sort keys %{$$Protein_Hash{$protein}{'peptides'}}){
				my $seq = $peptide;
				$seq =~ s/[\*\#\@]//g;
				$Seq =~ s/($seq)/\L$1/ig;
		}
		my $num = $Seq =~ s/([a-z])/$1/g;
		$$Protein_Hash{$protein}{'coverage'} = $num/length($Seq)*100;
		printf TXT "%.2f;", $num/length($Seq);

			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'occurrence'}</td>\n";
			#printf HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Protein_Hash{$protein}{'occurrence'}</B></center></Font></TD>\n";

			print TXT "$$Protein_Hash{$protein}{'occurrence'};";

		for (my $i=0; $i<scalar(@$folders); $i++){
			my $msdak = $$folders[$i];
			if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
				$HTMLTABLE .= "<td align=\"center\">0</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
				print TXT "0;";
			} else {
				$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</B></center></Font></TD>\n";
				print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'};";
			}
		}

		$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'total_nomod'}</td>\n";
		#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
		print TXT "$$Protein_Hash{$protein}{'total_nomod'};";

		   for (my $i=0; $i<scalar(@$folders); $i++){
			my $msdak = $$folders[$i];
					if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
				
				$HTMLTABLE .= "<td align=\"center\">0</td>\n";		
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
							print TXT "0;";
							print TXT "0;";
					}
			else
			{
				my $Seq = $$Protein_Hash{$protein}{'sequence'};
					for my $peptide (sort keys %{$$Sum_Hash{$msdak}{'protein'}{$protein}{'peptides'}}){
						my $seq = $peptide;
						$seq =~ s/[\*\#\@]//g;
						$Seq =~ s/($seq)/\L$1/ig;
					}
					$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'}</td>\n";
				#printf HTML "$totalcgi\n", $protein, "$save_dir/fractions/sum_$$folders[$i]", $database,$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};					
				print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};";
					my $num = $Seq =~ s/([a-z])/$1/g;
				$$Sum_Hash{$msdak}{'protein'}{$protein}{'coverage'} = $num/length($Seq)*100;
				printf TXT "%.2f;", $num/length($Seq);
			}
		}

		if (!defined($description)){
			$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
			#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
			print TXT ";";
		} else {
			$description =~ s/>/\&gt\;/g;
			$description =~ s/</\&lt\;/g;
			$HTMLTABLE .= "<td align=\"center\">$description</td>\n";
			#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
			print TXT "$description;";
		}
		$HTMLTABLE .= "<td align=\"center\">$mass</td>\n";
		#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
		print TXT "$mass;";
		printf TXT "%d\n", length($Seq);
	###### added by yanji
		$HTMLTABLE .= "<td align=\"center\">$tab_abundance_rounded</td>\n";
		#$HTMLTABLE .= "<td align=\"center\">$abundance</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		$HTMLTABLE .= "</tr>\n";
	}

	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "Multiple ID Summary w/ Groups", "ID.html", "Hide Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
	close TXT;
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	#close HTML;
}



sub gen_sumIDmodHtml{
	shift @_;
	
	my ($Sum_Hash, $Protein_Hash, $Peptide_Hash, $fprhash, $folders, $save_dir, $database, $grouphash) = @_;

	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";

	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	open (HTML, ">$save_dir/IDmod.html");
	#open (HTML, ">$save_dir/IDmod.html");
	#print HTML "<HTML>\n";
	#print HTML "<HEAD>\n";
	#print HTML "<TITLE>Multiple ID Summary without Groups for Modifications</TITLE>\n";
	open (TXT, ">$save_dir/modtext.txt");
	open (HTML2, ">$save_dir/IDmod_Redirect.html");
	#print_sumLegend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}	
	
	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	
	$redirect_path = "http://$server$remove/IDmod.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;		
	
	my $total_protein = scalar(keys %$Protein_Hash);
  	my $total_upeptide = 0;
  	for (sort {$$Peptide_Hash{$b}{'unique'} <=> $$Peptide_Hash{$a}{'unique'}} keys %$Peptide_Hash){
    		last if ($$Peptide_Hash{$_}{'unique'} == 0);
    		$total_upeptide++ if ($$Peptide_Hash{$_}{'unique'} == 1);
  	}
  	my $total_peptide = scalar(keys %$Peptide_Hash);
	my $total_SC = 0;
  	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
  	my $total_nomod = scalar(keys %peps_count);
	my $samples = scalar(@$folders);
	my $width = 70+($samples*10);
	
	
  	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=80%>\n";
  	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";
  	#print HTML "<TR>\n";
  	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

    my $unique_fpr = 0;
    $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
       
    $uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'final_protein_fpr'};
    $pep_fdr = $$fprhash{'final_peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};

#    printf HTML "<TD Align=center><Font Size=2>Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};
    #print HTML "</TD>\n"; 
  	#print HTML "<TR>\n";
  	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";
	
	#print HTML "<CENTER><A HREF=IDwGmod.html>Show Group Members</A></Center>";
	

	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGILINK = "";

	
	
	#print HTML "</TR></TD>\n";
  	#print HTML "<TR><TD Align=center>";
  	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  	#print HTML "</TD><TR>\n";
  	#print HTML "</TABLE>\n";

  	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = $width%>\n";
  	#print HTML "	<TR bgcolor=#EEEEEE>\n";
  	$HTMLTABLEHEADER .= "<thead >\n";
	$HTMLTABLEHEADER .= "<tr>\n";
	my $index_i = 0;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";

  	
  	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
  	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
  	#print HTML "		<TH width = 5%%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
### changed on 8/9/2012 

	
	print TXT "group;reference;totalpep;totalcoverage;";
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total SC</th>\n";
	#print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total SC</Font></B></TH>\n";
	for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_SC;";
		print TXT "$$folders[$i]_Coverage;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">SC $$folders[$i]</th>\n";
  		#print HTML "		<TH width = 4%><B><Font Color = blue size=1>SC $$folders[$i]</Font></B></TH>\n";
	}

	$index_i = $index_i + 1;
	$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[$index_i] + '</a>');\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
    #print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total Peptide</Font></B></TH>\n";
    for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_TP;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">TP $$folders[$i]</th>\n";
		#print HTML "            <TH width = 4%><B><Font Color = blue size=1>TP $$folders[$i]</Font></B></TH>\n";
		print TXT "$$folders[$i]_Coverage;";
	}
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";	
  	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
  	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
 	#print HTML "		<TH width = 7%><B><Font Color = blue size=1>Mass (KD)</Font></B></TH>\n";
 	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
  	#print HTML "	</TR>\n";

	print TXT "Description;KD;Length\n";
	$HTMLTABLEHEADER .= "</tr>\n";
	$HTMLTABLEHEADER .= "</thead>\n";
	
	
 	my %PrevProHash;
	my %array_abundance;
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{
		if($grouphash->{$protein}->{'order'}==1)
		{
			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};

			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			$mass = sprintf("%3.0f", $mass/1000);
		
	########### added by yanji
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;
	######	

			$HTMLTABLE .= "<tr class=\"gradeA\">\n";
			#print HTML "<TR BGColor = #ccffcc>\n";
		#added if loop DMD 5/19/05
			if ($groupnum < 10){
				$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
				print TXT "000$groupnum;";
			}
			elsif ($groupnum < 100){
				$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
				print TXT "00$groupnum;";
			}
			elsif($groupnum < 1000) {
				$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
				print TXT "0$groupnum;";
		}
		 else {
		 	$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
			print TXT "$groupnum;";
		}
		$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
		if (defined($2)){
			$HTMLTABLE .= "<td align=\"center\">$2</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
			print TXT "$2;";
		} 
		else {
			$HTMLTABLE .= "<td align=\"center\">$protein</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
			print TXT "$protein;";
		}

		print TXT "$$Protein_Hash{$protein}{'total_nomod'};";
		my $Seq = $$Protein_Hash{$protein}{'sequence'};
		for my $peptide (sort keys %{$$Protein_Hash{$protein}{'peptides'}}){
				my $seq = $peptide;
				$seq =~ s/[\*\#\@]//g;
				$Seq =~ s/($seq)/\L$1/ig;
		}
		my $num = $Seq =~ s/([a-z])/$1/g;
		$$Protein_Hash{$protein}{'coverage'} = $num/length($Seq)*100;
		printf TXT "%.2f;", $num/length($Seq);

		$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'occurrence'}</td>\n";
			#printf HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Protein_Hash{$protein}{'occurrence'}</B></center></Font></TD>\n";

			print TXT "$$Protein_Hash{$protein}{'occurrence'};";

		for (my $i=0; $i<scalar(@$folders); $i++){
			my $msdak = $$folders[$i];
			if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
				$HTMLTABLE .= "<td align=\"center\">0</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
				print TXT "0;";
			} else {
				$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</B></center></Font></TD>\n";
				
				print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'};";
			}
		}

		#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
		print TXT "$$Protein_Hash{$protein}{'total_nomod'};";

		   for (my $i=0; $i<scalar(@$folders); $i++){
			my $msdak = $$folders[$i];
					if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
						$HTMLTABLE .= "<td align=\"center\">0</td>\n";		
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
							print TXT "0;";
							print TXT "0;";
					}
			else
			{
				my $Seq = $$Protein_Hash{$protein}{'sequence'};
					for my $peptide (sort keys %{$$Sum_Hash{$msdak}{'protein'}{$protein}{'peptides'}}){
						my $seq = $peptide;
						$seq =~ s/[\*\#\@]//g;
						$Seq =~ s/($seq)/\L$1/ig;
					}
				$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'}</td>\n";
				#printf HTML "$totalcgi\n", $protein, "$save_dir/fractions/sum_$$folders[$i]_mod", $database,$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};	
				$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'}</td>\n";
				#printf HTML "$uniquecgi\n", $protein, $save_dir, $$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};
				print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};";
					my $num = $Seq =~ s/([a-z])/$1/g;
				$$Sum_Hash{$msdak}{'protein'}{$protein}{'coverage'} = $num/length($Seq)*100;
				printf TXT "%.2f;", $num/length($Seq);
			}
		}

		if (!defined($description)){
			$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
			#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
			print TXT ";";
		} else {
			$description =~ s/>/\&gt\;/g;
			$description =~ s/</\&lt\;/g;
			$HTMLTABLE .= "<td align=\"center\">$description</td>\n";
			#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
			print TXT "$description;";
		}
		$HTMLTABLE .= "<td align=\"center\">$mass</td>\n";
		#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
		print TXT "$mass;";
		printf TXT "%d\n", length($Seq);
	###### added by yanji
		$HTMLTABLE .= "<td align=\"center\">$tab_abundance_rounded</td>\n";
		#$HTMLTABLE .= "<td align=\"center\">$abundance</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
	}
	}
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	#close HTML;
	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "Multiple ID Summary without Groups for Modifications", "IDwGmod.html", "Show Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
   	close TXT;
}


sub gen_sumIDwGmodHtml{
	shift @_;
	
	my ($Sum_Hash, $Protein_Hash, $Peptide_Hash, $fprhash, $folders, $save_dir, $database, $grouphash) = @_;


	my $num_group = "";
	my $num_unique_pro = "";
	my $num_pro = "";
	my $num_pep = "";
	my $num_SC = "";
	my $uniqpro_fdr = "";
	my $pro_fdr = "";
	my $pep_fdr = "";
	my $num_CM_mod = "";


	my $utils = idsum2::CommonUtils->new();
	my $debug = 0;
	my $printoutfile = 0;

	my $server = $utils->server();
	#my $http = "HREF=http://$server/cgi-bin";
	my $http = "http://$server/cgi-bin";

	my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
	my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
	my $color = "#33DD33";

	open (HTML, ">$save_dir/IDwGmod.html");
	#open (HTML, ">$save_dir/IDwGmod.html");
	#print HTML "<HTML>\n";
	#print HTML "<HEAD>\n";
	#print HTML "<TITLE>Multiple ID Summary with Groups for Modifications</TITLE>\n";
	open (TXT, ">$save_dir/modwgtext.txt");
	open (HTML2, ">$save_dir/IDwGmod_Redirect.html");	
	#print_sumLegend(*HTML);
	#print HTML "</HEAD>\n\n";

	my $current_user = qx[whoami];
	chomp($current_user);
	if ($save_dir =~ /^\/home/) 
	{
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} 
	else 
	{
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}

	$remove = $save_dir;
	$remove =~ s|/var/www/html||g;
	
	$redirect_path = "http://$server/$current_user/$remove/IDwGmod.html";
	
	print HTML2 "<meta http-equiv=\"refresh\" content=\"0; url=$redirect_path\" />\n";
	close(HTML2);
	
	
	my ($groupnum, $subgroup_num);
	$groupnum = $grouphash->{(reverse sort {$grouphash->{$a}->{'group'} <=> $grouphash->{$b}->{'group'}} keys %$grouphash)[0]}{'group'};

	map {$subgroup_num++ if($grouphash->{$_}->{'order'}==1)} keys %$grouphash;		
	
	my $total_protein = scalar(keys %$Protein_Hash);
  	my $total_upeptide = 0;
  	for (sort {$$Peptide_Hash{$b}{'unique'} <=> $$Peptide_Hash{$a}{'unique'}} keys %$Peptide_Hash){
    		last if ($$Peptide_Hash{$_}{'unique'} == 0);
    		$total_upeptide++ if ($$Peptide_Hash{$_}{'unique'} == 1);
  	}
  	my $total_peptide = scalar(keys %$Peptide_Hash);
	my $total_SC = 0;
  	for (values %$Peptide_Hash){
		$total_SC += scalar(keys %{$$_{'outfiles'}});
	}
	my %peps_count;
	for my $peps (keys %$Peptide_Hash){ # added DMD 5/19/05
		$peps =~ s/C[\*\@\#]/C/g;
		$peps =~ s/M[\*\@\#]/M/g;
		$peps_count{$peps} = $peps if (!defined($peps_count{$peps}));
	}
  	my $total_nomod = scalar(keys %peps_count);
	my $samples = scalar(@$folders);
	my $width = 70+($samples*10);
  	#print HTML "<BODY>\n";
	#print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=80%>\n";
  	#print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";
  	#print HTML "<TR>\n";
	$num_group = $groupnum;
	#print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	$num_unique_pro = $subgroup_num;
	#print HTML "Unique proteins = $subgroup_num &nbsp&nbsp&nbsp ";
	$num_pro = $total_protein;
	#print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	$num_pep = $total_peptide;
	#print HTML "Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	$num_SC = $total_SC;
	#print HTML "SC = $total_SC</font></TD>\n"; #changed to total peptide DMD 5/19/05
	#print HTML " </TR>\n";

        my $unique_fpr = 0;
        $unique_fpr = $$fprhash{'protein_bad'} / $subgroup_num * 100;
    
    $uniqpro_fdr = $unique_fpr;
    $pro_fdr = $$fprhash{'final_protein_fpr'};
    $pep_fdr = $$fprhash{'final_peptide_fpr'};
    #printf HTML "<TD Align=center><Font Size=2>Unique Protein FDR = %.2f%%,&nbsp Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $unique_fpr, $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};

    #print HTML "</TD>\n"; 
  	#print HTML "<TR>\n";
	$num_CM_mod = $total_nomod;
	#print HTML "<TD Align=center><Font Size=2>Peptides without C and M modifications = $total_nomod</TD>\n";
	#print HTML " </TR>\n";
	#print HTML "<TR><TD align=center>\n";
	
	#print HTML "<CENTER><A HREF=IDmod.html>Hide Group Members</A></Center>";
	

	#print HTML "</TR></TD>\n";
  	#print HTML "<TR><TD Align=center>";
  	#print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  	#print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  	#print HTML "</TD><TR>\n";
  	#print HTML "</TABLE>\n";

	my $HTMLTABLE = "";
	my $HTMLTABLEHEADER = "";
	my $CGISPECIAL = "";

  	#print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = $width%>\n";
  	#print HTML "	<TR bgcolor=#EEEEEE>\n";
  	
  	
	$HTMLTABLEHEADER .= "<thead >\n";
	$HTMLTABLEHEADER .= "<tr>\n";
	my $index_i = 0;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Group</th>\n";
  	
  	#print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
  	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Reference</th>\n";
  	#print HTML "		<TH width = 5%%><B><Font Color = blue size=1>Reference</Font></B></TH>\n";
	print TXT "group;reference;totalpep;totalcoverage;";
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total SC</th>\n";
	#print HTML "          <TH width = 7%><B><Font Color = blue size=1>Total SC</Font></B></TH>\n";
	for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_SC;";
		print TXT "$$folders[$i]_Coverage;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_unique_beta.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">SC $$folders[$i]</th>\n";
  		#print HTML "		<TH width = 4%><B><Font Color = blue size=1>SC $$folders[$i]</Font></B></TH>\n";
	}
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_total.cgi?protein=' + aData[1] + '&dir=$save_dir&database=$database\">' + aData[$index_i] + '</a>');\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Total Peptide</th>\n";
        #print HTML "          <TH width = 7%><B><Font Color = blue  size=1>Total Peptide</Font></B></TH>\n";
        for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_TP;";
		$index_i = $index_i + 1;
		$CGILINK .= "\$('td:eq($index_i)', nRow).html('<a href=\"$http/id_shared.cgi?protein=' + aData[1] + '&dir=$save_dir\">' + aData[$index_i] + '</a>')\n";
		$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">TP $$folders[$i]</th>\n";
		#print HTML "            <TH width = 4%><B><Font Color = blue  size=1>TP $$folders[$i]</Font></B></TH>\n";
		print TXT "$$folders[$i]_Coverage;";
	}
	$index_i = $index_i + 1;
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Description</th>\n";
  	#print HTML "		<TH><B><Font Color = blue size=1>Description</Font></B></TH>\n";
  	$index_i = $index_i + 1;
 	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Mass (KD)</th>\n";
 	#print HTML "		<TH width = 7%><B><Font Color = blue  size=1>Mass (KD)</Font></B></TH>\n";
	$HTMLTABLEHEADER .= "<th style=\"border-bottom:3px solid black\">Abundance</th>\n";
	#print HTML "          <TH width = 2%><B><Font Color = blue size=1>Abundance</Font></B></TH>\n";
  	#print HTML "	</TR>\n";
	$HTMLTABLEHEADER .= "</tr>\n";
	$HTMLTABLEHEADER .= "</thead>\n";
	
	
	print TXT "Description;KD;Length\n";

	
 	my %PrevProHash;
	my %array_abundance;
	foreach my $protein (sort {$grouphash->{$a}->{'group'}<=>$grouphash->{$b}->{'group'} ||
							   $grouphash->{$a}->{'subgroup'}<=>$grouphash->{$b}->{'subgroup'} ||
							   $grouphash->{$a}->{'order'}<=>$grouphash->{$b}->{'order'}
						} keys %$grouphash)
	{

			my $groupnum = $grouphash->{$protein}->{'group'}."\.".$grouphash->{$protein}->{'subgroup'};

			my $description = $$Protein_Hash{$protein}{'annotation'};
			my $mass = $$Protein_Hash{$protein}{'MW'};
			$mass = sprintf("%3.0f", $mass/1000);
		
			my $abundance = $$Protein_Hash{$protein}{'abundance'};
			my $tab_abundance_rounded = sprintf "%0.2f", $abundance;

		$HTMLTABLE .= "<tr class=\"gradeA\">\n";
		if($grouphash->{$protein}->{'order'}==1)
		{			 
			#print HTML "<TR BGColor = #ccffcc>\n";
		}
		else
		{
			#print HTML "<TR>\n";
		}
		#added if loop DMD 5/19/05
			if ($groupnum < 10){
				$HTMLTABLE .= "<td align=\"center\">000$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
				print TXT "000$groupnum;";
			}
			elsif ($groupnum < 100){
				$HTMLTABLE .= "<td align=\"center\">00$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
				print TXT "00$groupnum;";
			}
			elsif($groupnum < 1000) {
				$HTMLTABLE .= "<td align=\"center\">0$groupnum</td>\n";
				#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
				print TXT "0$groupnum;";
		}
		 else {
		 	$HTMLTABLE .= "<td align=\"center\">$groupnum</td>\n";
			#print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
			print TXT "$groupnum;";
		}
		$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
		if (defined($2)){
			$HTMLTABLE .= "<td align=\"center\">$2</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
			print TXT "$2;";
		} 
		else {
			$HTMLTABLE .= "<td align=\"center\">$protein</td>\n";
			#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
			print TXT "$protein;";
		}

		print TXT "$$Protein_Hash{$protein}{'total_nomod'};";
		my $Seq = $$Protein_Hash{$protein}{'sequence'};
		for my $peptide (sort keys %{$$Protein_Hash{$protein}{'peptides'}}){
				my $seq = $peptide;
				$seq =~ s/[\*\#\@]//g;
				$Seq =~ s/($seq)/\L$1/ig;
		}
		my $num = $Seq =~ s/([a-z])/$1/g;
		$$Protein_Hash{$protein}{'coverage'} = $num/length($Seq)*100;
		printf TXT "%.2f;", $num/length($Seq);

			$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'occurrence'}</td>\n";
			#printf HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Protein_Hash{$protein}{'occurrence'}</B></center></Font></TD>\n";

			print TXT "$$Protein_Hash{$protein}{'occurrence'};";

		for (my $i=0; $i<scalar(@$folders); $i++){
			my $msdak = $$folders[$i];
			if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
				$HTMLTABLE .= "<td align=\"center\">0</td>\n";
				#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
				print TXT "0;";
			} else {
				$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</td>\n";
				#			print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'}</B></center></Font></TD>\n";
				print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'occurrence'};";
			}
		}

		$HTMLTABLE .= "<td align=\"center\">$$Protein_Hash{$protein}{'total_nomod'}</td>\n";
		#printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
		print TXT "$$Protein_Hash{$protein}{'total_nomod'};";

		   for (my $i=0; $i<scalar(@$folders); $i++){
			my $msdak = $$folders[$i];
					if (!defined($$Sum_Hash{$msdak}{'protein'}{$protein})){
						$HTMLTABLE .= "<td align=\"center\">0</td>\n";
						#print HTML "<TD><Font Size = 2.5 Color = blue><center><B>0</B></center></Font></TD>\n";
						print TXT "0;";
						print TXT "0;";
					}
			else
			{
				my $Seq = $$Protein_Hash{$protein}{'sequence'};
					for my $peptide (sort keys %{$$Sum_Hash{$msdak}{'protein'}{$protein}{'peptides'}}){
						my $seq = $peptide;
						$seq =~ s/[\*\#\@]//g;
						$Seq =~ s/($seq)/\L$1/ig;
					}
					$HTMLTABLE .= "<td align=\"center\">$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'}</td>\n";
				#printf HTML "$totalcgi\n", $protein, "$save_dir/fractions/sum_$$folders[$i]_mod", $database,$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};	
				#printf HTML "$uniquecgi\n", $protein, $save_dir, $$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};
				print TXT "$$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};";
					my $num = $Seq =~ s/([a-z])/$1/g;
				$$Sum_Hash{$msdak}{'protein'}{$protein}{'coverage'} = $num/length($Seq)*100;
				printf TXT "%.2f;", $num/length($Seq);
			}
		}

		if (!defined($description)){
			$HTMLTABLE .= "<td align=\"center\">&nbsp</td>\n";
			#print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
			print TXT ";";
		} else {
			$description =~ s/>/\&gt\;/g;
			$description =~ s/</\&lt\;/g;
			$HTMLTABLE .= "<td align=\"center\">$description</td>\n";
			#print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
			print TXT "$description;";
		}
		$HTMLTABLE .= "<td align=\"center\">$mass</td>\n";
		#print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
		print TXT "$mass;";
		printf TXT "%d\n", length($Seq);
	###### added by yanji
		$HTMLTABLE .= "<td align=\"center\">$tab_abundance_rounded</td>\n";
		#$HTMLTABLE .= "<td align=\"center\">$abundance</td>\n";
		#printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
		$HTMLTABLE .= "</tr>\n";
	}

  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	#close HTML;
  	
	my $final_html = generateHTMLStructure($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, "Multiple ID Summary with Groups for Modifications", "IDmod.html", "Hide Group Members");
	
  	#print HTML "</TABLE>\n";
  	#print HTML "</BODY>\n";
  	#print HTML "</HTML>\n";
  	
   	print HTML $final_html;
   	close HTML;
	close TXT;
}

sub generateHTMLStructure {
	#shift @_;
	my ($num_group, $num_unique_pro, $num_pro, $num_pep, $num_SC, $uniqpro_fdr, $pro_fdr, $pep_fdr, $num_CM_mod, $HTMLTABLEHEADER, $HTMLTABLE, $CGILINK, $TITLE, $OtherLink, $OtherLinkTITLE) = @_;
	
	my $html_page = "";
	$html_page .= "<!DOCTYPE html>\n";
	$html_page .= "<html lang=\"en\">\n";
  	$html_page .= "<head>\n";
    $html_page .= "<meta charset=\"utf-8\">\n";
    $html_page .= "<title>Protein Identification Summary</title>\n";
    $html_page .= "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
    $html_page .= "<meta name=\"description\" content=\"\">\n";
    $html_page .= "<meta name=\"author\" content=\"\">\n";

    $html_page .= "<!-- Le styles -->\n";
    #$html_page .= "<link href=\"css/bootstrap.css\" rel=\"stylesheet\">\n";
    $html_page .= "<link href=\"/webconfig/css/bootstrap.css\" rel=\"stylesheet\">\n";
    #$html_page .= "<link href=\"css/bootstrap-responsive.css\" rel=\"stylesheet\">\n";
    $html_page .= "<link href=\"/webconfig/css/bootstrap-responsive.css\" rel=\"stylesheet\">\n";
   
    #$html_page .= "<link href=\"dataTable/jquery.dataTables.css\" rel=\"stylesheet\">\n";
    $html_page .= "<link href=\"/webconfig/dataTable/jquery.dataTables.css\" rel=\"stylesheet\">\n";
    $html_page .= "<style type=\"text/css\">\n";
      $html_page .= "body {\n";
        $html_page .= "padding-top: 60px;\n";
        $html_page .= "padding-bottom: 40px;\n";
      $html_page .= "}\n";
      $html_page .= ".sidebar-nav {\n";
        $html_page .= "padding: 9px 0;\n";
      $html_page .= "}\n";

	#$html_page .= ".dataTables_filter {\n";
	#$html_page .= "display: none; \n";
	#$html_page .= "}\n";

	$html_page .= "</style>\n";

    $html_page .= "<!-- HTML5 shim, for IE6-8 support of HTML5 elements -->\n";
    $html_page .= "<!--[if lt IE 9]>\n";
      $html_page .= "<script src=\"http://html5shim.googlecode.com/svn/trunk/html5.js\"></script>\n";
    $html_page .= "<![endif]-->\n";

    $html_page .= "<!-- Fav and touch icons -->\n";
    $html_page .= "<link rel=\"apple-touch-icon-precomposed\" sizes=\"144x144\" href=\"ico/apple-touch-icon-144-precomposed.png\">\n";
    $html_page .= "<link rel=\"apple-touch-icon-precomposed\" sizes=\"114x114\" href=\"ico/apple-touch-icon-114-precomposed.png\">\n";
      $html_page .= "<link rel=\"apple-touch-icon-precomposed\" sizes=\"72x72\" href=\"ico/apple-touch-icon-72-precomposed.png\">\n";
                   $html_page .=  "<link rel=\"apple-touch-icon-precomposed\" href=\"ico/apple-touch-icon-57-precomposed.png\">\n";
                     $html_page .=  "<link rel=\"shortcut icon\" href=\"ico/favicon.png\">\n";
  	$html_page .= "</head>\n";


  	$html_page .= "<body>\n";

    $html_page .= "<div class=\"navbar navbar-inverse navbar-fixed-top\">\n";
      $html_page .= "<div class=\"navbar-inner\">\n";
        $html_page .= "<div class=\"container-fluid\">\n";
          $html_page .= "<a class=\"btn btn-navbar\" data-toggle=\"collapse\" data-target=\".nav-collapse\">\n";
            $html_page .= "<span class=\"icon-bar\"></span>\n";
            $html_page .= "<span class=\"icon-bar\"></span>\n";
            $html_page .= "<span class=\"icon-bar\"></span>\n";
          $html_page .= "</a>\n";
          $html_page .= "<a class=\"brand\" href=\"http://10.4.1.28\">Junmin Peng Lab&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</a>\n";
          $html_page .= "<div class=\"nav-collapse collapse\">\n";
            $html_page .= "<p class=\"navbar-text pull-right\">\n";
              $html_page .= "<!--Logged in as <a href=\"#\" class=\"navbar-link\">Username</a>-->\n";
            $html_page .= "</p>\n";
            $html_page .= "<ul class=\"nav\">\n";
              $html_page .= "<li><a href=\"http://10.4.1.28\">Home</a></li>\n";
              $html_page .= "<li class=\"active\"><a href=\"$OtherLink\">$OtherLinkTITLE</a></li>\n";
              #$html_page .= "<li><a href=\"Help.php\">Help</a></li>\n";
            $html_page .= "</ul>\n";
          $html_page .= "</div><!--/.nav-collapse -->\n";
        $html_page .= "</div>\n";
      $html_page .= "</div>\n";
    $html_page .= "</div>\n";

    $html_page .= "<div class=\"container-fluid\">\n";
      $html_page .= "<div class=\"row-fluid\">\n";
        $html_page .= "<!--<div class=\"span3\">\n";
          $html_page .= "<div class=\"well sidebar-nav\">\n";
            $html_page .= "<ul class=\"nav nav-list\">\n";
              $html_page .= "<li class=\"nav-header\">Sidebar</li>\n";
              $html_page .= "<li class=\"active\"><a href=\"index.php\">Home</a></li>\n";
              #$html_page .= "<li><a href=\"Help.php\">Help</a></li>\n";
            $html_page .= "</ul>\n";
          $html_page .= "</div>\n";
        $html_page .= "</div>-->\n";
        $html_page .= "<div class=\"span12\">\n";
          $html_page .= "<div class=\"hero-unit\">\n";
            $html_page .= "<h1>$TITLE</h1> <br />\n";
	    $html_page .= "<p>Groups = " . $num_group . "  Unique proteins = " . $num_unique_pro . "  Proteins = " . $num_pro . "  Peptides = " . $num_pep . "   SC = " . $num_SC . "</p><br />\n";
	    $html_page .= "<p>Unique Protein FDR = " . $uniqpro_fdr . "%,  Protein FDR = " . $pro_fdr . "%,  Peptide FDR = " . $pep_fdr . "% </p> <br />\n";
	    $html_page .= "<p>Peptides without C and M modifications = " . $num_CM_mod . "</p> <br />\n";
	    $html_page .= "<p><A HREF=$OtherLink>$OtherLinkTITLE</A></p><br />\n";
            $html_page .= "<h2>Table</h2>\n";
            $html_page .= "<div id=\"demo\">\n";
		$html_page .= "<h4 small>\n";
		$html_page .= "<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" font=10 class=\"display\" id=\"example\" width=\"100%\">\n";
        $html_page .= "\n";
        $html_page .= $HTMLTABLEHEADER . "\n";
        
        $html_page .= "<tbody>\n";
		
		$html_page .= $HTMLTABLE . "\n";
        
        $html_page .= "</tbody>\n";
        
        #$html_page .= "<tfoot>\n";
        #       $html_page .= "<tr>\n";
        #                $html_page .= "<th>Group</th>\n";
        #                $html_page .= "<th>Reference</th>\n";
        #                $html_page .= "<th>Total SC</th>\n";
        #                $html_page .= "<th>SC tmt1</th>\n";
        #                $html_page .= "<th>SC tmt2</th>\n";
        #                $html_page .= "<th>Total Peptide</th>\n";
		#	$html_page .= "<th>TP tmt1</th>\n";
        #                $html_page .= "<th>TP tmt2</th>\n";
		#	$html_page .= "<th>Description</th>\n";
		#	$html_page .= "<th>Mass (KD)</th>\n";
        #                $html_page .= "<th>Abundance</th>\n";
			

         #       $html_page .= "</tr>\n";
        #$html_page .= "</tfoot>\n";
		$html_page .= "</table>\n";
		$html_page .= "</h4 small>\n";

            
          $html_page .= "</div>\n";
        $html_page .= "</div><!--/span-->\n";
      $html_page .= "</div><!--/row-->\n";

      $html_page .= "<hr>\n";

      $html_page .= "<footer>\n";
        $html_page .= "<p>For encountered problems and questions please contact Timothy Shaw tim.shaw\@stjude.org</p>\n";
        $html_page .= "<p>&copy; St Jude Research Hospital 2014</p>\n";
      $html_page .= "</footer>\n";

    $html_page .= "</div><!--/.fluid-container-->\n";

    $html_page .= "<!-- Le javascript\n";
    $html_page .= "================================================== -->\n";
    $html_page .= "<!-- Placed at the end of the document so the pages load faster -->\n";
    $html_page .= "<script src=\"/webconfig/media/js/jquery.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/media/js/jquery.data.Tables.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-transition.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-alert.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-modal.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-dropdown.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-scrollspy.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-tab.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-tooltip.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-popover.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-button.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-collapse.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-carousel.js\"></script>\n";
    $html_page .= "<script src=\"/webconfig/js/bootstrap-typeahead.js\"></script>\n";

    
    $html_page .= "<script type=\"text/javascript\" language=\"javascript\" src=\"/webconfig/dataTable/jquery.dataTables.js\"></script>\n";
     $html_page .= "<script type=\"text/javascript\" charset=\"utf-8\">\n";
                        $html_page .= "\$(document).ready(function() {\n";
                                $html_page .= "\$('#example').dataTable({\n";
$html_page .= "\"aLengthMenu\": [\n";
            $html_page .= "[10, 50, 100, 200, 1000, -1],\n";
            $html_page .= "[10, 50, 100, 200, 1000, \"All\"]\n";
        $html_page .= "],\n";
        $html_page .= "\"sDom\": '<\"wrapper\"lfptip>',\n";
     $html_page .= "\"aaSorting\": [[ 0, \"asc\" ]],\n";
     $html_page .= "\"fnRowCallback\": function( nRow, aData, iDisplayIndex ) {\n";

			# NEED TO DEFINE THE HTML LINKS HERE!!!
            $html_page .= $CGILINK;
          	#$html_page .= "\$('td:eq(5)', nRow).html('<a href=\"http://spiderscluster.stjude.org/cgi-bin/id_total.cgi?protein=' + aData[1] + '&dir=/var/www/html/jmertz/2014/Mib1_Interactome_2014/sum_all&database=/home/jmertz/2014/Mib1_Interactome_2014/database/rat140320_ft_mc2_c57TMT.fasta\">' +\n";
               #$html_page .= "aData[5] + '</a>');\n";

    $html_page .= "return nRow;\n";
    $html_page .= "},\n";

	$html_page .= "});\n";
    $html_page .= "} );\n";


 	$html_page .= "// Assigning BootStrap CSS classes to the DataTables components\n";
	
 	$html_page .= "// Adding arrows to table headers used to show sorting direction\n";
 	$html_page .= "\$('th').each(function(){ \n";
 	$html_page .= "if(\$(this).text() != ''){\n";
 	$html_page .= "\$(this).append('<i class=\"icon-chevron-up\"></i> <i class=\"icon-chevron-down\"></i>');\n";
 	$html_page .= "}\n";
 	$html_page .= "});\n";

    $html_page .= "</script>\n";
 	$html_page .= "</body>\n";
	$html_page .= "</html>\n";


	while ($html_page =~ /<tr class="gradeA">\n<tr class="gradeA">/) {
		$html_page =~ s/<tr class="gradeA">\n<tr class="gradeA">/<tr class="gradeA">/g;
	}

	return $html_page;
}
1;



