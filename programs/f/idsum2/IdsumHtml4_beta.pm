#/usr/bin/perl -wT

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: idsum2::IdsumHtml4_beta

package idsum2::IdsumHtml4_beta;
use strict;
use Cwd;
use Storable;
use Data::Dumper;
use idsum2::CommonUtils;

my $utils = idsum2::CommonUtils->new();
my $debug = 0;
my $printoutfile = 0;
my $outpath = getcwd();
my $server = $utils->server();
my $http = "HREF=http://$server/cgi-bin";

my $totalcgi = "<TD><CENTER><A TARGET = \"win1\" $http/id_total.cgi?protein=%s&dir=%s&database=%s>%d</A></CENTER></TD>";
my $uniquecgi = "<TD><CENTER><A TARGET = \"win2\" $http/id_unique_beta.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
my $sharedcgi = "<TD><CENTER><A TARGET = \"win3\" $http/id_shared.cgi?protein=%s&dir=%s>%d</A></CENTER></TD>";
my $color = "#33DD33";

sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);
  
  return $self;
}

#-----------------------------------------------------------------

sub group_proteins{
	shift @_;
	my ($protein_hash, $peptide_hash, $percent, $param_hash) = @_;
	my $bypass = $$param_hash{'bypass_grouping'};
	my $bypassfil = $$param_hash{'bypass_filtering'};
	my %used;
	my $group_num = 1;
	my $subgroup_num = 0;
	if (scalar(keys %$protein_hash) == 0){
		print "No proteins were identified for this group !!!\n";
		return (0,0);
	}
#	print "     Grouping Proteins:";
####################################
  for my $protein (sort {$$protein_hash{$b}{'occurrence'} <=> $$protein_hash{$a}{'occurrence'} || $$protein_hash{$b}{'total_nomod'} <=> $$protein_hash{$a}{'total_nomod'}
						||
		$$protein_hash{$b}{'unique_nomod'} <=> $$protein_hash{$a}{'unique_nomod'}
						||
	  $b cmp $a} keys %$protein_hash){
##########################################
#for my $protein (sort {$$protein_hash{$b}{'abundance'} <=> $$protein_hash{$a}{'abundance'}} keys %$protein_hash){

		next if (defined($used{$protein}));
		if (scalar (keys %{$$protein_hash{$protein}{'peptides'}}) == 0){
			delete ($$protein_hash{$protein});next;
		}
		print "\rGrouping Proteins $group_num: $protein                     ";
		if ($$protein_hash{$protein}{'shared'} == 0){
			$$protein_hash{$protein}{'group'} = $group_num;
			$$protein_hash{$protein}{'subgroup'} = 1;
			$used{$protein} = $group_num;
			$subgroup_num++;
		} else {
			my %shared; $shared{$protein} = 1;
			#print "$protein\n";
			get_shared_proteins(\%shared, $protein_hash, $peptide_hash, $protein);
	#		print Dumper($$protein_hash{$protein}) if ($debug);
			if ($bypass == 1){
				for my $pro (keys %shared){
					$$protein_hash{$pro}{'group'} = $group_num;
					$used{$pro} = $group_num;
				}
				$subgroup_num += subgroup_proteins(\%shared, $protein_hash);
			} else {
				open (OUT, ">clus_test");
				my $number = 1;
###### changed by yanji
				for my $pro (sort {$$protein_hash{$b}{'total_nomod'} <=> $$protein_hash{$a}{'total_nomod'} || $$protein_hash{$b}{'unique_nomod'} <=> $$protein_hash{$a}{'unique_nomod'} || $b cmp $a} keys %shared){
    				#for my $pro (sort {$$protein_hash{$b}{'occurrence'} <=> $$protein_hash{$a}{'occurrence'} || $$protein_hash{$b}{'total_nomod'} <=> $$protein_hash{$a}{'total_nomod'} || $$protein_hash{$b}{'unique_nomod'} <=> $$protein_hash{$a}{'unique_nomod'} || $b cmp $a} keys %shared){
###### end of change
			print OUT ">$pro\n";
	    		print OUT "$$protein_hash{$pro}{'sequence'}\n";
					$shared{$pro} = $number;
					$number++;
					$$protein_hash{$pro}{'group'} = $group_num;
					$used{$pro} = $group_num;
				}
	    	close OUT;
			
    		exit if (system ("/usr/local/bin/clustalw2 clus_test > clus_output"));
    		open (IN, "<clus_output");
				my %hash;
    		while(<IN>){
					last if (/Guide/);
      		next if (!/Sequences/);
      		my @array = split(/ /, $_);
					$array[1] =~ s/[\(\)]+//g;
					my @runs = split (/:/, $array[1]);
					$hash{$runs[0]}{$runs[1]} = $array[scalar(@array)-1];
    		}
				$subgroup_num += subgroup_proteins(\%shared, $protein_hash, \%hash, $percent);
			}
		}
		$group_num++;
	}
	system ("rm clus_test.aln clus_test.dnd clus_test clus_output") if (-e "clus_test.aln");
	$group_num -= 1;
	print "\n";

  my ($prev_group, $prev_subgroup) = (0,0);
	if ($bypassfil == 0){
		for my $protein (keys %$protein_hash){
######### changed by Xusheng on 11/3/2011 ############################
#			if ($protein =~ /^SW\:TRYP\_PIG/ || $protein =~ /^Random__/ || $protein =~ /^Decoy__/){
                       if ($protein =~ /^SW\:TRYP\_PIG/ || $protein =~ /Random__/ || $protein =~ /Decoy__/){
				delete $$protein_hash{$protein};
			}
		}
	}
	#regroup excluding deleted proteins
  
	my ($total_g, $total_sg) = (0,0);
	($prev_group, $prev_subgroup, $group_num, $subgroup_num) = (0,0,0,0);
	printf "Regrouping Proteins %d:                     ", $group_num+1;
######################### changed by yanji on 09192012#################################  
for my $protein (sort {$$protein_hash{$a}{'group'} <=> $$protein_hash{$b}{'group'}
#for my $protein (sort {$$protein_hash{$b}{'occurrence'} <=> $$protein_hash{$a}{'occurrence'} || $$protein_hash{$a}{'group'} <=> $$protein_hash{$b}{'group'}
#####################
            ||
      $$protein_hash{$a}{'subgroup'} <=> $$protein_hash{$b}{'subgroup'}
            ||
      $$protein_hash{$b}{'total_nomod'} <=> $$protein_hash{$a}{'total_nomod'}
            ||
      $$protein_hash{$b}{'unique_nomod'} <=> $$protein_hash{$a}{'unique_nomod'}
            ||
      $b cmp $a} keys %$protein_hash){
#for my $protein (sort {$$protein_hash{$b}{'abundance'} <=> $$protein_hash{$a}{'abundance'}} keys %$protein_hash){
			printf "\rRegrouping Proteins %d: $protein                     ", $group_num+1;
			my ($group, $subgroup) = ($$protein_hash{$protein}{'group'}, $$protein_hash{$protein}{'subgroup'});
			#print "$protein $$protein_hash{$protein}{'total_nomod'}  $$protein_hash{$protein}{'unique_nomod'}\n";
			#exit if ($group > 5);
			if ($prev_group != $group && $prev_subgroup != $subgroup){
				$total_g++; 
				$total_sg++; 
				$group_num++; 
				$subgroup_num = 1;
				$$protein_hash{$protein}{'group'} = $group_num; 
				$$protein_hash{$protein}{'subgroup'} = $subgroup_num;
				$prev_group = $group; 
				$prev_subgroup = $subgroup;
			} elsif ($prev_group != $group){
				$total_g++; $total_sg++; $group_num++; $subgroup_num = 1;
				$$protein_hash{$protein}{'group'} = $group_num; $$protein_hash{$protein}{'subgroup'} = $subgroup_num;
			} elsif ($prev_subgroup != $subgroup){
				$total_sg++; $subgroup_num++;
				$$protein_hash{$protein}{'group'} = $group_num; $$protein_hash{$protein}{'subgroup'} = $subgroup_num;
			} else {
				$$protein_hash{$protein}{'group'} = $group_num; $$protein_hash{$protein}{'subgroup'} = $subgroup_num;
			}
			$prev_group = $group; $prev_subgroup = $subgroup;
	}
	print "\n";

	return ($total_g, $total_sg);
}

sub subgroup_proteins{
	my ($group_hash, $protein_hash, $clustal_hash, $percent) = @_;
	my $subgroup = 1;
	my %notused = %$group_hash;

	for my $protein (sort {defined($$protein_hash{$b}{'total_nomod'}) <=> defined($$protein_hash{$a}{'total_nomod'})
#      	for my $protein (sort {defined($$protein_hash{$b}{'occurrence'}) <=> defined($$protein_hash{$a}{'occurrence'}) || defined($$protein_hash{$b}{'total_nomod'}) <=> defined($$protein_hash{$a}{'total_nomod'})
                                                          ||
                         defined($$protein_hash{$b}{'unique_nomod'}) <=> defined($$protein_hash{$a}{'unique_nomod'})
                                                          ||
                                                      $b cmp $a} keys %$group_hash){
		next if (!defined($notused{$protein}));
		my $first = $$group_hash{$protein};
		delete $notused{$protein}; my %hash = %notused;
		$$protein_hash{$protein}{'subgroup'} = $subgroup;
		for my $protein2 (sort {defined($$protein_hash{$b}{'total_nomod'}) <=> defined($$protein_hash{$a}{'total_nomod'})
                                                          ||
                         defined($$protein_hash{$b}{'unique_nomod'}) <=> defined($$protein_hash{$a}{'unique_nomod'})
                                                          ||
                                                      $b cmp $a} keys %$group_hash){
			next if (!defined($notused{$protein2}));
###############################
#			if($protein_hash{$protein2}{'unique'})
#			{
#				$$protein_hash{$protein2}{'subgroup'} = $subgroup;
#			}
################ added by xusheng 
			my $second = $$group_hash{$protein2};
			if (defined($percent)){
				if ($$clustal_hash{$first}{$second} > $percent){
					delete $notused{$protein2};
					$$protein_hash{$protein2}{'subgroup'} = $subgroup;
				}
			} else {
				delete $notused{$protein2};
        $$protein_hash{$protein2}{'subgroup'} = $subgroup;
			}
		}
		$subgroup++;
	}
	$subgroup -= 1;

	return ($subgroup);
}

                                                                                                                                                             
sub get_shared_proteins{
  my ($sharedpros, $protein_hash, $peptide_hash, $protein) = @_;
  
	for my $peptide (keys %{$$protein_hash{$protein}{'peptides'}}){
#		print $peptide,"\n";
		if (!defined($$peptide_hash{$peptide}{'unique'})){

########### changed by xusheng on 8/30/2012 ############
#			print "$protein\n";
#			print Dumper($$protein_hash{$protein});
#			print Dumper($$peptide_hash{$peptide});
##############
		}
		else { 
			next if ($$peptide_hash{$peptide}{'unique'} == 1);
		}
  		for my $shared (keys %{$$peptide_hash{$peptide}{'proteins'}}){
			if (!defined($$protein_hash{$protein})){
				print Dumper($$peptide_hash{$peptide});exit;
				delete $$peptide_hash{$peptide}{'proteins'}{$protein}; next;
			}
      			next if (defined($$sharedpros{$shared}));
      			$$sharedpros{$shared} = 1;
     	 		get_shared_proteins($sharedpros, $protein_hash, $peptide_hash, $shared);
  		}
	}
}

sub gen_IDHtml{
  shift @_;
  my ($Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $groupnum, $subgroup_num, $groups, $mod) = @_;

  $groups = 0 if (!defined($groups));
	if (defined($mod)){
  	if ($groups){
    	open (HTML, ">$save_dir/IDwGmod.html");
	
    	print HTML "<HTML>\n";
    	print HTML "<HEAD>\n";
    	print HTML "<TITLE>IDMOD Summary w/ Groups</TITLE>\n";
  	} else {
    	open (HTML, ">$save_dir/IDmod.html");
    	print HTML "<HTML>\n";
    	print HTML "<HEAD>\n";
    	print HTML "<TITLE>IDMOD summary w/o Groups</TITLE>\n";
  	}
	} else {
  	if ($groups){
   	 	open (HTML, ">$save_dir/IDwG.html");
   	 	print HTML "<HTML>\n";
   	 	print HTML "<HEAD>\n";
   	 	print HTML "<TITLE>ID Summary w/ Groups</TITLE>\n";
  	} else {
  	  	open (HTML, ">$save_dir/ID.html");
   	 	print HTML "<HTML>\n";
   	 	print HTML "<HEAD>\n";
    		print HTML "<TITLE>ID summary w/o Groups</TITLE>\n";
  	}
	}
	print_Legend(*HTML);
  print HTML "</HEAD>\n\n";

  my $total_protein = scalar(keys %$Protein_Hash);
  my $total_upeptide = 0;
########## changed by xusheng on 8/21/2012
#  foreach my $peptide (sort {$$Peptide_Hash{$b}{'unique'} <=> $$Peptide_Hash{$a}{'unique'}} keys %$Peptide_Hash){
foreach my $peptide (keys %$Peptide_Hash){
	next if(!defined($$Peptide_Hash{$peptide}{'unique'}));
###########
    last if ($$Peptide_Hash{$peptide}{'unique'} == 0);
    $total_upeptide++ if ($$Peptide_Hash{$peptide}{'unique'} == 1);
  }
  my $total_peptide = scalar(keys %$Peptide_Hash);

############### For debug ################

#	foreach my $peptide (keys %$Peptide_Hash)
#	{
#		print $peptide,"aa\n";
#	}

####################
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

  #$report = $report."Protein Identification Summary";

  print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	print HTML "Subgroups = $subgroup_num &nbsp&nbsp&nbsp ";
	print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	print HTML "Total Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	print HTML "Total SC = $total_SC </font></TD>\n"; #changed to total peptide DMD 5/19/05
	print HTML " </TR>\n";

	#$report = $report."Groups = $groupnum   Subgroups = $subgroup_num   Proteins = $total_protein   Total Peptides = $total_peptide   Total SC = $total_SC\n";

	print HTML "<TR>\n";
	my ($overallrate, $rate1, $rate2, $rate3) = (0,0,0,0);
#	$overallrate = ($$fprhash{'r'}/($$fprhash{'t'}-$$fprhash{'r'}))*100 if ($$fprhash{'t'} != 0);
#	$rate1 = ($$fprhash{'r1'}/($$fprhash{'t1'}-$$fprhash{'r1'}))*100 if ($$fprhash{'t1'} != 0);
#	$rate2 = ($$fprhash{'r2'}/($$fprhash{'t2'}-$$fprhash{'r2'}))*100 if ($$fprhash{'t2'} != 0);
#	$rate3 = ($$fprhash{'r3'}/($$fprhash{'t3'}-$$fprhash{'r3'}))*100 if ($$fprhash{'t3'} != 0);
#	printf HTML "<TD Align=center><Font Size=2>False Positive Rate: Overall = %.2f%% &nbsp&nbsp&nbsp ", $overallrate;
#	printf HTML "3pep = %.2f%% &nbsp&nbsp&nbsp ", $rate3;
#	printf HTML "2pep = %.2f%% &nbsp&nbsp&nbsp ", $rate2;
#	printf HTML "1pep = %.2f%% &nbsp&nbsp&nbsp ", $rate1;

      printf HTML "<TD Align=center><Font Size=2>Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $$fprhash{'final_protein_fpr'}, $$fprhash{'final_peptide_fpr'};
	#$report = $report."False Positive Rate: Overall = %.2f%%   3pep = %.2f%%   2pep = %.2f%%   1pep = %.2f%%\n";

	print HTML "</font></TD>\n"; #changed to total peptide DMD 5/19/05
	print HTML " </TR>\n";
  print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Total Peptides without C and M modifications = $total_nomod</TD>\n";
	print HTML " </TR>\n";
	print HTML "<TR><TD align=center>\n";

	#$report = $report."Total Peptides w/o C and M modifications = $total_nomod\n";
  if ($groups){
		if (defined($mod)){
###### changed by yanji, open links under local direcotry
    	#print HTML "<CENTER><A HREF=http://$server/$save_dir/IDmod.html>Hide Group Members</A></Center>";
			print HTML "<CENTER><A HREF=IDmod.html>Hide Group Members</A></Center>";
			} else {
    	#print HTML "<CENTER><A HREF=http://$server/$save_dir/ID.html>Hide Group Members</A></Center>";
	print HTML "<CENTER><A HREF=ID.html>Hide Group Members</A></Center>";
		}
  } else {
		if (defined($mod)){
    	#print HTML "<CENTER><A HREF=http://$server/$save_dir/IDwGmod.html>Show Group Members</A></Center>";
	print HTML "<CENTER><A HREF=IDwGmod.html>Show Group Members</A></Center>";
		} else {
    	#print HTML "<CENTER><A HREF=http://$server/$save_dir/IDwG.html>Show Group Members</A></Center>";
	print HTML "<CENTER><A HREF=IDwG.html>Show Group Members</A></Center>";
###### end of change
		}
  }
	print HTML "</TR></TD>\n";
  print HTML "<TR><TD Align=center>";
  print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  print HTML "</TD><TR>\n";
  print HTML "</TABLE>\n";

  print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = 95%>\n";
  print HTML "	<TR bgcolor=#EEEEEE>\n";
  print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";
  print HTML "		<TH width = 5%><B><Font Color = blue>Reference</Font></B></TH>\n";
	if (defined($$fprhash{'params'})){
  	print HTML "		<TH width = 2%><B><Font Color = blue>AX</Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}

#  print HTML "          <TH width = 2%><B><Font Color = blue>Abundance</font></Font></B></TH>\n";

  print HTML "		<TH width = 2%><B><Font Color = blue>SC</font></Font></B></TH>\n";
  print HTML "		<TH width = 2%><B><Font Color = blue>Total<br><font size=1>Peptides</font></Font></B></TH>\n";
  print HTML "		<TH width = 2%><B><Font Color = blue>Unique<br><font size=1>Peptides</font></Font></B></TH>\n";
  print HTML "		<TH width = 2%><B><Font Color = blue>Shared<br><font size=1>Peptides</font></Font></B></TH>\n";
  print HTML "		<TH><B><Font Color = blue>Description</Font></B></TH>\n";
  print HTML "		<TH width = 2%><B><Font Color = blue>Mass<br><font size=1>(kD)</font></Font></B></TH>\n";
print HTML "          <TH width = 2%><B><Font Color = blue>Abundance<br><font size=1>index</font></Font></B></TH>\n";
	if (defined($$fprhash{'params'})){
  	print HTML "		<TH width = 2%><B><Font Color = blue>N<font size=1>obsbl</font></Font></B></TH>\n" if ($$fprhash{'params'}{'abundance_index'} == 1);
	}
  print HTML "	</TR>\n";
	
  my %PrevProHash;
  my %array_abundance;
 
	my ($prev_group, $prev_subgroup) = (0,0);
	for my $protein (sort {$$Protein_Hash{$a}{'group'} <=> $$Protein_Hash{$b}{'group'}}  keys %$Protein_Hash){
		my ($group, $subgroup) = ($$Protein_Hash{$protein}{'group'}, $$Protein_Hash{$protein}{'subgroup'});
		if ($prev_group != $group){
			my $flag = 1;
			foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
	#		foreach my $protein (keys %array_abundance)
			{
				if($$Protein_Hash{$protein}{'annotation'} !~ /Uncharacterized/)
				{
					
					if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
					{
						#print "$$Protein_Hash{$protein}{'unique_nomod'}\n";
						if($flag != 1)
						{
							$subgroup++;
							$$Protein_Hash{$protein}{'subgroup'} = $subgroup;
						}
						gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, undef, $fprhash);
						delete $array_abundance{$protein}; 
					}
					$flag=0;					
				}				
			}
			if($flag==1)
			{
				foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
#				foreach my $protein (keys %array_abundance)
                        	{
					if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                                        {
						#print "$$Protein_Hash{$protein}{'unique_nomod'}\n";
                                                if($flag != 1)
                                                {
                                                	$subgroup++; 
						       $$Protein_Hash{$protein}{'subgroup'} = $subgroup + 1;
                                                }
                                       		gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, undef, $fprhash);
						delete $array_abundance{$protein};
					}
                                        $flag=0;
                                }
			}
			if($groups==1)
			{
                               foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
			#	foreach my $protein ( keys %array_abundance)
                               {
					gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, 1, $fprhash);
			       }
			}
			undef %array_abundance;
			$array_abundance{$protein} = $$Protein_Hash{$protein}{'occurrence'};
			$prev_group = $group; $prev_subgroup = $subgroup;
		}
		else
		{

			$array_abundance{$protein} = $$Protein_Hash{$protein}{'occurrence'};
		}
	}
########### for the last group ############
	my $flag=1;
        foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
#	foreach my $protein ( keys %array_abundance)   
     {
		my ($group, $subgroup) = ($$Protein_Hash{$protein}{'group'}, $$Protein_Hash{$protein}{'subgroup'});
               if($$Protein_Hash{$protein}{'annotation'} !~ /Uncharacterized/)
               {
                     if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                     {
                         if($flag != 1)
                         {
                                $$Protein_Hash{$protein}{'subgroup'} = $$Protein_Hash{$protein}{'subgroup'} + 1;
                         }

                     	gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, undef, $fprhash);
		     }
                     $flag=0;
                }
        }
        if($flag==1)
        {
             foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
#	     foreach my $protein (keys %array_abundance)
             {
		my ($group, $subgroup) = ($$Protein_Hash{$protein}{'group'}, $$Protein_Hash{$protein}{'subgroup'});
                     if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                     {
                               if($flag != 1)
                               {
                                     $$Protein_Hash{$protein}{'subgroup'} = $$Protein_Hash{$protein}{'subgroup'} + 1;
                                }

                  		 gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, undef, $fprhash);
		     }
                   $flag=0;
              }
        }
        undef %array_abundance;
	if($groups == 0){gen_IDHtml(0, $Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $groupnum, $subgroup_num, 1, $mod)};

######################

   
}
=head
############## For ID.html page #################################
			if($groups == 0)
			{
				print $groups,"\t",$prev_group,"\t",$$Protein_Hash{$protein}{'group'},"\n";

				gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, undef, $fprhash);
				if($$Protein_Hash{$protein}{'annotation'} =~ /Uncharacterized/)
				{
					if($prev_group<($group-1))
					{
		#				my $temp_pro = shift (@array_protein);
						while(my $temp_pro = shift (@array_protein))
						{
							print @array_protein,"aaaaaaa\n";
							if($$Protein_Hash{$temp_pro}{'group'}==$prev_group)
							{
			#					gen_FirstRow($Protein_Hash, $temp_pro, $save_dir, $database, 1, $fprhash);
								$flag=1;
							}
							else
							{
								gen_FirstRow($Protein_Hash, $temp_pro, $save_dir, $database, undef, $fprhash);
								$prev_group = $$Protein_Hash{$protein}{'group'};
							}
						}
						undef @array_protein;
					}
					if($prev_group == ($group-1))
					{
						$flag=1;
				#	$prev_group = $group; $prev_subgroup = $subgroup;
					}			
					
					push (@array_protein,$protein);	
				
				}
				else
				{
					$prev_group = $group; $prev_subgroup = $subgroup;
					gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, undef, $fprhash);
					if($groups == 1)
                			{
						for(my $i=0;$i<scalar(@array_protein);$i++)
						{
							print $$Protein_Hash{$array_protein[$i]}{'unique_nomod'},"\n";
                        				if($$Protein_Hash{$array_protein[$i]}{'unique_nomod'})
                        				{
                               					 $$Protein_Hash{$array_protein[$i]}{'subgroup'} = $$Protein_Hash{$array_protein[$i]}{'subgroup'} + 1;
                        				}
                       					 gen_FirstRow($Protein_Hash, $array_protein[$i], $save_dir, $database, 1, $fprhash);
						}
                       			 	$flag=0;
                			}			
				}
			} 
			elsif ($groups == 1){  #############IDwG.html ########################
                       	 	if($$Protein_Hash{$protein}{'unique_nomod'})
                        	{
                               	 	$$Protein_Hash{$protein}{'subgroup'}++;
                        	}
 #   				gen_FirstRow($Protein_Hash, $protein, $save_dir, $database, 1, $fprhash);
			}
		}
	}

	
  	print HTML "</TABLE>\n";
  	print HTML "</BODY>\n";
 	 print HTML "</HTML>\n";
 	 close HTML;
#  if($groups == 0){gen_IDHtml(0, $Protein_Hash, $Peptide_Hash, $fprhash, $save_dir, $database, $groupnum, $subgroup_num, 1, $mod)};
}
=cut
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
  print HTML "      write(\"<Font Color = blue><B><U>Show/Hide Group Members Links</U></B></Font><BR>\");\n";
  print HTML "      write(\"Follow the links to show/hide protein groups.<BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"<Font Color = blue><B><U>SC (Spectral Count)</U></B></Font><BR>\");\n";
  print HTML "      write(\"Number of times protein is sequenced (number of outfiles).<BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"<Font Color = blue><B><U>Total Peptide</U></B></Font><BR>\");\n";
  print HTML "      write(\"Total number of peptide sequences (w/o C and M modifications) that match \");\n";
  print HTML "      write(\"with the given protein.  Follow the link to view \");\n";
	print HTML "      write(\"the highlighted protein sequence and the corresponding peptide table. <BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"<Font Color = blue><B><U>Unique Peptide</U></B></Font><BR>\");\n";
  print HTML "      write(\"Total number of peptide sequences (w/o C and M modifications) that match \");\n";
  print HTML "      write(\"<Font Color = red><B>ONLY</B></Font> with the given protein.  \");\n";
	print HTML "      write(\"Follow the link to view a detailed table showing the peptides and their scores. <BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"<Font Color = blue><B><U>Shared Peptide</U></B></Font><BR>\");\n";
  print HTML "      write(\"Total number of peptide sequences (w/o C and M modifications) that match \");\n";
  print HTML "      write(\"with more than one protein. \");\n";
	print HTML "      write(\"Follow the link to view a detailed table showing the peptides and their scores. <BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"</Font>\");\n";
  print HTML "      write(\"</body></html>\");\n";
  print HTML "    }\n";
  print HTML "    descript.document.close();\n";
  print HTML "  }\n";
	print HTML "-->\n";
	print HTML "</Script>\n";
}

sub gen_FirstRow{
	my ($Protein_Hash, $protein, $save_dir, $database, $sharedrow, $fprhash) = @_;
	my $groupnum = $$Protein_Hash{$protein}{'group'}."\.".$$Protein_Hash{$protein}{'subgroup'};
	my $description = $$Protein_Hash{$protein}{'annotation'};
	my $mass = $$Protein_Hash{$protein}{'MW'};
	my $SC = $$Protein_Hash{$protein}{'occurrence'};
	my $abundance = 0;
	$abundance = $$Protein_Hash{$protein}{'abundance'};
	if (defined($$fprhash{'params'})){
  	if ($$fprhash{'params'}{'abundance_index'} == 1){
			my $top = 0;
			if ($$fprhash{'params'}{'PAI_top'} eq "TP"){
				$top = $$Protein_Hash{$protein}{'total_nomod'};
			} else {
				$top = $SC;
			}
			#if ($protein =~ /Ub_WT_YEAST/){
			#	print Dumper($$Protein_Hash{$protein});exit;
			#	print "top=$top bottom=$$Protein_Hash{$protein}{'emPAI_ftpeps'}\n";
			#}
			my $PAI = $top/$$Protein_Hash{$protein}{'emPAI_ftpeps'};
#			$abundance = (10**$PAI) - 1;
		}
	}
	
	$mass = sprintf("%3.0f", $mass/1000);
	
	if (!defined($sharedrow)){
	  print HTML "<TR BGColor = #ccffcc>\n";
	} else {
	  print HTML "<TR>\n";
	}
	#added if loop DMD 5/19/05
	if ($groupnum < 10){
		print HTML "<TD><Font Size = 1.5 Color = blue><center><B>000$groupnum</B></center></Font></TD>\n";
	} elsif ($groupnum < 100){
		print HTML "<TD><Font Size = 1.5 Color = blue><center><B>00$groupnum</B></center></Font></TD>\n";
	} elsif($groupnum < 1000) {
		print HTML "<TD><Font Size = 1.5 Color = blue><center><B>0$groupnum</B></center></Font></TD>\n";
	} else {
		print HTML "<TD><Font Size = 1.5 Color = blue><center><B>$groupnum</B></center></Font></TD>\n";
	}
	$protein =~ s/(\|ref\|)([A-Za-z0-9\_\.]+)(\|)/$1$2$3/;
	if (defined($2)){
		print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$2</B></center></Font></TD>\n";
	} else {
		print HTML "<TD><Font Size = 2.5 Color = blue><center><B>$protein</B></center></Font></TD>\n";
	}
#	if ($abundance != 0){
#		printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
#	}
	printf HTML "<TD><Font Color = blue><Center>%d</Center></Font></TD>\n", $SC;

###### added by yanji, changed path
        my $current_user = qx[whoami];
        chomp($current_user);

	if ($save_dir =~ /^\/home/) {
		$save_dir =~ s/^\/home/\/var\/www\/html/;
	} else {
		$save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
	}
#        my @dir_list = split/\//, $save_dir;
#        my $last_dir = pop @dir_list;
	
	# consider multi fractions
#	my $second_last_dir = pop @dir_list;
#	if ($second_last_dir eq "fractions") {
#		my $origin_dir = pop @dir_list;
#		$save_dir = "/var/www/html/".$current_user."/".$origin_dir."/".$second_last_dir."/".$last_dir;	
#	} else {
#        	$save_dir = "/var/www/html/".$current_user."/".$last_dir;
#	}
###### end addition
	printf HTML "$totalcgi\n", $protein, $save_dir, $database, $$Protein_Hash{$protein}{'total_nomod'};
	if (!$$Protein_Hash{$protein}{'unique'}){
		print HTML "<TD><Font Size = 4 Color = Red><Center>$$Protein_Hash{$protein}{'unique_nomod'}</Center></Font></TD>\n";
	} else {
		printf HTML "$uniquecgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'unique_nomod'};
	}
	if (!$$Protein_Hash{$protein}{'shared'}){
		print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'shared_nomod'}</Center></Font></TD>\n";
	} else {
		printf HTML "$sharedcgi\n", $protein, $save_dir, $$Protein_Hash{$protein}{'shared_nomod'};
	}
	if (!defined($description)){
		print HTML "<TD><Font Size = 2>&nbsp</Font></TD>\n";
	} else {
		$description =~ s/>/\&gt\;/g;
		$description =~ s/</\&lt\;/g;
		print HTML "<TD><Font Size = 2>$description</Font></TD>\n";
	}
	print HTML "<TD><Font Size = 3 ><Center>$mass</Center></Font></TD>\n";
        printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
    
#	if ($abundance != 0){
		#print HTML "<TD><Font Size = 3 ><Center>$$Protein_Hash{$protein}{'emPAI_ftpeps'}</Center></Font></TD>\n";
#	}
	return;
}

sub gen_SumIDHtml{
  	shift @_;
  	my ($Sum_Hash, $Protein_Hash, $Peptide_Hash, $fprhash, $folders, $save_dir, $database, $groupnum, $subgroup_num, $groups, $mod) = @_;
  	$groups = 0 if (!defined($groups));
	if (defined($mod)){
  		if ($groups){
    			open (HTML, ">$save_dir/IDwGmod.html");
    			print HTML "<HTML>\n";
    			print HTML "<HEAD>\n";
    			print HTML "<TITLE>Multiple ID Summary w/ Groups</TITLE>\n";
    			open (TXT, ">$save_dir/gmodtext.txt");
  		} else {
  			open (HTML, ">$save_dir/IDmod.html");
   			print HTML "<HTML>\n";
    			print HTML "<HEAD>\n";
    			print HTML "<TITLE>Multiple ID summary w/o Groups</TITLE>\n";
    			open (TXT, ">$save_dir/modtext.txt");
  		}
	} 
	else {
  		if ($groups){
    			open (HTML, ">$save_dir/IDwG.html");
		    	print HTML "<HTML>\n";
		    	print HTML "<HEAD>\n";
		    	print HTML "<TITLE>Multiple ID Summary w/ Groups</TITLE>\n";
		    	open (TXT, ">$save_dir/gtext.txt");
	  	} else {
  			open (HTML, ">$save_dir/ID.html");
   			print HTML "<HTML>\n";
 		   	print HTML "<HEAD>\n";
		    	print HTML "<TITLE>Multiple ID summary w/o Groups</TITLE>\n";
		    	open (TXT, ">$save_dir/text.txt");
  		}
	}
	print_SumLegend(*HTML);
	print HTML "</HEAD>\n\n";

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
	my $width = 90+($samples*10);
  	print HTML "<BODY>\n";
	print HTML "<TABLE CELLPADDING=2 ALIGN=CENTER WIDTH=90%>\n";
  	print HTML "<TR><TD Align=center><Font Size=4><B>Protein Identification Summary</font></B></Font></TD></TR>\n";
  	print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Groups = $groupnum &nbsp&nbsp&nbsp ";
	print HTML "Subgroups = $subgroup_num &nbsp&nbsp&nbsp ";
	print HTML "Proteins = $total_protein &nbsp&nbsp&nbsp ";
	print HTML "Total Peptides = $total_peptide &nbsp&nbsp&nbsp ";
	print HTML "Total SC = $total_SC</font></TD>\n"; #changed to total peptide DMD 5/19/05
	print HTML " </TR>\n";
        printf HTML "<TD Align=center><Font Size=2>Protein FDR = %.2f%%,&nbsp  Peptide FDR = %.2f%% ", $$fprhash{'protein_fpr'}, $$fprhash{'peptide_fpr'};
        print HTML "</TD>\n"; 
  	print HTML "<TR>\n";
	print HTML "<TD Align=center><Font Size=2>Total Peptides w/o C and M modifications = $total_nomod</TD>\n";
	print HTML " </TR>\n";
	print HTML "<TR><TD align=center>\n";
	if (defined($mod)){
		if ($groups){
###### changed by yanji
			#print HTML "<CENTER><A HREF=http://$server/$save_dir/IDmod.html>Hide Group Members</A></Center>";
			print HTML "<CENTER><A HREF=IDmod.html>Hide Group Members</A></Center>";
		} else {
			#print HTML "<CENTER><A HREF=http://$server/$save_dir/IDwGmod.html>Show Group Members</A></Center>";
			print HTML "<CENTER><A HREF=IDwGmod.html>Show Group Members</A></Center>";
		}
	} else {
		if ($groups){
			#print HTML "<CENTER><A HREF=http://$server/$save_dir/ID.html>Hide Group Members</A></Center>";
			print HTML "<CENTER><A HREF=ID.html>Hide Group Members</A></Center>";
		} else {
			#print HTML "<CENTER><A HREF=http://$server/$save_dir/IDwG.html>Show Group Members</A></Center>";
			print HTML "<CENTER><A HREF=IDwG.html>Show Group Members</A></Center>";
###### end of change
		}
	}
	print HTML "</TR></TD>\n";
  	print HTML "<TR><TD Align=center>";
  	print HTML "<Button style = 'width:70;align:middle;bgcolor=#EEEEEE;' onclick=mainopen()>";
  	print HTML "<Font Color = blue Size = 2><B>Legend</B></Font></Button>";
  	print HTML "</TD><TR>\n";
  	print HTML "</TABLE>\n";

  	print HTML "<TABLE BORDER=1 CELLPADDING=2 ALIGN=CENTER WIDTH = $width%>\n";
  	print HTML "	<TR bgcolor=#EEEEEE>\n";
  	print HTML "		<TH width = 2%><B><Font Color = blue size=1>Group</Font></B></TH>\n";

  	print HTML "		<TH width = 5%%><B><Font Color = blue>Reference</Font></B></TH>\n";
### changed on 8/9/2012 
# print HTML "          <TH width = 5%%><B><Font Color = blue>Abundance</Font></B></TH>\n";
########################
#  print HTML "		<TH width = 2%><B><Font Color = blue>Total<br><font size=1>Peptides</font></Font></B></TH>\n";
	print TXT "group;reference;totalpep;totalcoverage;";
	print HTML "          <TH width = 2%><B><Font Color = blue>Total<br><font size=1>SC</font></Font></B></TH>\n";
	for (my $i=0; $i<scalar(@$folders); $i++){
#		print TXT "$$folders[$i]_TP;";
		print TXT "$$folders[$i]_SC;";
		print TXT "$$folders[$i]_Coverage;";
#  		print HTML "		<TH width = 2%><B><Font Color = blue>TP<br><font size=1>$$folders[$i]</font></Font></B></TH>\n";
  		print HTML "		<TH width = 2%><B><Font Color = blue>SC<br><font size=1>$$folders[$i]</font></Font></B></TH>\n";
	}

        print HTML "          <TH width = 2%><B><Font Color = red>Total<br><font size=1>Peptides</font></Font></B></TH>\n";
        for (my $i=0; $i<scalar(@$folders); $i++){
		print TXT "$$folders[$i]_TP;";
		print HTML "            <TH width = 2%><B><Font Color = red>TP<br><font size=1>$$folders[$i]</font></Font></B></TH>\n";
		print TXT "$$folders[$i]_Coverage;";
	}

  	print HTML "		<TH><B><Font Color = blue>Description</Font></B></TH>\n";
 	print HTML "		<TH width = 2%><B><Font Color = blue>Mass<br><font size=1>(kD)</font></Font></B></TH>\n";
###### added by yanji
	print HTML "          <TH width = 2%><B><Font Color = blue>Abundance<br><font size=1>index</font></Font></B></TH>\n";
######
  	print HTML "	</TR>\n";

	print TXT "Description;KD;Length\n";
  	my %PrevProHash;
	my %array_abundance;
	my ($prev_group, $prev_subgroup) = (0,0);
  	for my $protein (sort {$$Protein_Hash{$a}{'group'} <=> $$Protein_Hash{$b}{'group'}
						||
			$$Protein_Hash{$a}{'subgroup'} <=> $$Protein_Hash{$b}{'subgroup'}
						||
			$$Protein_Hash{$b}{'total_nomod'} <=> $$Protein_Hash{$a}{'total_nomod'}
						||
	   	  	$$Protein_Hash{$b}{'unique_nomod'} <=> $$Protein_Hash{$a}{'unique_nomod'}
						||
		 	$b cmp $a} keys %$Protein_Hash)
		{
		my ($group, $subgroup) = ($$Protein_Hash{$protein}{'group'}, $$Protein_Hash{$protein}{'subgroup'});
                if ($prev_group != $group){
                        my $flag = 1;
                        foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
                        {
                                if($$Protein_Hash{$protein}{'annotation'} !~ /Uncharacterized/)
                                {
                                        if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                                        {
                                                if($flag != 1)
                                                {
                                                        $subgroup++;
                                                        $$Protein_Hash{$protein}{'subgroup'} = $subgroup;
                                                }
						gen_SumFirstRow($Sum_Hash, $Protein_Hash, $folders, $protein, $save_dir, $database);                 
                                                delete $array_abundance{$protein};
                                        }
                                        $flag=0;
                                }
                        }
                        if($flag==1)
                        {
                                foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
#                               foreach my $protein (keys %array_abundance)
                                {
                                        if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                                        {
                                                if($flag != 1)
                                                {
                                                        $subgroup++;
                                                       $$Protein_Hash{$protein}{'subgroup'} = $subgroup + 1;
                                                }
						gen_SumFirstRow($Sum_Hash, $Protein_Hash, $folders, $protein, $save_dir, $database);
                                                delete $array_abundance{$protein};
                                        }
                                        $flag=0;
                                }
                        }
                        if($groups==1)
                        {
                               foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
                               {
					gen_SumFirstRow($Sum_Hash, $Protein_Hash, $folders, $protein, $save_dir, $database, 1);
                               }
                        }
                        undef %array_abundance;
                        $array_abundance{$protein} = $$Protein_Hash{$protein}{'occurrence'};
                        $prev_group = $group; $prev_subgroup = $subgroup;
                }
                else
                {

                        $array_abundance{$protein} = $$Protein_Hash{$protein}{'occurrence'};
                }
        }
########### for the last group ############
        my $flag=1;
        foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
	{
                my ($group, $subgroup) = ($$Protein_Hash{$protein}{'group'}, $$Protein_Hash{$protein}{'subgroup'});
               if($$Protein_Hash{$protein}{'annotation'} !~ /Uncharacterized/)
               {
                     if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                     {
                         if($flag != 1)
                         {
                                $$Protein_Hash{$protein}{'subgroup'} = $$Protein_Hash{$protein}{'subgroup'} + 1;
                         }
			gen_SumFirstRow($Sum_Hash, $Protein_Hash, $folders, $protein, $save_dir, $database);
                     }
                     $flag=0;
                }
        }
        if($flag==1)
        {
             foreach my $protein (reverse sort {$array_abundance{$a}<=>$array_abundance{$b}} keys %array_abundance)
             {
                my ($group, $subgroup) = ($$Protein_Hash{$protein}{'group'}, $$Protein_Hash{$protein}{'subgroup'});
                     if($flag || $$Protein_Hash{$protein}{'unique_nomod'}>0)
                     {
                               if($flag != 1)
                               {
                                     $$Protein_Hash{$protein}{'subgroup'} = $$Protein_Hash{$protein}{'subgroup'} + 1;
                                }
			gen_SumFirstRow($Sum_Hash, $Protein_Hash, $folders, $protein, $save_dir, $database);
                     }
                   $flag=0;
              }
        }
        undef %array_abundance;


  	print HTML "</TABLE>\n";
  	print HTML "</BODY>\n";
  	print HTML "</HTML>\n";
  	close HTML;
  	if($groups == 0){gen_SumIDHtml(0, $Sum_Hash, $Protein_Hash, $Peptide_Hash,$fprhash,$folders, $save_dir, $database, $groupnum, $subgroup_num, 1, $mod)};
}

sub print_SumLegend{
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
  print HTML "      write(\"<Font Color = blue><B><U>Show/Hide Group Members Links</U></B></Font><BR>\");\n";
  print HTML "      write(\"Follow the links to show/hide protein groups.<BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"<Font Color = blue><B><U>Total Peptide</U></B></Font><BR>\");\n";
  print HTML "      write(\"Total number of peptide sequences (w/o C and M modifications) that match \");\n";
  print HTML "      write(\"with the given protein.  Follow the link to view \");\n";
	print HTML "      write(\"the highlighted protein sequence and the corresponding peptide table. <BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"<Font Color = blue><B><U>SC (Spectral Count)</U></B></Font><BR>\");\n";
  print HTML "      write(\"Number of times protein is sequenced (number of outfiles).<BR>\");\n";
  print HTML "      write(\"<BR>\");\n";
  print HTML "      write(\"</body></html>\");\n";
  print HTML "    }\n";
  print HTML "    descript.document.close();\n";
  print HTML "  }\n";
	print HTML "-->\n";
	print HTML "</Script>\n";
}

sub gen_SumFirstRow{
	my ($Sum_Hash, $Protein_Hash, $folders, $protein, $save_dir, $database, $sharedrow) = @_;
	my $groupnum = $$Protein_Hash{$protein}{'group'}."\.".$$Protein_Hash{$protein}{'subgroup'};
	my $description = $$Protein_Hash{$protein}{'annotation'};
	my $mass = $$Protein_Hash{$protein}{'MW'};
	$mass = sprintf("%3.0f", $mass/1000);
	
########### added by yanji
	my $abundance = $$Protein_Hash{$protein}{'abundance'};
######	

	if (!defined($sharedrow)){
	  print HTML "<TR BGColor = #ccffcc>\n";
	} 
	else {
	  print HTML "<TR>\n";
	}
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

###### added by yanji, changed path
        my $current_user = qx[whoami];
        chomp($current_user);

        if ($save_dir =~ /^\/home/) {
                $save_dir =~ s/^\/home/\/var\/www\/html/;
        } else {
                $save_dir = "\/var\/www\/html\/".$current_user.$save_dir;
        }
##############	

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
			printf HTML "$uniquecgi\n", $protein, $save_dir, $$Sum_Hash{$msdak}{'protein'}{$protein}{'total_nomod'};
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
	printf HTML "<TD><Font Color = blue><Center>%.2f</Center></Font></TD>\n", $abundance;
######	
	return;
}
1;
