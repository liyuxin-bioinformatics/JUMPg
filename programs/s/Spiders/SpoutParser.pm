#!/usr/local/bin/perl

## Release date: 01/31/2015
## Release version: version 11.1.1
## Module name: Spiders::SpoutParser

package Spiders::SpoutParser;
        
use strict;
use warnings;
use Storable;
use vars qw($VERSION @ISA @EXPORT);
use Excel::Writer::XLSX;

$VERSION     = 2.01;
@ISA	 = qw(Exporter);
@EXPORT      = ();

my ($Hydrogen_mass, $HydrogenMinus_mass, $Cterm_mass)=(1.007825032,1.00335,19.017806);
my ($now_date, $now_time) = &gettime();

BEGIN
{
}

sub new{
    my ($class,%arg) = @_;
    my $self = {
    };	
	$self->{'totalScan'}=$self->{'MS1scanNumber'}=$self->{'MS2scanNumber'}=$self->{'MS1MS2ratio'}=$self->{'chargeDistribution'}=$self->{'ppiDistribution'}=0;
    bless $self, $class;
    return $self;
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

sub AA_mass{
my ($AA)=@_;
my (%AA_mass);
		$AA_mass{'G'} = 57.02146372057;	
		$AA_mass{'D'} = 115.02694302383;
		$AA_mass{'A'} = 71.03711378471;
		$AA_mass{'Q'} = 128.05857750528;
		$AA_mass{'S'} = 87.03202840427;
		$AA_mass{'K'} = 128.094963014;
		$AA_mass{'P'} = 97.05276384885;
		$AA_mass{'E'} = 129.04259308797;
		$AA_mass{'V'} = 99.06841391299;
		$AA_mass{'M'} = 131.04048491299;
		$AA_mass{'T'} = 101.04767846841;
		$AA_mass{'H'} = 137.05891185845;
		$AA_mass{'C'} = 103.00918478471;
		$AA_mass{'F'} = 147.06841391299;
		$AA_mass{'I'} = 113.08406397713;		
		$AA_mass{'L'} = 113.08406397713;
		$AA_mass{'J'} = 113.08406397713;
		$AA_mass{'R'} = 156.1011110236;
		$AA_mass{'N'} = 114.04292744114;	
		$AA_mass{'Y'} = 163.06332853255;
		$AA_mass{'W'} = 186.07931294986;
return $AA_mass{$AA} ;
}


sub enzymeStat
{
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


sub readSpsearchParam
{
        my ($self, $seqPara,$paramF)=@_;
        open(SQT,"$paramF") or die 'Cannot open $paramF\n';

        my @array1=('peptide_mass_tolerance','fragment_ion_tolerance','ion_series','max_num_differential_AA_per_mod','num_output_lines','remove_precursor_peak','ion_cutoff_percentage',
        'match_peak_count','match_peak_allowed_error','match_peak_tolerance','protein_mass_filter');
        $seqPara->{'optionalParams'}=\@array1;

        my @array2=('xcorr','deltacn','deltacnstar','spscore','sprank','lSideMass',   'rSideMass',  'PeptideEvalue', 'TagSeq',  'WeightedEvalue',  'TagsNum');
        $seqPara->{'hitOptionalParams'}=\@array2;

        my ($DB_local_path);
        while (<SQT>){
                #last if (/SEQUEST_ENZYME_INFO/);
                next if (/^#/ || /^\s+\Z/);
                chomp;  s/[;#].*//; s/\s+\Z//;          #Q: remove empty lines
                my ($key, $value) = split(" = ", $_);
                next if (!defined($value));

                #push(@$seqPara{'optionalParams'},$key);
                $$seqPara{$key}=$value;
        }

        #Tryptic KR P
        $$seqPara{'enzyme_info'} =~ /^(\w+)\s(\w+)/; $$seqPara{'enzyme_info'} = $1;
        $$seqPara{'enzyme'}=lc($$seqPara{'enzyme_info'});
        $$seqPara{'misClv'}=$$seqPara{'max_mis_cleavage'};
        $$seqPara{'frmTol'}=$$seqPara{'frag_mass_tolerance'};

        close(SQT);
}

sub parseSpout
{
        my ($self, $seqParaHash, $outInforHash, $run)=@_;

        #deal with out files
        opendir(DIR,$run); my @out = grep {/\.spout\Z/} readdir(DIR);
        closedir(DIR);

        if (scalar(@out) == 0){ die "     THERE ARE NO .spout FILES IN DIRECTORY $run !!!!\n";}
        my ($orignum, $count) = (scalar(@out), 0);

        #extract general params from an out file
        genParOut($seqParaHash,$run,$out[0]);

		my %temphash; # used for checking duplications
        #read outfiles
        for my $outfile (@out)
        {
                $count++;
                print "\r  Gathering information from $count of $orignum outfiles     ";

                readOutfile($seqParaHash, $outInforHash, $run, $outfile,\%temphash);
        }

}


sub genParOut {
	my ($paraHash,$run,$outfile)=@_;
	my $fullpath = "$run\/$outfile";                #can have a problem
	open (IN, $fullpath);
	my $mod;
	while (<IN>) {
		chomp;
		if (/^JUMP version/) {
			$$paraHash{'search_engine'}=$_;
		} elsif (/^Database\s=/) {
			#Database = /home/xwang4/JUMP_database/ratmouse/ratmouse_con.fasta.mdx
			my @tmp=split(/\s=\s/,$_);
			$$paraHash{'first_database_name'}=$tmp[1];
		} elsif (/^ion series ABCDVWXYZ/) {
			my @tmp=split(/\:/,$_);
			$$paraHash{'ion_series'}=$tmp[1];
			$mod=<IN>;
			$mod=<IN>;
			chomp($mod);
		}
	}
	
	#modification
	if (defined($mod)) {
		$mod =~ s/\s*\/\s*/\//g; 
		$mod =~ s/(\s)\s+/$1/g; 
		$mod =~ s/^\s+//; 
		$mod =~ s/\s+\Z//;
		#dynamic_M@=15.99492  add_C_Cysteine=57.02146
		my @tmp=split(/\s/,$mod);
		for my $term (@tmp) {
			if ($term =~ /dynamic_(\w+)([\@\#\*\^\~\$\%\&\?\!\(\)\{\}\[\]\:\;\'\<\>])\=([0-9\.]+)/) {	#dynamic modification
				my ($AA,$symble,$massdiff)=($1,$2,$3);
				$$paraHash{'mod'}=1;
				my @tmp2=split(//,$AA);
				for (my $i=0; $i<scalar(@tmp2); $i++) {
					if (defined($$paraHash{'modification'}{'dynamic'}{'count'})) {
						$$paraHash{'modification'}{'dynamic'}{'count'}++;
					} else {
						$$paraHash{'modification'}{'dynamic'}{'count'}=1;
					}
					$$paraHash{'modification'}{'dynamic'}{'AA'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$tmp2[$i];
					$$paraHash{'modification'}{'dynamic'}{'symble'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$symble;
					$$paraHash{'modification'}{'dynamic'}{'massdiff'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$massdiff;
					$$paraHash{'modification'}{'dynamic'}{'mass'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=AA_mass($tmp2[$i])+$massdiff;
				}
			} elsif ($term =~ /add_(\w)_(\w+)\=([0-9\.]+)/) {	#static modification
				$$paraHash{'mod'}=1;
				if (defined($$paraHash{'modification'}{'static'}{'count'})) {
					$$paraHash{'modification'}{'static'}{'count'}++;
				} else {
					$$paraHash{'modification'}{'static'}{'count'}=1;
				}
				$$paraHash{'modification'}{'static'}{'AA'}[$$paraHash{'modification'}{'static'}{'count'}]=$1;
				$$paraHash{'modification'}{'static'}{'massdiff'}[$$paraHash{'modification'}{'static'}{'count'}]=$3;
				$$paraHash{'modification'}{'static'}{'mass'}[$$paraHash{'modification'}{'static'}{'count'}]=$3+AA_mass($1);
			} elsif ($term =~ /add_Nterm_peptide\=([0-9\.]+)/) {	#Nterm mod
				$$paraHash{'mod'}=1;
				$$paraHash{'modification'}{'Nterm'}{'massdiff'}=$1;
				$$paraHash{'modification'}{'Nterm'}{'mass'}=$1+$Hydrogen_mass;
			} elsif ($term =~ /add_Cterm_peptide\=([0-9\.]+)/) {	#Cterm mod
				$$paraHash{'mod'}=1;
				$$paraHash{'modification'}{'Cterm'}{'massdiff'}=$1;
				$$paraHash{'modification'}{'Cterm'}{'mass'}=$1+$Cterm_mass;
			}
		}
	} else {
		print "No modification infor provided.\n";
	}
	close IN;
}


sub readOutfile
{
        my ($paraHash,$outInfor, $run, $outfile,$temphash)=@_;
        my $maxHitCondiered=5;                  #max number of hits to be considered

		
        my $fullpath = "$run\/$outfile";                #can have a problem
        open (IN, $fullpath);
        #Rat_B_100ng_Q.10000.1.3.spout
        $outfile =~ /\.(\d+)\.(\d+)\.(\d)\.spout/;
        my ($scan,$prOrder,$charge) = ($1,$2,$3);
        $$outInfor{$outfile}{'scan'}=$scan;
        $$outInfor{$outfile}{'prOrder'}=$prOrder;
        $$outInfor{$outfile}{'charge'}=$charge;
		
        #print "$scan,$charge\n";
        while (<IN>)
        {
                if (/^Precursor mass\s*=\s*([0-9\.]+)[\s\t]+Percentage of precursor peak intensity\s*=\s*([\d\.\%]+)/)
                {
                        $$outInfor{$outfile}{'MHmass'}=$1;
                        $$outInfor{$outfile}{'PPIpct'}=$2;
                }
                elsif (/^MS2 signal noise ratio\s*=\s*([0-9\.]+)[\,\s\t]+Tag number\s*=\s*([\d\.\%]+)/)
                {
                        $$outInfor{$outfile}{'MS2SN'}=$1;
                        $$outInfor{$outfile}{'TagNum'}=$2;
                }
		elsif (/^Precursor matched peptides = (.*?), Peptides mass tolerance =/)
		{
			$$outInfor{$outfile}{'Precursor_matched_peptides'}=$1;
		}
		elsif (/^PTM sites = (.*?),  PTM score = (.*?)$/)
		{
			$$outInfor{$outfile}{'PTM_sites'}=$1;
			$$outInfor{$outfile}{'PTM_score'}=$2;
		}
                elsif (/^-+\n/)
                {
                        last;   #here comes the hits
                }
        }
        #unless (defined($$outInfor{$outfile}{'MHmass'})) {die "MHmass not found in $fullpath\n";}
        #unless (defined($$outInfor{$outfile}{'PPIpct'})) {die "PPIpct not found in $fullpath\n";}

        my $i=0; my $count;
        $$outInfor{$outfile}{'hitShown'}=0;
        while (<IN>)
        {
                last if (/^\s*$/);                                      #?
                last if (/All identified peptide/);
                chomp; #Remove return line
                s/^\s+//; #Remove leading spaces
                my @vals = split /\s+/,$_;

                if(scalar(@vals)>11){       die "Error parsing file $outfile at line:  Some values are missing\n"; }
                elsif (scalar(@vals)==11)               #hit lines
                {
                        $i++; $$outInfor{$outfile}{'hitShown'}++; #hit counts
                        #Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  WeightedEvalue  TagsNum      Reference                       Peptide
                        #-------------------------------------------------------------------------------------------------------------------------------------------------
                        #1      2127.994    983.4356    1145.5659    0.6093    S           1.89           1        ##Decoy__sp|P51650-2|SSDH_RAT   K.ELYED@IGYSKGGYCVYVL.-
                        #Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq    TagRank  TagNum  Jscore      Reference               Peptide
                        #-------------------------------------------------------------------------------------------------------------------------------------------------
                        #1      1361.666    302.0983    1060.5754    6.06    G       17    1      6.06       tr|C9JT76|C9JT76_HUMAN  K.ADDGPCKAIMKR.F
                        my $peptide = $vals[10];
												
						
						my $origpeptide=$peptide;
                        #K.ELYEDIGYSKGGYCVYVL.-
						$peptide =~ s/([A-Z\-])\.([A-Z\@\#\*\^\~\$\%\&\?\!\(\)\{\}\[\]\:\;\'\<\>]+)\.([A-Z\-])/$2/;
                        $$outInfor{$outfile}{'prAA'}[$i]=$1;
                        $$outInfor{$outfile}{'nxAA'}[$i]=$3;

############# not keep the same IDs for the scan when second search applied by adjust the mono-isotopic value 						
						next if (defined($temphash->{$scan}->{$peptide}));
						
                        if ( defined($$outInfor{$outfile}{'fullPeptide'}[$i-1]) and $peptide eq $$outInfor{$outfile}{'fullPeptide'}[$i-1]) #same peptide
                        {
                                $i--; $$outInfor{$outfile}{'hitShown'}--; #hit counts

                                #alternativ proteins
                                $count++;
                                $$outInfor{$outfile}{'altPro'}[$i][$count]=$vals[9];
                                $$outInfor{$outfile}{'red'}[$i]=$count;

                                next;
                        }

                        $count=0;       #for alternative proteins

                        #print "\n$outfile\n$peptide,$$outInfor{$outfile}{'prAA'}[$i],$$outInfor{$outfile}{'nxAA'}[$i]\n";

                        $$outInfor{$outfile}{'protein'}[$i]=$vals[9];
			if ( $vals[9] =~ /Decoy/ ) { $$outInfor{$outfile}{'random'}[$i]=1; } else { $$outInfor{$outfile}{'random'}[$i]=0; }
                        $$outInfor{$outfile}{'origpeptide'}[$i]=$origpeptide;


                        $$outInfor{$outfile}{'hitRank'}[$i]=$i;
                        $$outInfor{$outfile}{'rank'}[$i]=$i;
                        #$$outInfor{$outfile}{'sprank'}[$i]='';
                        $$outInfor{$outfile}{'expmass'}[$i]=$vals[1];
                        #$$outInfor{$outfile}{'xcorr'}[$i]='';

                        $$outInfor{$outfile}{'lSideMass'}[$i]=$vals[2];
                        $$outInfor{$outfile}{'rSideMass'}[$i]=$vals[3];
                        $$outInfor{$outfile}{'PeptideEvalue'}[$i]=$vals[4];
                        $$outInfor{$outfile}{'TagSeq'}[$i]=$vals[5];
                        $$outInfor{$outfile}{'TagRank'}[$i]=$vals[6];						
                        $$outInfor{$outfile}{'WeightedEvalue'}[$i]=$vals[8];
                        $$outInfor{$outfile}{'TagsNum'}[$i]=$vals[7];
                        $$outInfor{$outfile}{'fullPeptide'}[$i]=$peptide;

                        #miss cleavage
                        my $testpeptide = $$outInfor{$outfile}{'fullPeptide'}[$i];
                        #$testpeptide =~ s/.\.([A-Z\@\#\*\^\~\$]+)\../$1/;
                        $$outInfor{$outfile}{'misClv'}[$i] = ($testpeptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;


                        $$outInfor{$outfile}{'deltacn'}[$i]=0;
                        #$$outInfor{$outfile}{'deltacnstar'}[$i]=0;
                        if ($i>1)
                        {
                                if ($$outInfor{$outfile}{'WeightedEvalue'}[$i-1] !=0 )
                                {
                                        $$outInfor{$outfile}{'deltacn'}[$i-1]=($$outInfor{$outfile}{'WeightedEvalue'}[$i-1]-$$outInfor{$outfile}{'WeightedEvalue'}[$i])/$$outInfor{$outfile}{'WeightedEvalue'}[$i-1];
                                }
                        }


                        $$outInfor{$outfile}{'matchIons'}[$i]=0;
                        $$outInfor{$outfile}{'totalIons'}[$i]=0;
                        $$outInfor{$outfile}{'red'}[$i]=0;

                        #consider modifications:
                        if (defined($$paraHash{'modification'}{'dynamic'}{'count'}))            #dynamic mod:
                        {
	                        $$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]=0;
        	                my @tmp=split(//,$peptide);

                	        my $j=0;
                        	while ($j<scalar(@tmp))
	                        #for (my $j=0; $j<scalar(@tmp); $j++)
        	                {
                	                for (my $k=1; $k<=$$paraHash{'modification'}{'dynamic'}{'count'}; $k++)
                        	        {
                                	        if ($tmp[$j] eq $$paraHash{'modification'}{'dynamic'}{'AA'}[$k])
                                        	{
                                                	if ($j+1<scalar(@tmp) and $tmp[$j+1] eq $$paraHash{'modification'}{'dynamic'}{'symble'}[$k])
	                                                {
        	                                                $j++;
                	                                        $$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]++;
                        	                                $$outInfor{$outfile}{'dynamicMod'}{'modIndex'}[$i][$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]]=$k;
                                	                        $$outInfor{$outfile}{'dynamicMod'}{'pos'}[$i][$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]]=$j-$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]+1;
                                        	                $$outInfor{$outfile}{'dynamicMod'}{'mass'}[$i][$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]]=$$paraHash{'modification'}{'dynamic'}{'mass'}[$k];

                                                	        #print "\n$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]\:$i,$j,$k";
                                                        	last;
	                                                }
        	                                }
                	                }
                        	        $j++;
                        	}

	                        $peptide =~ s/\W//g;
        	        }

                        if (defined($$paraHash{'modification'}{'static'}{'count'}))             #static mod:
                        {       #       print "\nstatic mod";
	                        $$outInfor{$outfile}{'staticMod'}{'found'}[$i]=0;
        	                my @tmp=split(//,$peptide);
                	        for (my $j=0; $j<scalar(@tmp); $j++)
                        	{
                                	for (my $k=1; $k<=$$paraHash{'modification'}{'static'}{'count'}; $k++)
	                                {
        	                                my $a=$$paraHash{'modification'}{'static'}{'AA'}[$k];
                	                        if ($tmp[$j] eq $a)
                        	                {
                                	                $$outInfor{$outfile}{'staticMod'}{'found'}[$i]++;
                                        	        $$outInfor{$outfile}{'staticMod'}{'modIndex'}[$i][$$outInfor{$outfile}{'staticMod'}{'found'}[$i]]=$k;
                                                	$$outInfor{$outfile}{'staticMod'}{'pos'}[$i][$$outInfor{$outfile}{'staticMod'}{'found'}[$i]]=$j+1;
	                                                $$outInfor{$outfile}{'staticMod'}{'mass'}[$i][$$outInfor{$outfile}{'staticMod'}{'found'}[$i]]=$$paraHash{'modification'}{'static'}{'mass'}[$k];

        	                                        #print "\n$$outInfor{$outfile}{'staticMod'}{'found'}[$i]\:$i,$j,$k";
                	                                last;
                        	                }
                                	}
	                        }
        	        }

                        $$outInfor{$outfile}{'peptide'}[$i]=$peptide;

                        my @AA_KR=$peptide =~ /([KR])/g;
                        my @AA_KRP=$peptide =~ /([KR]P)/g;
                        $$outInfor{$outfile}{'num_tol_term'}[$i]=scalar(@AA_KR)-scalar(@AA_KRP);

                }
                else            #alternative proteins
                {
                        $count++;
                        $$outInfor{$outfile}{'altPro'}[$i][$count]=$vals[0];
                        $$outInfor{$outfile}{'red'}[$i]=$count;
                }
        }

        if ($$outInfor{$outfile}{'hitShown'}==1) {$$outInfor{$outfile}{'deltacn'}[1]=1;}


        close(IN);
}

sub printTable
{
	my ($self, $outInfor, $run, $output)=@_;

        opendir(DIR,$run); my @outs = grep {/\.spout\Z/} readdir(DIR);
        closedir(DIR);

	open(OUT,">$output");

	print OUT "outfile\tPrecursor_matched_peptides\tmatchedType\tscore\tdeltaScore\tpeptide\tmatchedTagSeq\tmatchedTagLength\n";
	#foreach my $out (%$outInfor)
	foreach my $out (@outs)
	{
		print OUT "$out\t$$outInfor{$out}{Precursor_matched_peptides}\t";
		if ( $$outInfor{$out}{hitShown}>0 )
		{
			my $matchedType;
			if ( $$outInfor{$out}{protein}[1] =~ m/Decoy/ ) { $matchedType='decoy'; }
			else { $matchedType='target'; }

			my ($matchedTagSeq,$matchedTagLength);
			if ( $$outInfor{$out}{TagSeq}[1] eq 'N/A' )
			{
				$matchedTagSeq='NA';
				$matchedTagLength=0;
			}
			else
			{
				$matchedTagSeq=$$outInfor{$out}{TagSeq}[1];
				$matchedTagLength=length($matchedTagSeq);
			}

			print OUT "$matchedType\t$$outInfor{$out}{WeightedEvalue}[1]\t$$outInfor{$out}{deltacn}[1]\t$$outInfor{$out}{peptide}[1]\t$matchedTagSeq\t$matchedTagLength\n";
		}
		#else { print OUT "unmatched","\tNA"x5,"\n"; }
		else { print OUT "unmatched","\tNA"x4,"\t0\n"; }
	}

	close OUT;
}

sub printPepXML
{
        my ($self, $seqPara, $outInfor, $run, $outXML, $topHit)=@_;
		#my ($now_date, $now_time) = ($self->{'now_date'},$self->{'now_time'});

        open(XML,">$outXML");

        print XML '<?xml version="1.0" encoding="Accoding to Out2XML(TPP v2.9 GALE rev.3, Build 200611091255)"?>',"\n";
        print XML '<msms_pipeline_analysis date="',"$now_date $now_time",'" ','summary_xml="',$outXML,'">',"\n";
        print XML '<msms_run_summary base_name="',$run,'" raw_data_type="raw" raw_data="">',"\n";

        print XML '<sample_enzyme name="',$$seqPara{'enzyme'},'">',"\n";
        print XML '<specificity cut="',enzymeStat($$seqPara{'enzyme'},'cut');
        unless (enzymeStat($$seqPara{'enzyme'},'nocuts') eq '')
        {
                print XML '" no_cut="',enzymeStat($$seqPara{'enzyme'},'nocuts');
        }
        print XML '" sense="',enzymeStat($$seqPara{'enzyme'},'sense'),'"/>',"\n";
        print XML '</sample_enzyme>',"\n";

        print XML '<search_summary base_name="',$run,'" search_engine="',$$seqPara{search_engine},'" precursor_mass_type="monoisotopic" fragment_mass_type="monoisotopic" out_data_type="spout" out_data=".spout" search_id="">',"\n";
        print XML '<search_database local_path="',$$seqPara{'first_database_name'},'" type="AA"/>',"\n";
        if (defined($$seqPara{'modification'}{'dynamic'}{'count'}))
        {
                for (my $i=1; $i<=$$seqPara{'modification'}{'dynamic'}{'count'}; $i++)
                {
                        print XML '<aminoacid_modification aminoacid="',$$seqPara{'modification'}{'dynamic'}{'AA'}[$i],'" massdiff="',$$seqPara{'modification'}{'dynamic'}{'massdiff'}[$i],'" mass="',$$seqPara{'modification'}{'dynamic'}{'mass'}[$i],'" variable="Y" symbol="',$$seqPara{'modification'}{'dynamic'}{'symble'}[$i],'"/>',"\n";
                }
        }
        if (defined($$seqPara{'modification'}{'static'}{'count'}))
        {
                for (my $i=1; $i<=$$seqPara{'modification'}{'static'}{'count'}; $i++)
                {
                        print XML '<aminoacid_modification aminoacid="',$$seqPara{'modification'}{'static'}{'AA'}[$i],'" massdiff="',$$seqPara{'modification'}{'static'}{'massdiff'}[$i],'" mass="',$$seqPara{'modification'}{'static'}{'mass'}[$i],'" variable="N"/>',"\n";
                }
        }
        if (defined($$seqPara{'modification'}{'Nterm'}{'mass'}))
        {
                print XML '<terminal_modification terminus="n" massdiff="',$$seqPara{'modification'}{'Nterm'}{'massdiff'},'" mass="',$$seqPara{'modification'}{'Nterm'}{'mass'},'" variable="N" protein_terminus="N"/>',"\n";
        }
        if (defined($$seqPara{'modification'}{'Cterm'}{'mass'}))
        {
                print XML '<terminal_modification terminus="c" massdiff="',$$seqPara{'modification'}{'Cterm'}{'massdiff'},'" mass="',$$seqPara{'modification'}{'Cterm'}{'mass'},'" variable="N" protein_terminus="C"/>',"\n";
        }

        #OPTIONAL PARAMETERS
        foreach my $tmpTitle (@{$seqPara->{'optionalParams'}})
        {#print "$tmpTitle,$$seqPara{$tmpTitle}\n";
                if (defined($$seqPara{$tmpTitle}))
                {
                        print XML "\<parameter name=\"$tmpTitle\" value=\"$$seqPara{$tmpTitle}\"\/\>\n";
                }
        }

        print XML '</search_summary>',"\n";

        opendir(DIR,$run); my @out = grep {/\.spout\Z/} readdir(DIR);
        closedir(DIR);

        my %tmpoutfiles;
        for my $outfile (@out)
        {
                #Rat_B_100ng_Q.10000.10000.3.spout
                $outfile =~ /\.(\d+)\.(\d+)\.(\d+)\.\w+/;
                $tmpoutfiles{$outfile}=$1;
        }
        @out=sort {$tmpoutfiles{$a} <=> $tmpoutfiles{$b}} keys %tmpoutfiles;


        #@out=sort(@out);
        my $count=0;
        for my $outfile (@out)
        {
                $count++;
                $outfile =~ s/(.*?)(\.spout)/$1$2/;
				
               	next if(!defined($$outInfor{$outfile}{'PPIpct'}));
                next if(!defined($$outInfor{$outfile}{'prOrder'}));
                #print XML "\<spectrum_query spectrum=\"$1\" start_scan=\"$$outInfor{$outfile}{'scan'}\" end_scan=\"$$outInfor{$outfile}{'scan'}\" precursor_neutral_mass=\"",$$outInfor{$outfile}{'MHmass'}-$Hydrogen_mass,"\" precursor_peak_intensity_percentage=\"$$outInfor{$outfile}{'PPIpct'}\" precursor_peak_intensity_order=\"$$outInfor{$outfile}{'prOrder'}\" assumed_charge=\"$$outInfor{$outfile}{'charge'}\" index=\"$count\"\>\n";
                print XML "\<spectrum_query spectrum=\"$1\" start_scan=\"$$outInfor{$outfile}{'scan'}\" end_scan=\"$$outInfor{$outfile}{'scan'}\" precursor_neutral_mass=\"",$$outInfor{$outfile}{'MHmass'}-$Hydrogen_mass,"\" precursor_peak_intensity_percentage=\"$$outInfor{$outfile}{'PPIpct'}\" precursor_peak_intensity_order=\"$$outInfor{$outfile}{'prOrder'}\" precursor_matched_peptides=\"$$outInfor{$outfile}{'Precursor_matched_peptides'}\" assumed_charge=\"$$outInfor{$outfile}{'charge'}\" index=\"$count\"\>\n";
                print XML '<search_result>',"\n";

                next if(!defined($$outInfor{$outfile}{'hitShown'}));
                for (my $i=1; $i<=$$outInfor{$outfile}{'hitShown'} and $i<=$topHit; $i++ )
                {
						$$outInfor{$outfile}{'peptide'}[$i]="" if(!defined($$outInfor{$outfile}{'peptide'}[$i]));
						$$outInfor{$outfile}{'prAA'}[$i]="" if(!defined($$outInfor{$outfile}{'prAA'}[$i]));
						$$outInfor{$outfile}{'nxAA'}[$i]="" if(!defined($$outInfor{$outfile}{'nxAA'}[$i]));
						$$outInfor{$outfile}{'protein'}[$i]="" if(!defined($$outInfor{$outfile}{'protein'}[$i]));
						$$outInfor{$outfile}{'red'}[$i]=0 if(!defined($$outInfor{$outfile}{'red'}[$i]));
						$$outInfor{$outfile}{'matchIons'}[$i]=0 if(!defined($$outInfor{$outfile}{'matchIons'}[$i]));
						$$outInfor{$outfile}{'totalIons'}[$i]=0 if(!defined($$outInfor{$outfile}{'totalIons'}[$i]));
						$$outInfor{$outfile}{'expmass'}[$i]=0 if(!defined($$outInfor{$outfile}{'expmass'}[$i]));
						$$outInfor{$outfile}{'MHmass'}=0 if(!defined($$outInfor{$outfile}{'MHmass'}));
						$$outInfor{$outfile}{'num_tol_term'}[$i]=0 if(!defined($$outInfor{$outfile}{'num_tol_term'}[$i]));
						$$outInfor{$outfile}{'misClv'}[$i]=0 if(!defined($$outInfor{$outfile}{'misClv'}[$i]));						
						
                        unless (defined($$outInfor{$outfile}{'peptide'}[$i]) and defined($$outInfor{$outfile}{'prAA'}[$i]) and defined($$outInfor{$outfile}{'nxAA'}[$i])) {print "\n$outfile\n";}
                        print XML "\<search_hit hit_rank=\"",$i, '" peptide="',$$outInfor{$outfile}{'peptide'}[$i],'" peptide_prev_aa="',$$outInfor{$outfile}{'prAA'}[$i],'" peptide_next_aa="',$$outInfor{$outfile}{'nxAA'}[$i],
                        '" protein="',$$outInfor{$outfile}{'protein'}[$i],'" num_tot_proteins="',$$outInfor{$outfile}{'red'}[$i]+1,'" num_matched_ions="',$$outInfor{$outfile}{'matchIons'}[$i],'" tot_num_ions="',
                        $$outInfor{$outfile}{'totalIons'}[$i],'" calc_neutral_pep_mass="',$$outInfor{$outfile}{'expmass'}[$i] - $Hydrogen_mass,'" massdiff="',$$outInfor{$outfile}{'MHmass'}-$$outInfor{$outfile}{'expmass'}[$i],
                        '" num_tol_term="',$$outInfor{$outfile}{'num_tol_term'}[$i],'" num_missed_cleavages="',$$outInfor{$outfile}{'misClv'}[$i],'" is_rejected="',0,"\"\>\n";

                        #alternative_protein
                        for (my $j=1; $j<=$$outInfor{$outfile}{'red'}[$i]; $j++)
                        {
                                unless (defined($$outInfor{$outfile}{'altPro'}[$i][$j])) {print "$outfile,red=$$outInfor{$outfile}{'red'}[$i],i=$i,j=$j\n";next;}
                                print XML '<alternative_protein protein="',$$outInfor{$outfile}{'altPro'}[$i][$j],'"/>',"\n";
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
                                                print XML '<mod_aminoacid_mass position="',$$outInfor{$outfile}{'dynamicMod'}{'pos'}[$i][$j],'" mass="',$$outInfor{$outfile}{'dynamicMod'}{'mass'}[$i][$j],'"/>',"\n";
                                        }
                                }

                                if (defined($$outInfor{$outfile}{'staticMod'}{'found'}[$i]))
                                {#print "\nmark,$$outInfor{$outfile}{'staticMod'}{'found'}[$i]";
                                        for (my $j=1; $j<=$$outInfor{$outfile}{'staticMod'}{'found'}[$i]; $j++)
                                        {
                                                print XML '<mod_aminoacid_mass position="',$$outInfor{$outfile}{'staticMod'}{'pos'}[$i][$j],'" mass="',$$outInfor{$outfile}{'staticMod'}{'mass'}[$i][$j],'"/>',"\n";
                                        }
                                }



                                print XML '</modification_info>',"\n";
                        }

                        foreach my $tmpTitle (@{$seqPara->{'hitOptionalParams'}})
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

sub set_MS1scanNumber
{
	my ($self,$MS1scanNumber)=@_;
	$self->{'MS1scanNumber'}=$MS1scanNumber;	
}

sub set_MS2scanNumber
{
        my ($self,$MS2scanNumber)=@_;
        $self->{'MS2scanNumber'}=$MS2scanNumber;
}

=head
sub set_MS1MS2ratio
{
	my ($self,$MS1MS2ratio)=@_;
        $self->{'MS1MS2ratio'}=$MS1MS2ratio;
}
=cut

sub set_chargeDistribution
{
        my ($self,$chargeDistribution)=@_;
        $self->{'chargeDistribution'}=$chargeDistribution;
}

sub set_ppiDistribution
{
        my ($self,$ppiDistribution)=@_;
        $self->{'ppiDistribution'}=$ppiDistribution;
}

sub set_totalScan
{
	my ($self,$totalScan)=@_;
        $self->{'totalScan'}=$totalScan;
}

sub num2letter
{
        my ($num)=@_;
        return chr($num+64);
}

sub printXLSX
{
    my ($self, $seqPara, $outInfor, $run, $FileName, $topMS2_array)=@_;
#	my ($now_date, $now_time) = ($self->{'now_date'},$self->{'now_time'});



    # Create a new Excel file
    #my $FileName = "$run\/Report\.xlsx";
#	  my $FileName = "Report\.xlsx";
    my $workbook = Excel::Writer::XLSX->new($FileName);

    # set tmp path
    #$workbook->set_tempdir( '/data/XLSXtmp' );

    # Add 'Summary' worksheet
    my $worksheet0 = $workbook->add_worksheet('Summary Table');
    $worksheet0->set_column('A:A',35);
	
    # Add 'MS2 quality' worksheet
    my $worksheet3 = $workbook->add_worksheet('MS2 Qualtity');
    $worksheet3->set_column('A:A',35);
	
    # Add a worksheet
    my $worksheet1 = $workbook->add_worksheet('PSM');
    $worksheet1->set_column('E:F', 35);

    # Define the format and add it to the worksheet
    my $titleFormat = $workbook->add_format();
    $titleFormat->set_bold();
    my $dcmFormat2 = $workbook->add_format();
    $dcmFormat2->set_num_format('0.00');
    my $dcmFormat4 = $workbook->add_format();
    $dcmFormat4->set_num_format('0.0000');
    my $ctrFormat = $workbook->add_format();
    $ctrFormat->set_align('center');

    # 'Summary' worksheet
	if ($self->{'MS1scanNumber'} and $self->{'MS2scanNumber'}) 
	{
		$self->{'MS1MS2ratio'}=$self->{'MS2scanNumber'}/$self->{'MS1scanNumber'};
	}

    my $i0=0;
    $worksheet0->write($i0,0,"Summary table",$titleFormat); $i0++;
	$i0++;
    $worksheet0->write($i0,0,'Software: JUMP v2.0',$titleFormat);$i0++;
    $worksheet0->write($i0,0,"Date and time: $now_date $now_time",$titleFormat);$i0++;
	$i0++;
    $worksheet0->write($i0,0,'Input file',$titleFormat); $worksheet0->write($i0,1,$run);$i0++;
    $worksheet0->write($i0,0,'Database',$titleFormat);  $worksheet0->write($i0,1,$$seqPara{'first_database_name'});$i0++;
	$i0++;
    $worksheet0->write($i0,0,"Total scan",$titleFormat); $worksheet0->write($i0,1, $self->{'totalScan'}); $i0++;
    $worksheet0->write($i0,0,"Number of MS1 scan",$titleFormat); $worksheet0->write($i0,1, $self->{'MS1scanNumber'}); $i0++;
    $worksheet0->write($i0,0,"Number of MS2 scan",$titleFormat); $worksheet0->write($i0,1, $self->{'MS2scanNumber'}); $i0++;
    $worksheet0->write($i0,0,"MS2/MS1 ratio",$titleFormat); $worksheet0->write($i0,1, $self->{'MS1MS2ratio'},$dcmFormat2); $i0++;
	for (my $tmp=1; $tmp<=10; $tmp++)
	{
		$worksheet0->write($i0,0,"TOP$tmp",$titleFormat); $worksheet0->write($i0,1, $$topMS2_array[$tmp-1]); $i0++;
	}

	$i0++;
    $worksheet0->write($i0,0,"Charge distribution",$titleFormat); $worksheet0->write($i0,1, $self->{'chargeDistribution'}); $i0++;
    $worksheet0->write($i0,0,"PPI distribution",$titleFormat); $worksheet0->write($i0,1, $self->{'ppiDistribution'}); $i0++;

    opendir(DIR,$run); my @out = grep {/\.spout\Z/} readdir(DIR);
    closedir(DIR);

	
	
    my %tmpoutfiles;
    for my $outfile (@out)
    {
            #Rat_B_100ng_Q.10000.10000.3.spout
            $outfile =~ /\.(\d+)\.(\d+)\.(\d+)\.\w+/;
            $tmpoutfiles{$outfile}=$1;
		$$outInfor{$outfile}{'PPI'}=$2;
    }
    @out=sort {$tmpoutfiles{$a} <=> $tmpoutfiles{$b}} keys %tmpoutfiles;

	# calculate q values
	my %smallhash;
	my %pepHash;
	while ( my ($outfile, $hash) = each %$outInfor ) 
	{
		next if (!defined($$outInfor{$outfile}{'WeightedEvalue'}[1]));
		$smallhash{$outfile}{peptide}=$$outInfor{$outfile}{'peptide'}[1]; 		
		$smallhash{$outfile}{WeightedEvalue}=$$outInfor{$outfile}{'WeightedEvalue'}[1]; 
		$smallhash{$outfile}{protein}=$$outInfor{$outfile}{'protein'}[1]; 
	}
	_calculate_qValue(\%smallhash,'WeightedEvalue','target');
	buildNomodHash(\%smallhash,\%pepHash,'WeightedEvalue');
	_calculate_qValue(\%pepHash,'WeightedEvalue','target');
    #@out=sort(@out);
    my $i=0;
    $worksheet1->write($i,0,"Summary of Peptide Spectra Matches \(PSMs\)",$titleFormat); 
    $worksheet3->write($i,0,"Summary of MS2 quality",$titleFormat); 
	$i++;
    $i++;
    $worksheet1->write($i,0,"Program: Jump",$titleFormat);     $worksheet3->write($i,0,"Program: Jump",$titleFormat); $i++;
    $worksheet1->write($i,0,"Date and time: $now_date $now_time",$titleFormat); $worksheet3->write($i,0,"Date and time: $now_date $now_time",$titleFormat); $i++;
    $worksheet1->write($i,0,"Data set: $run",$titleFormat); $worksheet3->write($i,0,"Data set: $run",$titleFormat); $i++;
    $i++;

    my $j=0;
    $worksheet1->write($i,$j,"Order",$titleFormat); $worksheet3->write($i,$j,"Order",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Scan \#",$titleFormat); $worksheet3->write($i,$j,"Scan \#",$titleFormat); $j++;
    $worksheet1->write($i,$j,"PPI \#",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Charge",$titleFormat); $worksheet3->write($i,$j,"MS2 S/N",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Peptide Sequence",$titleFormat); $worksheet3->write($i,$j,"Tag Number",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Protein",$titleFormat); $worksheet3->write($i,$j,"Matched Tag rank",$titleFormat);$j++;
    $worksheet1->write($i,$j,"Jscore",$titleFormat); $worksheet3->write($i,$j,"Matched Tag Sequence",$titleFormat); $j++;
    $worksheet1->write($i,$j,"DeltaCN",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Tag Sequence",$titleFormat); $j++;
    $worksheet1->write($i,$j,"exp. MH\+",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Theo MH\+",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Mass shift",$titleFormat); $j++;
    $worksheet1->write($i,$j,"Miscleavage",$titleFormat); $j++;
    $worksheet1->write($i,$j,"PSM Q value",$titleFormat); $j++;	
    $worksheet1->write($i,$j,"Peptide Q value",$titleFormat); $j++;		
   # $worksheet1->write($i,$j,"",$titleFormat); $j++;




   
   
    my $count=0; my (%histogram,%wev); $wev{'decoy'}[0]=$wev{'target'}[0]=0; $wev{'min'}=100; $wev{'max'}=0;
    for (my $tmp=1; $tmp<=5; $tmp++) { $histogram{'charge'}[$tmp]=0; }
    for (my $tmp=0; $tmp<=$$seqPara{'misClv'}; $tmp++) { $histogram{'misClv'}[$tmp]=0; }
    my $maxTagL=0; my $tagL=0; #my %pp;
	my (@precursorMass,$pcmI);$pcmI=0;
	my %petidehash; #$petidehash{D}=$petidehash{T}=0;
    for my $outfile (sort {$smallhash{$b}{WeightedEvalue} <=> $smallhash{$a}{WeightedEvalue}} keys %smallhash)
    {
        if ($$outInfor{$outfile}{'hitShown'})
		{
			$$outInfor{$outfile}{'prAA'}[1] ="" if(!defined($$outInfor{$outfile}{'prAA'}[1]));
			$$outInfor{$outfile}{'fullPeptide'}[1] ="" if(!defined($$outInfor{$outfile}{'fullPeptide'}[1]));
			$$outInfor{$outfile}{'nxAA'}[1] ="" if(!defined($$outInfor{$outfile}{'nxAA'}[1]));
			
			$$outInfor{$outfile}{'prAA'}[1] ="" if(!defined($$outInfor{$outfile}{'prAA'}[1]));
			$$outInfor{$outfile}{'prAA'}[1] ="" if(!defined($$outInfor{$outfile}{'prAA'}[1]));
			
            $i++;$count++;$j=0;
            $worksheet1->write($i,$j,$count,$ctrFormat);  $worksheet3->write($i,$j,$count,$ctrFormat); $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'scan'},$ctrFormat); $worksheet3->write($i,$j,$$outInfor{$outfile}{'scan'},$ctrFormat); $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'PPI'},$ctrFormat); $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'charge'},$ctrFormat); $worksheet3->write($i,$j,$$outInfor{$outfile}{'MS2SN'}); $j++;
            $worksheet1->write($i,$j,"$$outInfor{$outfile}{'prAA'}[1]\.$$outInfor{$outfile}{'fullPeptide'}[1]\.$$outInfor{$outfile}{'nxAA'}[1]"); unless ($$outInfor{$outfile}{'TagNum'} eq '0') {$worksheet3->write($i,$j,$$outInfor{$outfile}{'TagNum'})}; $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'protein'}[1]); unless ($$outInfor{$outfile}{'TagRank'}[1] eq '0') {$worksheet3->write($i,$j,$$outInfor{$outfile}{'TagRank'}[1])}; $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'WeightedEvalue'}[1]); unless ($$outInfor{$outfile}{'TagSeq'}[1] eq 'N/A') {$worksheet3->write($i,$j,$$outInfor{$outfile}{'TagSeq'}[1])} ;$j++;		
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'deltacn'}[1],$dcmFormat2); $j++;
            unless ($$outInfor{$outfile}{'TagSeq'}[1] eq 'N/A')
            {$worksheet1->write($i,$j,$$outInfor{$outfile}{'TagSeq'}[1]);} $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'MHmass'}); $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'expmass'}[1]); $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'MHmass'}-$$outInfor{$outfile}{'expmass'}[1],$dcmFormat4); $j++;
            $worksheet1->write($i,$j,$$outInfor{$outfile}{'misClv'}[1],$ctrFormat); $j++;
            $worksheet1->write($i,$j,$smallhash{$outfile}{'qvalue'},$dcmFormat4); $j++;

		# peptide FDR
		if ($smallhash{$outfile}{protein} =~ /Decoy/) { $petidehash{D}{$$outInfor{$outfile}{'fullPeptide'}[1]}=''; }
		else { $petidehash{T}{$$outInfor{$outfile}{'fullPeptide'}[1]}='';}
		my $D_pep=scalar(keys %{$petidehash{D}}); my $T_pep=scalar(keys %{$petidehash{T}});
            $worksheet1->write($i,$j,$D_pep/($T_pep),$dcmFormat4); $j++;
			$worksheet1->write($i,$j,$pepHash{"$$outInfor{$outfile}{'prAA'}[1]\.$$outInfor{$outfile}{'fullPeptide'}[1]\.$$outInfor{$outfile}{'nxAA'}[1]"}{'qvalue'},$dcmFormat4); $j++;
            #$worksheet1->write($i,$j,$$outInfor{$outfile}{'peptide'}[1]); $j++;



	#precursorMass
			$precursorMass[$pcmI]=($$outInfor{$outfile}{'expmass'}[1]-$HydrogenMinus_mass)/$$outInfor{$outfile}{'charge'} + $HydrogenMinus_mass;
			$pcmI++;

            #uniq peptide and proteins
            #unless (defined($pp{'peptide'}{$$outInfor{$outfile}{'fullPeptide'}[1]})) {$pp{'peptide'}{$$outInfor{$outfile}{'fullPeptide'}[1]}=1;}
            #unless (defined($pp{'protein'}{$$outInfor{$outfile}{'protein'}[1]})) {$pp{'protein'}{$$outInfor{$outfile}{'protein'}[1]}=1;}

            #histogram
            $histogram{'charge'}[$$outInfor{$outfile}{'charge'}]++;
            $histogram{'misClv'}[$$outInfor{$outfile}{'misClv'}[1]]++;




            if ($$outInfor{$outfile}{'TagSeq'}[1] eq 'N/A') {$tagL=0;}
            else {$tagL=length($$outInfor{$outfile}{'TagSeq'}[1]);}
            if (defined( $histogram{'tagL'}[$tagL])) {$histogram{'tagL'}[$tagL]++;}
            else {$histogram{'tagL'}[$tagL]=1;}
            if ($tagL>$maxTagL) {$maxTagL=$tagL;}

            #weighted E value
            if ($$outInfor{$outfile}{'WeightedEvalue'}[1]>$wev{'max'}) {$wev{'max'}=$$outInfor{$outfile}{'WeightedEvalue'}[1];}
            elsif ($$outInfor{$outfile}{'WeightedEvalue'}[1]<$wev{'min'}) {$wev{'min'}=$$outInfor{$outfile}{'WeightedEvalue'}[1];}

            if ($$outInfor{$outfile}{'protein'}[1] =~ /Decoy/) {$wev{'decoy'}[0]++; $wev{'decoy'}[$wev{'decoy'}[0]]=$$outInfor{$outfile}{'WeightedEvalue'}[1];}
            else {$wev{'target'}[0]++; $wev{'target'}[$wev{'target'}[0]]=$$outInfor{$outfile}{'WeightedEvalue'}[1];}
        }

    }

	#precursorMass
	@precursorMass=sort {$a <=> $b} @precursorMass;
	my $pcmRangeCut=0.0005;
	my (%pcmHash,@pcmArray,$pcmSlot); my $pcmHb=0; 
	$pcmSlot=10;
	if(scalar(@precursorMass)==0)
	{
		$pcmArray[0]=0; 
	}
	else
	{
		for (my $k=0; $k<=int($precursorMass[scalar(@precursorMass)-1]/$pcmSlot); $k++) 
		{
			$pcmArray[$k]=0; 
		}
	}	
	for (my $k=0; $k<scalar(@precursorMass); $k++)
	{
		$pcmArray[int($precursorMass[$k]/$pcmSlot)]++;
		if ( $precursorMass[$k]<$pcmHb )
		{
			
		}
		else
		{
			$pcmHb=$precursorMass[$k]+$pcmRangeCut;
			$pcmHash{$precursorMass[$k]}=1;
		}
	}

        # 'Summary' worksheet
      $worksheet0->write($i0,0,'Number of outfiles',$titleFormat); $worksheet0->write($i0,1,scalar(@out));$i0++;
	$worksheet0->write($i0,0,'Number of unmatched',$titleFormat); $worksheet0->write($i0,1,scalar(@out)-$count);
	my $tmpPct=fstr((scalar(@out)-$count)*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)");  $i0++;
	$worksheet0->write($i0,0,'Number of unique precursor ions',$titleFormat); $worksheet0->write($i0,1,scalar(keys %pcmHash)); 
	$tmpPct=fstr(scalar(keys %pcmHash)*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)");  $i0++;
        $worksheet0->write($i0,0,'Number of PSMs',$titleFormat); $worksheet0->write($i0,1,$count);
	$tmpPct=fstr($count*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)"); $i0++;
=head
        my $td_line=$i0; $i0+=6;
        $worksheet0->write($i0,0,'Number of unique peptide',$titleFormat); $worksheet0->write($i0,1,scalar(keys %{$pp{'peptide'}}));$i0++;
        $worksheet0->write($i0,0,'Number of unique protein',$titleFormat); $worksheet0->write($i0,1,scalar(keys %{$pp{'protein'}}));$i0++;
=cut

        #my $chargeChart = $workbook->add_chart( type => 'column', name => 'Charge distribution');
        my $chargeChart = $workbook->add_chart( type => 'column', name => 'Charge distribution', embedded => 1); $chargeChart->set_title( name => 'Charge distribution' );
        my $misClvChart = $workbook->add_chart( type => 'column', name => 'Misscleavage distribution', embedded => 1);$misClvChart->set_title( name => 'Misscleavage distribution' );

        # Add a worksheet
        my $worksheet2 = $workbook->add_worksheet('charts');
        #$worksheet2->insert_chart( 'E2', $chargeChart );
        #$worksheet2->insert_chart( 'E20', $misClvChart );

        my @markRows;
        $i=$j=0; $markRows[1]=$i+1;
        $worksheet2->write($i,0,'Charge',$titleFormat);$i++;
        for (my $tmp=1; $tmp<=5; $tmp++,$i++) { $worksheet2->write($i,1,$tmp);  $worksheet2->write($i,2,$histogram{'charge'}[$tmp]);}
        $i=0; $markRows[2]=$i+1; my $colOffset=15;
        $worksheet2->write($i,$colOffset,'Miss Cleavage',$titleFormat);$i++;
        for (my $tmp=0; $tmp<=$$seqPara{'misClv'}; $tmp++,$i++) { $worksheet2->write($i,$colOffset+1,$tmp);  $worksheet2->write($i,$colOffset+2,$histogram{'misClv'}[$tmp]);}


        # Configure the chart.
        $chargeChart->add_series(
                categories => '=charts!$B$2:$B$6',
                values     => '=charts!$C$2:$C$6',
                name       => 'charge',
        );
        my $a1=num2letter($colOffset+2); my $a2=num2letter($colOffset+3);
        my $b1=$markRows[2]+1; my $b2=$markRows[2]+$$seqPara{'misClv'}+1;
         $misClvChart->add_series(
                 categories => "\=charts\!\$$a1\$$b1\:\$$a1\$$b2",
                 values     => "\=charts\!\$$a2\$$b1\:\$$a2\$$b2",
                 name       => 'miss cleavage',
         );


        $a1=num2letter($colOffset+5);
        $worksheet2->insert_chart( "E$markRows[1]", $chargeChart );
        $worksheet2->insert_chart( "$a1$markRows[2]", $misClvChart );

        #tag length dstr
        $i=$i+20;
        $worksheet2->write($i,0,'Tag sequence length',$titleFormat);$i++;$markRows[5]=$i;
        for (my $tmp=0; $tmp<=$maxTagL; $tmp++,$i++)
        {
                $worksheet2->write($i,1,$tmp);
                if (defined($histogram{'tagL'}[$tmp])) { $worksheet2->write($i,2,$histogram{'tagL'}[$tmp]);}
                else {$worksheet2->write($i,2,'0');}
        }

        #graphing
        my $tagLengthChart = $workbook->add_chart( type => 'column', name => 'Tag sequence length', embedded => 1); $tagLengthChart->set_title( name => 'Tag sequence length distributiom' );
        $a1='B'; $a2='C';
        $b1=$markRows[5]+1; $b2=$markRows[5]+$maxTagL+1;
        $tagLengthChart->add_series(
                 categories => "\=charts\!\$$a1\$$b1\:\$$a1\$$b2",
                 values     => "\=charts\!\$$a2\$$b1\:\$$a2\$$b2",
                 name       => 'tag length',
        );
        $worksheet2->insert_chart( "E$markRows[5]", $tagLengthChart );


        #weighted evalue dstr
        my $binN=100; my $binSlot=($wev{'max'}-$wev{'min'}+0.001)/$binN;
        for (my $k=0; $k<$binN; $k++) { $wev{'decoydstr'}[$k]= $wev{'targetdstr'}[$k]=0;}
        for (my $k=1; $k<=$wev{'decoy'}[0]; $k++){$wev{'decoydstr'}[int($wev{'decoy'}[$k]/$binSlot)]++;}
        for (my $k=1; $k<=$wev{'target'}[0]; $k++){$wev{'targetdstr'}[int($wev{'target'}[$k]/$binSlot)]++;}

        #to cut the tail
        my $cutPct=0.99; $wev{'targetCml'}[0]=$wev{'targetdstr'}[0];
        for (my $k=1; $k<$binN; $k++) { $wev{'targetCml'}[$k]=$wev{'targetCml'}[$k-1]+$wev{'targetdstr'}[$k]; }
        $wev{'decoyCml'}[0]=$wev{'decoydstr'}[0];
        for (my $k=1; $k<$binN; $k++) { $wev{'decoyCml'}[$k]=$wev{'decoyCml'}[$k-1]+$wev{'decoydstr'}[$k]; }

        my $markBin=0;
        for (my $k=0; $k<$binN; $k++) {$wev{'targetCml'}[$binN-1] = 1 if($wev{'targetCml'}[$binN-1]==0); if  ($wev{'targetCml'}[$k]/$wev{'targetCml'}[$binN-1]>=$cutPct) {$markBin=$k; last;} }
        #my $rest=$wev{'targetCml'}[$binN-1]-$wev{'targetCml'}[$markBin];
        my @binBd; for (my $k=0; $k<$binN; $k++) {$binBd[$k]=$k*$binSlot;}
        my $cutValue="\>$binBd[$markBin]";

        my $wevChart = $workbook->add_chart( type => 'line', name => 'Jscore distribution', embedded => 1);$wevChart->set_title( name => 'Jscore distribution' );

        #$i=$i+20;
        $i=$markRows[5]-1;
        $worksheet2->write($i,$colOffset,'Jscore',$titleFormat);$i++;
        $worksheet2->write($i,$colOffset+2,'decoy',$titleFormat);$worksheet2->write($i,$colOffset+3,'target',$titleFormat);$i++;
        $markRows[3]=$i+1;
        for (my $k=0; $k<=$markBin; $k++,$i++)
        {
                $worksheet2->write($i,$colOffset+1,$binBd[$k]);
                $worksheet2->write($i,$colOffset+2,$wev{'decoydstr'}[$k]);
                $worksheet2->write($i,$colOffset+3,$wev{'targetdstr'}[$k]);
        }
        $worksheet2->write($i,$colOffset+1,$cutValue);
        $worksheet2->write($i,$colOffset+2,$wev{'decoyCml'}[$binN-1]-$wev{'decoyCml'}[$markBin]);
        $worksheet2->write($i,$colOffset+3,$wev{'targetCml'}[$binN-1]-$wev{'targetCml'}[$markBin]);
        $i++;

        $b1=$markRows[3]; $b2=$b1+$markBin+1;
        $a1=num2letter($colOffset+2); $a2=num2letter($colOffset+3); my $a3=num2letter($colOffset+4);
        $wevChart->add_series(
                 categories =>  "\=charts\!\$$a1\$$b1\:\$$a1\$$b2",
                 values     =>  "\=charts\!\$$a2\$$b1\:\$$a2\$$b2",
                 name       => 'decoy',
         );
        $wevChart->add_series(
                  categories =>  "\=charts\!\$$a1\$$b1\:\$$a1\$$b2",
                  values     =>  "\=charts\!\$$a3\$$b1\:\$$a3\$$b2",
                  name       => 'target',
          );

        $a1=num2letter($colOffset+5);
        $worksheet2->insert_chart( "$a1$markRows[3]", $wevChart );

	#precursor mass
	$i=$i+10;#$markRows[6]=$i;
	$worksheet2->write($i,0,'Precursor mass',$titleFormat);$i++;$markRows[6]=$i;
	my $number = (scalar(@precursorMass)==0) ? 1 : scalar(@precursorMass);
	for (my $k=int($precursorMass[0]/$pcmSlot); $k<=int($precursorMass[scalar(@precursorMass)-1]/$pcmSlot); $k++,$i++)
	{
		$worksheet2->write($i,1,$k*$pcmSlot);
		$worksheet2->write($i,2,$pcmArray[$k]);
	}

	my $pcmChart = $workbook->add_chart( type => 'line', name => 'Precursor mass', embedded => 1); $pcmChart->set_title( name => 'Precursor mass distribution' );
	$a1='B'; $a2='C';
        $b1=$markRows[6]+1; $b2=$markRows[6]+int($precursorMass[scalar(@precursorMass)-1]/$pcmSlot)-int($precursorMass[0]/$pcmSlot)+1;
        $pcmChart->add_series(
                 categories => "\=charts\!\$$a1\$$b1\:\$$a1\$$b2",
                 values     => "\=charts\!\$$a2\$$b1\:\$$a2\$$b2",
                 name       => 'precursor mass',
        );
        $worksheet2->insert_chart( "E$markRows[6]", $pcmChart );
 
        #FDR calculation
        #0.1 slot wev dstr
        $binSlot=0.1; $binN=($wev{'max'}-$wev{'min'}+0.001)/$binSlot;
        for (my $k=0; $k<$binN; $k++) { $wev{'decoydstr'}[$k]= $wev{'targetdstr'}[$k]=0;}
        for (my $k=1; $k<=$wev{'decoy'}[0]; $k++){$wev{'decoydstr'}[int($wev{'decoy'}[$k]/$binSlot)]++;}
        for (my $k=1; $k<=$wev{'target'}[0]; $k++){$wev{'targetdstr'}[int($wev{'target'}[$k]/$binSlot)]++;}
        #cumulative
        for (my $k=0; $k<$binN; $k++) { $wev{'targetCml'}[$k]=$wev{'decoyCml'}[$k]=0;$wev{'fdr'}[$k]=1; }
        $wev{'targetCml'}[0]=$wev{'targetdstr'}[0];
        for (my $k=1; $k<$binN; $k++) { $wev{'targetCml'}[$k]=$wev{'targetCml'}[$k-1]+$wev{'targetdstr'}[$k]; }
        $wev{'decoyCml'}[0]=$wev{'decoydstr'}[0];
        for (my $k=1; $k<$binN; $k++) { $wev{'decoyCml'}[$k]=$wev{'decoyCml'}[$k-1]+$wev{'decoydstr'}[$k]; }
        #reverse cml
        for (my $k=0; $k<$binN; $k++)
        {
                $wev{'targetRevCml'}[$k]=$wev{'target'}[0]-$wev{'targetCml'}[$k]+$wev{'targetdstr'}[$k];
                $wev{'decoyRevCml'}[$k]=$wev{'decoy'}[0]-$wev{'decoyCml'}[$k]+$wev{'decoydstr'}[$k];
        }
        #fdr
        for (my $k=$binN-1; $k>=0; $k--)
        {
                $wev{'fdr'}[$k]=$wev{'decoyRevCml'}[$k]/($wev{'targetRevCml'}[$k]);
        }
=head
        for (my $k=0; $k<$binN; $k++)
        {
                print $wev{'min'}+$k*$binSlot,"\t",$wev{'targetRevCml'}[$k],"\t",$wev{'decoyRevCml'}[$k],"\t",$wev{'fdr'}[$k],"\n";
        }
=cut

        $worksheet0->write($i0,0,'Number of target PSMs',$titleFormat); $worksheet0->write($i0,1,$wev{'target'}[0]);
	$tmpPct=fstr($wev{'target'}[0]*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)"); $i0++;
        $worksheet0->write($i0,0,'Number of decoy PSMs',$titleFormat); $worksheet0->write($i0,1,$wev{'decoy'}[0]);
	$tmpPct=fstr($wev{'decoy'}[0]*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)");$i0++;
	$worksheet0->write($i0,0,'Number of expected target PSMs',$titleFormat); $worksheet0->write($i0,1,$wev{'target'}[0]-$wev{'decoy'}[0]);$i0++;

        my (%pp, $fdrCut, $wevCut);

        $fdrCut=0.05; $wevCut=0; undef(%pp);
        for (my $k=0; $k<$binN; $k++) { if ($wev{'fdr'}[$k]<=$fdrCut) {$wevCut=$wev{'min'}+$k*$binSlot;$markBin=$k; last;} }
	$i0++;
        $worksheet0->write($i0,0,'5% FDR',$titleFormat);$i0++;
        $worksheet0->write($i0,0,'Number of validated target',$titleFormat); $worksheet0->write($i0,1,$wev{'targetRevCml'}[$markBin]);
	$tmpPct=fstr($wev{'targetRevCml'}[$markBin]*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)"); $i0++;	
        $worksheet0->write($i0,0,'Number of validated decoy',$titleFormat); $worksheet0->write($i0,1,$wev{'decoyRevCml'}[$markBin]);
	$tmpPct=fstr($wev{'decoyRevCml'}[$markBin]*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)"); $i0++;
        for my $outfile (@out)
        {
                if ($$outInfor{$outfile}{'hitShown'} and $$outInfor{$outfile}{'WeightedEvalue'}[1]>=$wevCut)
                {
                        unless (defined($pp{'peptide'}{$$outInfor{$outfile}{'fullPeptide'}[1]})) {$pp{'peptide'}{$$outInfor{$outfile}{'fullPeptide'}[1]}=1;}
                        unless (defined($pp{'protein'}{$$outInfor{$outfile}{'protein'}[1]})) {$pp{'protein'}{$$outInfor{$outfile}{'protein'}[1]}=1;}
                }
        }
        $worksheet0->write($i0,0,'Number of unique peptide',$titleFormat); $worksheet0->write($i0,1,scalar(keys %{$pp{'peptide'}}));$i0++;
        $worksheet0->write($i0,0,'Number of unique protein',$titleFormat); $worksheet0->write($i0,1,scalar(keys %{$pp{'protein'}}));$i0++;

        $fdrCut=0.01; $wevCut=0;undef(%pp);
        for (my $k=0; $k<$binN; $k++) { if ($wev{'fdr'}[$k]<=$fdrCut) {$wevCut=$wev{'min'}+$k*$binSlot;$markBin=$k; last;} }
	$i0++;
        $worksheet0->write($i0,0,'1% FDR',$titleFormat);$i0++;
        $worksheet0->write($i0,0,'Number of validated target',$titleFormat); $worksheet0->write($i0,1,$wev{'targetRevCml'}[$markBin]);
	$tmpPct=fstr($wev{'targetRevCml'}[$markBin]*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)"); $i0++;
        $worksheet0->write($i0,0,'Number of validated decoy',$titleFormat); $worksheet0->write($i0,1,$wev{'decoyRevCml'}[$markBin]);
	$tmpPct=fstr($wev{'decoyRevCml'}[$markBin]*100/scalar(@out));$worksheet0->write($i0,2,"\($tmpPct\%\)"); $i0++;
        for my $outfile (@out)
        {
                if ($$outInfor{$outfile}{'hitShown'} and $$outInfor{$outfile}{'WeightedEvalue'}[1]>=$wevCut)
                {
                        unless (defined($pp{'peptide'}{$$outInfor{$outfile}{'fullPeptide'}[1]})) {$pp{'peptide'}{$$outInfor{$outfile}{'fullPeptide'}[1]}=1;}
                        unless (defined($pp{'protein'}{$$outInfor{$outfile}{'protein'}[1]})) {$pp{'protein'}{$$outInfor{$outfile}{'protein'}[1]}=1;}
                }
        }
        $worksheet0->write($i0,0,'Number of unique peptide',$titleFormat); $worksheet0->write($i0,1,scalar(keys %{$pp{'peptide'}}));$i0++;
        $worksheet0->write($i0,0,'Number of unique protein',$titleFormat); $worksheet0->write($i0,1,scalar(keys %{$pp{'protein'}}));$i0++;


        $workbook->close();

}

sub fstr {
  my ($value,$precision) = @_;
  $precision ||= 2;
  my $s = sprintf("%.${precision}f", $value);
  $s =~ s/\.?0*$//;
  $s
}

sub _calculate_qValue
{
        my ($smallhash,$scoreType,$fdrType)=@_;

        my @array;
        for my $outhash (sort {$$a{$scoreType} <=> $$b{$scoreType}} values %$smallhash)
        {
                push(@array, $outhash);
        }

        my @dCount; $dCount[scalar(@array)]=0;
        for (my $i=scalar(@array)-1; $i>=0; $i--)
        {
                my $tmp=$array[$i];
                $dCount[$i]=$dCount[$i+1];
                if ( $$tmp{protein} =~ /Decoy/ ) { $dCount[$i]++; }
        }

        for (my $i=0; $i<scalar(@array); $i++ )
        {
                my $tmp=$array[$i];
                if ($fdrType eq 'total') { $$tmp{qvalue}=2*$dCount[$i]/(scalar(@array)-$i); }
                elsif ($fdrType eq 'target') { $$tmp{qvalue}=$dCount[$i]/(scalar(@array)-$i-$dCount[$i]); }
        }

        my $preMax=0;
        for (my $i=scalar(@array)-1; $i>=0; $i--)
        {
                my $tmp=$array[$i];
                if ( $$tmp{qvalue}<$preMax ) { $$tmp{qvalue}=$preMax; }
                else { $preMax=$$tmp{qvalue}; }
        }
}

sub buildNomodHash
{
        my ($psmHash,$pepHash,$scoreType)=@_;

        while ( my ($outfile, $hash) = each %$psmHash )
        {
                my $pep=$$hash{peptide};

                if ( !defined($$pepHash{$pep}) ) { $$pepHash{$pep}{totalScan}=0; }
                $$pepHash{$pep}{totalScan}++;

                if ( !defined($$pepHash{$pep}{$scoreType}) || $$pepHash{$pep}{$scoreType}<$$hash{$scoreType} )
                {
             #           $$pepHash{$pep}{outfile}=$outfile;
                        $$pepHash{$pep}{$scoreType}=$$hash{$scoreType};
            #            $$pepHash{$pep}{dCn}=$$hash{dCn};
                        $$pepHash{$pep}{protein}=$$hash{protein};
                }
        }
}


=head
sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}
=cut
1;
