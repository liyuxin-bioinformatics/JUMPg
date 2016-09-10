#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.4/programs/f

use strict;
use warnings;
#use SpoutParser;
#use pepXML_parser;
use idsum2::pepXML_parser;
#use idsum2::pepXML_parser;
use Cwd;
use Statistics::R;
#use Clone 'clone';

if (scalar(@ARGV)!=1)
{
	die "Usage: perl spectrumQC.pl jump_qc.params\n";
}

# initialization
my (%parahash,%acpID,%runhash,%goodLdaBin);
parse_params($ARGV[0],\%parahash);
my $prefix='qc_';
my $outputDir=getcwd;
chomp($outputDir);
$outputDir.="\/$prefix";
$outputDir.=$parahash{output_folder};
#print "$outputDir\n";
if (! -e $outputDir) { system("mkdir $outputDir"); }
system(qq(cp $ARGV[0] $outputDir));

# input file checking
#print "Checking input files ...\n";
inputFileCheck(\%parahash);
#print "All input files are found\n";

# check successful rate
my $successRate=checkSuccessfulRate($parahash{confident_IDtxt});
if ($successRate<10){
	$parahash{confident_IDtxt}=0;
}else {
	#print "Pass minimum successful rate requirement= $successRate\%\n";
	print "Pass minimum successful rate requirement (10\%)\n";
}

print "\nPrepare matrix for quality score calculation:\n";

open(TMP,">\.jump_qc_tmp");
if ($parahash{confident_IDtxt} ne '0')
{
#=head
# build %acpID
# 1) print ID_scanInfor.txt
# 2) build %acpID based on ID_scanInfor.txt
print "Reading PSM filtering result ($parahash{confident_IDtxt}) ...\n";
#printScanInfor($parahash{confident_IDtxt},$outputDir,\%acpID);
printScanInfor($parahash{confident_IDtxt},"$outputDir\/confident_ID_scanInfor.txt",\%acpID);

# for each fraction, generate:
# mzXML_infor: from mzXML
# spout.tab: from pepXML
# tagInfor: from tags
# tag.pks.acp.tab: spout.tab + tagInfor + mzXML_infor + %acpID
#foreach my $run (keys %{$parahash{'runs'}})
for my $run (sort {$a cmp $b} keys %{$parahash{'runs'}})
{
	my (%outhash);

	print "Processing fraction $run:\n";
	$run =~ s/\/$//;

	my $run_outputDir=$outputDir; my $runBase;
	if (1)
	{
		my @t=split /\//,$run;
		$t[$#t] =~ s/\.\d+$//;
		$run_outputDir.="\/$t[$#t]"; $runBase=$t[$#t];
		if (! -e $run_outputDir) { system("mkdir $run_outputDir"); }
		print TMP "$run_outputDir\n";
	}

	# spout.tab: from pepXML
	#my $pepxml=pepXML_parser->new();
	my $pepxml=idsum2::pepXML_parser->new();
	my $pepXML=$run;
	$pepXML.='.pepXML';
	print "Reading pepXML file ($pepXML) ...\n";
	$pepxml->pepXML2Hashes(\%parahash,\%outhash,$pepXML);
	#$pepxml->printPepXML(\%parahash,\%outhash,'test.pepXML',10);
	
	# tagInfor: from tags
	parseTags(\%outhash,$run);

	# mzXML_infor: from mzXML
	parseMzXML(\%outhash,$run,$run_outputDir);

	# add accepted information
	attachAccepted(\%outhash,\%acpID);
	
	# print tag.pks.acp.tab
	my $output=$run_outputDir;
	$output.="\/$runBase\.tag.pks.acp.tab";
	printPSMinfor(\%outhash,$output,'fraction');	

	# hash pass to %runhash
	$runhash{$runBase}=\%outhash;
	#$runhash{$runBase}=clone(\%outhash);

	# record $run~$runBase pairs
	$parahash{'runpaths'}{$runBase}=$run;
}

# combine all fractions and print 'all.tag.pks.acp.tab'
my $output="$outputDir\/all.tag.pks.acp.tab";
open(OUT,">$output");
print OUT "outfile\tprecursor_matched_peptides\tmatchedType\tscore\tdeltaScore\tpeptide\tmatchedTagSeq\tmatchedTagLength\ttagNum\ttopTagSeq\ttopTagLength\ttopTagEvalue\trun\tscan\tPPI\tcharge\trt\tpreInt\tpreMz\tms2PeaksCount\tms2Int\taccepted\n";
close OUT;
for my $run (sort {$a cmp $b} keys %runhash)
{
	#print "$run\n";
	printPSMinfor($runhash{$run},$output,'combined');
}
#=cut

# LDA analysis
print "\nPerforming LDA analysis ...\n";
rLDA_analysis($outputDir,$parahash{'PSM_recoveray_rate'});

# attach LDA information to %runhash
attachLDA(\%runhash,"$outputDir\/all.tag.pks.acp.qs.tab");

# build %goodLdaBin by parsing 'threshold_table.txt'
parseThresholdTable(\%goodLdaBin,"$outputDir\/threshold_table.txt",$parahash{'PSM_recoveray_rate'});
} # end of "if ($parahash{confident_IDtxt} ne '0')"
else
{
for my $run (sort {$a cmp $b} keys %{$parahash{'runs'}})
{
	my (%outhash);

	print "Processing fraction $run:\n";
	$run =~ s/\/$//;

	my $run_outputDir=$outputDir; my $runBase;
	if (1)
	{
		my @t=split /\//,$run;
		$t[$#t] =~ s/\.\d+$//;
		$run_outputDir.="\/$t[$#t]"; $runBase=$t[$#t];
		if (! -e $run_outputDir) { system("mkdir $run_outputDir"); }
		print TMP "$run_outputDir\n";
	}

	# spout.tab: from pepXML
	#my $pepxml=pepXML_parser->new();
	my $pepxml=idsum2::pepXML_parser->new();
	my $pepXML=$run;
	$pepXML.='.pepXML';
	print "Reading pepXML file ($pepXML) ...\n";
	$pepxml->pepXML2Hashes(\%parahash,\%outhash,$pepXML);

	# hash pass to %runhash
	$runhash{$runBase}=\%outhash;

	# record $run~$runBase pairs
	$parahash{'runpaths'}{$runBase}=$run;
}
} # end of 'else'

# re-build %acpID using 'accepted_IDtxt'
# 1) print ID_scanInfor.txt
# 2) build %acpID based on ID_scanInfor.txt
# 3) update %runshash{$run}{$out}{accepted} using new %acpID
undef %acpID;
#print "\nReading ID.txt ($parahash{accepted_IDtxt}) ...\n";
printScanInfor($parahash{accepted_IDtxt},"$outputDir\/accepted_ID_scanInfor.txt",\%acpID);
for my $run (keys %runhash)
{
	attachAccepted($runhash{$run},\%acpID);
}

# build %scanhash
# %scanhash{$run}{$scan}{accepted/highquality}
my (%scanhash);
buildScanhash(\%scanhash,\%runhash,\%goodLdaBin,($parahash{confident_IDtxt} eq '0')?0:1);
print "\n";
printSummary(\%scanhash);

# for each fraction:
# parse .dtas file and accept a dta file for futher analysis only if: 
# 1) PPI==1
# 2) accepted==0 (previously not accepted spectra)
# 3) with good ladscore (defined in %goodLdaBin)
print "\nFiltering low quality spectra ...\n";
filterDta(\%parahash,\%runhash,\%goodLdaBin,$outputDir,\%scanhash);
filterTag(\%parahash,\%runhash,\%goodLdaBin,$outputDir,\%scanhash);

#----------------------------------------------------------------------------
sub printSummary {
	my ($scanhash)=@_;

	my ($ttl,$acp,$hq,$lq);
	$ttl=$acp=$hq=$lq=0;

	foreach my $run (keys %{$scanhash}) {
		# total MS2 counts
		$ttl+=scalar(keys %{$scanhash{$run}}); 
		# accepted MS2 counts
		map {$acp++ if ($$scanhash{$run}{$_}{accepted}==1)} keys %{$scanhash{$run}};
		# high quality but not accpeted MS2 counts
		map {$hq++ if ($$scanhash{$run}{$_}{accepted}==0 
			&& $$scanhash{$run}{$_}{highquality}==1)} keys %{$scanhash{$run}};
		# low quality and not accpeted MS2 counts
		map {$lq++ if ($$scanhash{$run}{$_}{accepted}==0 
			&& $$scanhash{$run}{$_}{highquality}==0)} keys %{$scanhash{$run}};
	}

	print "Starting with $ttl MS2 scans:\n";
	print "  There are $acp scans accepted (based on PSM filtering results)\n";
	print "  Of the ",$ttl-$acp," remaining scans:\n";
	print "    Low quality spectra (to be removed): $lq\n";
	print "    High quality spectra (to be kept): $hq\n";
}

sub inputFileCheck
{
	my ($parahash)=@_;

	# ID.txt
	if ($$parahash{confident_IDtxt} ne '0' and !-e $$parahash{confident_IDtxt}) { die "ID.txt $$parahash{confident_IDtxt} not found!!!\n"; }
	if (!-e $$parahash{accepted_IDtxt}) { die "ID.txt $$parahash{accepted_IDtxt} not found!!!\n"; }

	# foreach fraction, check:
	# mzXML, pepXML, tags, dtas
	foreach my $run (keys %{$$parahash{'runs'}})
	{
		$run =~ s/\/$//;

		# pepXML
		my $file=$run;
		$file.='.pepXML';
		if (!-e $file) { die "$file  not found!!!\n"; }

		# mzXML
		if ($$parahash{confident_IDtxt} ne '0')
		{
			$file=$run;
			$file.='.mzXML';
			if (-e $file) { }
			else
			{
				# ad_pl07.1.mzXML => ad_pl07.mzXML
				$file =~ s/\.\d+\.mzXML/\.mzXML/;
				if (!-e $file) { die "$file  not found!!!\n"; }
			}
		}

		# tags
		$file=$run;
		$file.='.tags';
		if (!-e $file) { die "$file  not found!!!\n"; }

		#dtas
		$file=$run;
		$file.='.dtas';
		if (!-e $file) { die "$file  not found!!!\n"; }
	}
}

sub filterTag
{
	my ($parahash,$runhash,$goodLdaBin,$outputDir,$scanhash)=@_;

	my $lda=($$parahash{confident_IDtxt} eq '0')?0:1;
	foreach my $run (keys %{$runhash})
	{
		# input path setup
		my $runpath=$$parahash{'runpaths'}{$run};
		my $tags=$runpath;
		$tags.='.tags';

		# output path setup
		my $run_outputDir=$outputDir;
		$run_outputDir.="\/$run\/$run\.filtered\.tags";

		# parse and filter .tags
		#print "Reading tags file ($tags) ...\n";
		open(IN,$tags) || die "cannot open $tags";
		open(OUT,">$run_outputDir");
		my @array; my $out='';
		while(<IN>)
		{
			chomp;
			if (/\.tag/)
			{
				if ($out ne '')
				{
					if ( accepctSpectrum($lda,$out,$$runhash{$run},$goodLdaBin,$scanhash) )
					{
						for (my $i=0; $i<=$#array;$i++)
						{
							my $line=$array[$i];
							print OUT "$line\n";
						}
					}
					undef @array;
				}

				push @array, $_;
				$out=$_;
				$out =~ s/\.tag//;
			}
			else
			{
				push @array, $_;
			}
		}

		# the last tag file
		if ($out ne '')
		{
			if ( accepctSpectrum($lda,$out,$$runhash{$run},$goodLdaBin,$scanhash) )
			{
				for (my $i=0; $i<=$#array;$i++)
				{
					#my $line=shift @array;
					my $line=$array[$i];
					print OUT "$line\n";
				}
			}
		}

		close IN;
		close OUT;
	}
}

sub buildScanhash
{
	my ($scanhash,$runhash,$goodLdaBin,$lda)=@_;
	foreach my $run (keys %{$runhash})
	{
		foreach my $out (keys %{$$runhash{$run}})
		{
			my ($runname,$scan,$ppi,$z)=split /\./,$out;

			if ( $$runhash{$run}{$out}{accepted}==1 ) { 
				$$scanhash{$runname}{$scan}{accepted}=1; 
			# this step is tricky because one single scan could have more than 
			# one outfiles, and maximally only one is correct; check whether
			# this value if exist to avoid erase an already accepted scan
			} elsif (!defined($$scanhash{$runname}{$scan}{accepted})) { 
				$$scanhash{$runname}{$scan}{accepted}=0; 
			}

			if ($lda)
			{
				if (!defined($$runhash{$run}{$out}{ladBin})) { $$scanhash{$runname}{$scan}{highquality}=0; }
				else {
					my $ladBin=$$runhash{$run}{$out}{ladBin};
					if (defined($$goodLdaBin{$ladBin})) { $$scanhash{$runname}{$scan}{highquality}=1; }
					else { $$scanhash{$runname}{$scan}{highquality}=0; }
				}
			}
			else { $$scanhash{$runname}{$scan}{highquality}=1; }
		}
	}
}

sub accepctSpectrum
{
	my ($lda,$out,$hash,$goodLdaBin,$scanhash)=@_;
	my ($runname,$scan,$ppi,$z)=split /\./,$out;
	#if ($ppi>1) { return 0; }  # only keep 1st ppi

	if (defined($$scanhash{$runname}{$scan}))
	{
		return 0 if ($$scanhash{$runname}{$scan}{accepted}==1 ); # rm accepted PSMs
		if ($lda)
		{
			return 0 if ($$scanhash{$runname}{$scan}{highquality}==0); # only keep good quality spectra
		}
		return 1;
	}
	#else { die "$out in dtas/tags but not found in scanhash!!!";  }
	else { return 1;  } # dta/tag missed in the search

=head
	if (defined($$hash{$out}))
	{
		return 0 if ($$hash{$out}{accepted}==1); # rm accepted PSMs
		if ($lda)
		{
			my $ladBin=$$hash{$out}{ladBin};
			return 0 if (!defined($$goodLdaBin{$ladBin})); # only keep good quality spectra
		}
		return 1;
	}
	else { die "$out in dtas/tags but not found in runhash!!!";  }
=cut
}

sub filterDta
{
	my ($parahash,$runhash,$goodLdaBin,$outputDir,$scanhash)=@_;

	my $lda=($$parahash{confident_IDtxt} eq '0')?0:1;
	foreach my $run (keys %{$runhash})
	{
		# input path setup
		my $runpath=$$parahash{'runpaths'}{$run};
		my $dtas=$runpath;
		$dtas.='.dtas';

		# output path setup
		my $run_outputDir=$outputDir;
		$run_outputDir.="\/$run\/$run\.filtered\.dtas";

		# parse and filter .dtas
		#print "Reading dtas file ($dtas) ...\n";
		open(IN,$dtas) || die "cannot open $dtas";
		open(OUT,">$run_outputDir");
		while(<IN>)
		{
			my $line0=$_;

			chomp;
			my @t=split / /,$_;
			my $out=$t[0];
			$out =~ s/\.dta//;

			my $line1=<IN>;
			my $line2=<IN>;

			if (accepctSpectrum($lda,$out,$$runhash{$run},$goodLdaBin,$scanhash))
			{
				print OUT "$line0$line1$line2";
			}
		}
		close IN;
		close OUT;
	}
}

sub parseThresholdTable
{
	my ($goodLdaBin,$input,$PSM_recoveray_rate)=@_;
	my $idlostCut=1-$PSM_recoveray_rate/100;

	open(IN,$input) || die "Cannot open $input...\n";
	while(<IN>)
	{
		next if (/^fail/);
		chomp;
		my @t=split /\t/,$_;
		my ($ladBin,$IDlostPct)=($t[0],$t[$#t]);
		if ($IDlostPct>$idlostCut) { $$goodLdaBin{$ladBin}=''; }
	}
	close IN;
}

sub attachLDA
{
	my ($runhash,$input)=@_;

	open(IN,$input) || die "Cannot open $input...\n";
	while(<IN>)
	{
		next if (/^outfile/);
		chomp;
		my @t=split /\t/,$_;
		my ($out,$run,$ldaScore,$ladBin)=($t[0],$t[12],$t[23],$t[24]);
		if (defined($$runhash{$run}{$out}))
		{
			$$runhash{$run}{$out}{ldaScore}=$ldaScore;
			$$runhash{$run}{$out}{ladBin}=$ladBin;
		}
		else { die "$run $out in all.tag.pks.acp.qs.tab but not found in runhash!\n"; }
	}
	close IN;
}

sub rLDA_analysis
{
	my ($outputDir,$PSM_recoveray_rate)=@_;

	my $inputFileBaseName = 'all.tag.pks.acp.tab';
	my $R = Statistics::R -> new();
	$R -> set('inputFile', $inputFileBaseName);
	$R -> set('inputDirectory', $outputDir);
	$R -> set('recoveryPct',$PSM_recoveray_rate);
	$R -> run(q`setwd(inputDirectory);
	suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
	library(MASS);
	## Data loading
	tb=read.table(inputFile,head=T,sep="\t")
	tb=tb[tb$PPI==1,]
	tb$topTagEvalue[is.na(tb$topTagEvalue)]=0
	# specify training set
	tb$trainingSet=rep('notUsed',nrow(tb))
	tb$trainingSet[tb$accepted==1]='positive'
	tb$trainingSet[tb$matchedType=='unmatched']='negative'
	# Generation of a data matrix for LDA analysis
	subData=tb[tb$trainingSet!='notUsed',]
	X = as.data.frame(cbind(subData$topTagLength,subData$topTagEvalue, sqrt(subData$ms2PeaksCount), log10(subData$ms2Int)))
	group = subData$trainingSet
	# LAD analysis
	calclda <- function(variables,loadings)
	{
		as.data.frame(variables)
		numsamples <- nrow(variables)
		ld <- numeric(numsamples)
		numvariables <- length(variables)
		for (i in 1:numsamples)
		{
			valuei <- 0
			for (j in 1:numvariables)
			{
				valueij <- variables[i,j]
				loadingj <- loadings[j]
				valuei <- valuei + (valueij * loadingj)
			}
			ld[i] <- valuei
		}
		return(ld)
	}
	if (length(table(group))>=2 & table(group)[1]>1 & table(group)[2]>1) 
	{
		# LDA modeling
		ldaModel = lda(group ~ ., X);
		scores = as.numeric(predict(ldaModel, X)$x);
		# apply same model to the whole dataset
		m=as.data.frame(cbind(tb$topTagLength,tb$topTagEvalue, sqrt(tb$ms2PeaksCount), log10(tb$ms2Int)))
		tb$ldaScore=calclda(m[,1:4], ldaModel$scaling[,1])
		# binning for ldaScore
		binN=100
		intv=(max(tb$ldaScore)-min(tb$ldaScore))/binN
		tb$binldaScore=floor((tb$ldaScore-min(tb$ldaScore))/intv)
		mp1=table(tb$binldaScore,tb$accepted)
		mp=matrix(nrow=nrow(mp1),ncol=6)
		colnames(mp)=c('fail','accept','successRate','identifiableSpectra','spectraUsed4search','IDlost')
		rownames(mp)=rownames(mp1)
		mp[,1:2]=mp1
		mp=as.data.frame(mp)
		mp$successRate=mp$accept/(mp$fail+mp$accept)
		#mp$identifiableSpectra=mp$fail*mp$successRate # based on successRate*fail
		mp$identifiableSpectra=mp$accept	# same as accept
		#spectraUsed4search
		#mp[,5]=mp[,1] # only fail
		mp[,5]=mp[,1]+mp[,2] # both fail and accept
		for (i in (nrow(mp)-1):1)
		{
			#mp[i,5]=mp[i,1]+mp[i+1,5] # only fail
			mp[i,5]=mp[i,1]+mp[i,2]+mp[i+1,5] # both fail and accept
		}
		#IDlost
		mp[,6]=mp$identifiableSpectra
		for (i in 2:nrow(mp))
		{
			mp[i,6]=mp[i-1,6]+mp[i,4]
		}
		#mp$spectraUsed4searchPct=mp$spectraUsed4search/sum(mp[,1]) # only fail
		mp$spectraUsed4searchPct=mp$spectraUsed4search/(sum(mp[,1])+sum(mp[,2])) # both fail and accept
		mp$IDlostPct=mp$IDlost/sum(mp[,4])
		bin.lda=mp
		spectraSkip=100-tail(bin.lda$spectraUsed4searchPct[bin.lda$IDlostPct<(1-recoveryPct/100)]*100,1)
		# figure
		w=2
		par(las=1)
		pdf(file=paste(inputDirectory,'/PSM_recovery_vs_spctraUsed.pdf',sep=''),width=5,height=5)
		mp=bin.lda
		plot((1-mp$spectraUsed4searchPct)*100,(1-mp$IDlostPct)*100,xlab='% spectra skipped for database search',ylab='% expected identified PSMs',col='red',type='l',lwd=w)
		abline(h=recoveryPct,col='grey',lty=2)
		abline(v=spectraSkip,col='grey',lty=2)
		dev.off()
		write.table(tb,file=paste(inputDirectory,'/all.tag.pks.acp.qs.tab',sep=''),quote=F,sep="\t",row.names=F)
		write.table(bin.lda,file=paste(inputDirectory,'/threshold_table.txt',sep=''),quote=F,sep="\t")
	}`);
	$R -> stop();
	print "LDA modeling in R has been finished\n";
}

sub attachAccepted
{
	my ($outhash,$acpID)=@_;

	foreach my $out (keys %{$outhash})
	{
		if (defined($$acpID{$out}))
		{
			$$outhash{$out}{accepted}=1;
		}
		else
		{
			$$outhash{$out}{accepted}=0;
		}
	}
}

sub parseMzXML
{
	my ($outhash,$run,$run_outputDir)=@_;

	my $mzXML=$run;
	$mzXML.='.mzXML';
	if (-e $mzXML) { }
	else
	{
		# ad_pl07.1.mzXML => ad_pl07.mzXML
		$mzXML =~ s/\.\d+\.mzXML/\.mzXML/;
		if (!-e $mzXML) { die "$mzXML  not found!!!\n"; }
	}

	print "Reading mzXML file ($mzXML) ...\n";

	# parse mzXML, build %scans
	my %scans;
	open(IN,$mzXML)  or die "Cannot open $mzXML!!!\n";
	my ($scan,$ms2,$rt,$peaksCount,$ms2Int);
	while(<IN>)
	{
        	chomp;
	        if (/\<scan num=\"(\d+)\"/) { $scan=$1;  }
        	elsif (/msLevel=\"2\"/) { $ms2=1; }
	        elsif (/retentionTime=\"(.*?)\"/) { $rt=$1; }
        	elsif (/peaksCount=\"(.*?)\"/) { $peaksCount=$1; }
	        elsif (/basePeakIntensity=\"(.*?)\"/)  {  $ms2Int=$1; }
        	elsif (/\<precursorMz precursorIntensity=\"(.*?)\" activationMethod=".*?" \>(.*?)\<\/precursorMz\>/)
	        {
        	        my ($preInt,$mz)=($1,$2);

                	next unless ($ms2==1);

	                $scans{$scan}{rt}=$rt;
        	        $scans{$scan}{preInt}=$preInt;
                	$scans{$scan}{mz}=$mz;
	                $scans{$scan}{peaksCount}=$peaksCount;
        	        $scans{$scan}{ms2Int}=$ms2Int;

                	$ms2=0;
	        }
	}
	close IN;

	# print mzXML_infor
	if (1)
	{
		my @t=split /\//,$run;
		my $runname=$t[$#t];
		open(OUT,">$run_outputDir\/$runname\.mzXMLINfor.txt");
		print OUT "fraction\tscan\trt\tpreInt\tpreMz\tms2PeaksCount\tms2Int\n";
		for my $scan (sort {$a<=>$b}  keys %scans)
		{
			print OUT "$runname\t$scan\t$scans{$scan}{rt}\t$scans{$scan}{preInt}\t$scans{$scan}{mz}\t$scans{$scan}{peaksCount}\t$scans{$scan}{ms2Int}\n";
		}
		close OUT
	}

	# attach mzXML_infor to %outhash
	foreach my $out (keys %{$outhash})
	{
		# ad_pl07.156438.1.2
		my ($runname,$scan,$ppi,$z)=split /\./,$out;
		if (defined($scans{$scan}))
		{
			$$outhash{$out}{rt}=$scans{$scan}{rt};
			$$outhash{$out}{preInt}=$scans{$scan}{preInt};
			$$outhash{$out}{preMz}=$scans{$scan}{mz};
			$$outhash{$out}{ms2PeaksCount}=$scans{$scan}{peaksCount};
			$$outhash{$out}{ms2Int}=$scans{$scan}{ms2Int};
		}
		else {die "$runname,$scan,$ppi,$z not exsit in mzXML infor\n";}
	}
}

sub printPSMinfor
{
	my ($hash,$output,$printType)=@_;

	if ($printType eq 'fraction')
	{
		open(OUT,">$output");
		print OUT "outfile\tprecursor_matched_peptides\tmatchedType\tscore\tdeltaScore\tpeptide\tmatchedTagSeq\tmatchedTagLength\ttagNum\ttopTagSeq\ttopTagLength\ttopTagEvalue\trun\tscan\tPPI\tcharge\trt\tpreInt\tpreMz\tms2PeaksCount\tms2Int\taccepted\n";
	}
	elsif ($printType eq 'combined')
	{
		open(OUT,">>$output");
	}
	else { die "Unexpected printType: $printType\n"; }

	foreach my $out (keys %$hash)
	{
		my ($score,$deltaScore)=('n/a','n/a');
		if (defined($$hash{$out}{WeightedEvalue}[1])) { $score=$$hash{$out}{WeightedEvalue}[1]; $deltaScore=$$hash{$out}{deltacn}[1]; }
		elsif (defined($$hash{$out}{xcorr}[1])) { $score=$$hash{$out}{xcorr}[1]; $deltaScore=$$hash{$out}{deltacn}[1]; }
		my $precursor_neutral_mass=$$hash{$out}{precursor_neutral_mass};
		my $precursor_matched_peptides=$$hash{$out}{precursor_matched_peptides};

		my ($calc_neutral_pep_mass,$peptide,$protein,$TD,$matchedType,$TagSeq,$tagL)=('n/a','n/a','n/a','n/a','unmatched','',0);
		if ($score ne 'n/a') # with matches
		{
			$calc_neutral_pep_mass=$$hash{$out}{calc_neutral_pep_mass}[1];
			$peptide=$$hash{$out}{peptide}[1];
			$protein=$$hash{$out}{protein}[1];
			$TagSeq=$$hash{$out}{TagSeq}[1];
			if ($TagSeq ne 'N/A') { $tagL=length($TagSeq); }
			if ($protein =~ m/Decoy/) {$TD=$matchedType='decoy'; $protein =~ s/\#\#//;} else {$TD=$matchedType='target'; }
		}

		#my ($tagN,$toptagseq,$toptagE)=(0,'n/a',0);
		my $tagN=(defined($$hash{$out}{tagN}))?$$hash{$out}{tagN}:0;
		my $toptagseq=(defined($$hash{$out}{toptagseq}))?$$hash{$out}{toptagseq}:'n/a';
		my $toptagL=($toptagseq eq 'n/a')?0:length($toptagseq);
		my $toptagE=(defined($$hash{$out}{toptagE}))?$$hash{$out}{toptagE}:0;
		#my $toptagseq=$$hash{$out}{toptagseq};
		#my $toptagE=$$hash{$out}{toptagE};

		# ad_pl07.156438.1.2
		my ($runname,$scan,$ppi,$z)=split /\./,$out;
		my $rt=$$hash{$out}{rt};
		my $preInt=$$hash{$out}{preInt};
		my $preMz=$$hash{$out}{preMz};
		my $ms2PeaksCount=$$hash{$out}{ms2PeaksCount};
		my $ms2Int=$$hash{$out}{ms2Int};
		my $accepted=$$hash{$out}{accepted};

		#print OUT "$out\t$score\t$deltaScore\t$precursor_neutral_mass\t$calc_neutral_pep_mass\t$peptide\t$protein\t$TD\n";
		print OUT "$out\t$precursor_matched_peptides\t$matchedType\t$score\t$deltaScore\t$peptide\t$TagSeq\t$tagL\t$tagN\t$toptagseq\t$toptagL\t$toptagE\t$runname\t$scan\t$ppi\t$z\t$rt\t$preInt\t$preMz\t$ms2PeaksCount\t$ms2Int\t$accepted\n";
	}

	close OUT;
}

sub parseTags
{
	my ($outhash,$run)=@_;

	my $tagsfile=$run;
	$tagsfile.='.tags';
	print "Reading tags file ($tagsfile) ...\n";

	open(IN,$tagsfile)  or die "Cannot open $tagsfile!!!\n";
	my ($toptag,$tagN,$out)=(0,0,'');
	while(<IN>)
	{
		chomp;
		if (/tag/)
		{
			if ($out ne '') {$$outhash{$out}{tagN}=$tagN;}

			# ad_pl07.100.1.2.tag
			s/\.tag//;
			$out=$_;
			if (!defined($$outhash{$out})) { die "$out in tags but not found in outhash!!!"; }
			$toptag=1; $tagN=0;
		}
		else
		{
			if ($toptag)
			{
				my @t=split /\t/,$_;
				$$outhash{$out}{toptagseq}=$t[1];
				$$outhash{$out}{toptagE}=$t[5];
				$toptag=0;
			}
			$tagN++;
		}
	}
	if ($out ne '') {$$outhash{$out}{tagN}=$tagN;}
	close IN;
}

sub printScanInfor
{
	my ($idtxt,$output,$hash)=@_;
	#my %hash;

	open(IN,$idtxt)  or die "Cannot open $idtxt!!!\n";
	while(<IN>)
	{
		if (/^Database/) {  next; }
		if (/^Peptide/) {  next; }

		chomp;
		my ($Peptide,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos)=split(/;/,$_);

		my @t=split /\//,$Outfile;
		$Outfile=$t[$#t];
		$Outfile =~ s/\.spout$//;
		$Outfile =~ s/\.out$//;

		next if (defined($$hash{$Outfile}));

		$$hash{$Outfile}{XCorr}=$XCorr;
		$$hash{$Outfile}{dCn}=$dCn;
		$$hash{$Outfile}{ppm}=$ppm;
		$$hash{$Outfile}{Peptide}=$Peptide;
		$$hash{$Outfile}{Protein}=$Protein;
	}
	close IN;

	#my $output=$outputDir;
	#$output.="\/ID_scanInfor.txt";
	open(OUT,">$output");
	print OUT "outfile\tXCorr\tdCn\tppm\tpeptide\tprotein\n";
	foreach my $Outfile (keys %{$hash})
	{
		print OUT "$Outfile\t$$hash{$Outfile}{XCorr}\t$$hash{$Outfile}{dCn}\t$$hash{$Outfile}{ppm}\t$$hash{$Outfile}{Peptide}\t$$hash{$Outfile}{Protein}\n";
	}
	close OUT;
}

sub parse_params
{
        my ($par,$parahash)=@_;
        open(IN,$par);
        while (<IN>)
        {
                s/^\s+//;
                next if (/^#/);
                chomp;
                next unless (/ = /);
                my ($key,$value) = split(' = ',$_);
                if ($key =~ /^run(\d+)/)
                {
                        $$parahash{'runs'}{$value}='';
                }
                else
                {
                        $$parahash{$key}=$value;
                }
        }
        close IN;
}

sub checkSuccessfulRate {
	my ($idtxt)=@_;

	my $r=-1;
	if (-e $idtxt) {
		$idtxt =~ s/ID.txt$/jump_f.log/;
		unless (-e $idtxt){
			return $r;
		}

		open(IN, $idtxt) || die "cannot open $idtxt\n";
		while (<IN>) {
			chomp;
			next unless (/Originally starting with (\d+) MS2 scans/);
			my $total=$1;
			my $line=<IN>;
			$line =~ m/Accepting (\d+) outfiles: (\d+) targets and (\d+) decoys /;
			my $acpT=$2;

			$r=$acpT*100/$total;
			return $r;
		}
		close IN;
	}else {
		return $r;
	}
}
