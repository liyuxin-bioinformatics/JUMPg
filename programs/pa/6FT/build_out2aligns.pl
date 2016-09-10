#!/usr/bin/perl -w
use strict;
use Storable;

if (scalar(@ARGV)<=5)
{
	die "usage: perl build_out2aligns.pl reads.fas out2read.hash pep2readPos.txt output tab_files... [e.g. vsCDS.tab vsTranscriptome.tab ...]\n";
}

#read fas, build %reads:
#%reads{$readID}{'length'}=$read_length
#               {level: 1/2/3/etc.}{score/gene/align_length/mismatch/gap/qStart/qEnd}
my (%reads);
open(IN,"$ARGV[0]");
while(<IN>)
{
        chomp;
        my $id=$_;
        $id=substr($id,1,length($id)-1);

        my $seq=<IN>;
        chomp($seq);

	$reads{$id}{'length'}=length($seq);
}
close IN;

#read tab, build %reads
# %reads{$readID}{$level}{score/gene/align_length/mismatch/gap/qStart/qEnd/sStart/sEnd/strand}
my $input_offset=4;
for (my $i=$input_offset; $i<scalar(@ARGV); $i++)
{
	parseTabFile($ARGV[$i],\%reads,$i-$input_offset+1);
}

#load hash
my $out2read=retrieve($ARGV[1]);

#build out2align
# %out2align{$outfile}{score/level/readID/gene/align_length/mismatch/gap/qStart/qEnd/length}
my %out2align=build_out2align(\%reads,$out2read,scalar(@ARGV)-$input_offset);

# read pep2readPos.txt; build %%pep2readPos{$pep}{$readID}{frame/pos/strand}
my (%pep2readPos);
open(IN,"$ARGV[2]");
while(<IN>)
{
	chomp;
	next if (/^peptide/);
	my ($pep,$readID,$pos,$frame)=split /\t/,$_;
	my ($pStrand,$offset)=split //,$frame;
	$pep2readPos{$pep}{$readID}{frame}=$frame;
	$pep2readPos{$pep}{$readID}{pos}=$pos;
	$pep2readPos{$pep}{$readID}{pStrand}=$pStrand;#print "$pStrand,";

	my ($start,$end)=split /\-/,$pos;
	if ($pStrand eq 'F')
	{
		$pep2readPos{$pep}{$readID}{pStart}=$start;
		$pep2readPos{$pep}{$readID}{pEnd}=$end;
	}
	elsif ($pStrand eq 'R')
	{
		my $readLength=$reads{$readID}{'length'};
		$pep2readPos{$pep}{$readID}{pStart}=$readLength-$end;
		$pep2readPos{$pep}{$readID}{pEnd}=$readLength-$start;
	}
	else { die "unexpected pStrand: $pStrand\n"; }
}
close IN;

# calculate peptide 2 subject position: pep2subPos
pep2subPos(\%out2align,\%pep2readPos,$out2read);

#output
open(OUT,">$ARGV[3]");
open(RM,">unmapped_outfiles");
print OUT "outfile\tXCorr\tdCn\tpeptide\tscanCount\treadCount\treadID\tucscID\tmismatch\tgap\treadLength\talignLength\tqStart\tqEnd\tlevels\tsStart\tsEnd\tsStrand\tb5\tb3\tstrand\tframe\n";
foreach my $out (keys %out2align)
{
	my $hash=\%{$out2align{$out}};
	if (defined($$hash{'score'}) and $$hash{'score'}){
		#print OUT "$out\t$$hash{readID}\t$$hash{gene}\t$$hash{mismatch}\t$$hash{gap}\t$$hash{length}\t$$hash{align_length}\t$$hash{qStart}\t$$hash{qEnd}\n";
		print OUT "$out\t$$out2read{$out}{XCorr}\t$$out2read{$out}{dCn}\t$$out2read{$out}{peptide}\t$$out2read{$out}{scanCount}\t",scalar(keys %{$$out2read{$out}{readID}}),"\t$$hash{readID}\t$$hash{gene}\t$$hash{mismatch}\t$$hash{gap}\t$$hash{length}\t$$hash{align_length}\t$$hash{qStart}\t$$hash{qEnd}\t$$hash{level}\t$$hash{sStart}\t$$hash{sEnd}\t$$hash{strand}\t$$hash{b5}\t$$hash{b3}\t$$hash{finalStrand}\t$$hash{frame}\n";
	}
	else 
	{
		print RM "$out\n";
		print OUT "$out\t$$out2read{$out}{XCorr}\t$$out2read{$out}{dCn}\t$$out2read{$out}{peptide}\t$$out2read{$out}{scanCount}\t",scalar(keys %{$$out2read{$out}{readID}}),"\t\t\t\t\t\t\t\t\t0\t\t\t\t\t\t\t\n";
		
	}
}
close OUT;
close RM;
#-------------------------------------------------------------------
sub pep2subPos
{
	my ($out2align,$pep2readPos,$out2read)=@_;
	foreach my $out (keys %{$out2align})
	{
		my $hash=\%{$$out2align{$out}};
		next if (!defined($$hash{'score'}) or $$hash{'score'}==0); # unmapped reads
		if ($$hash{gap}>0) # with gap
		{
			$$hash{b5}='';
			$$hash{b3}='';
			$$hash{finalStrand}='';
			$$hash{frame}='';
			next;
		}

		my $pep=$$out2read{$out}{peptide};
		$pep =~ s/[\*\@\$\^]//g;
		#my $nomod=$pep;
		my $readID=$$hash{readID};

		my ($pStart,$pEnd,$pStrand)=($$pep2readPos{$pep}{$readID}{pStart},$$pep2readPos{$pep}{$readID}{pEnd},$$pep2readPos{$pep}{$readID}{pStrand});
		if (!defined($pStrand)) { die "pStrand not defined: $out, $pep\n"; }

		# determine final strand
		my $strand='';
		if ( $$hash{strand} eq '+' and $pStrand eq 'F' ) { $strand='+'; }
		elsif ( $$hash{strand} eq '+' and $pStrand eq 'R' ) { $strand='-'; }
		elsif ( $$hash{strand} eq '-' and $pStrand eq 'F' ) { $strand='-'; }
		elsif ( $$hash{strand} eq '-' and $pStrand eq 'R' ) { $strand='+'; }

		# determine position: $b5, $b3
		my ($b5, $b3);
		if ( $$hash{strand} eq '+' )
		{
			$b5 = ($$hash{sStart}-1) + ($$hash{qStart}-1) + $pStart;
			$b3 = ($$hash{sStart}-1) + ($$hash{qStart}-1) + $pEnd;
		}
		else
		{
			$b5 = ($$hash{sStart}-1) + ($$hash{length}-$$hash{qEnd}) + ($$hash{length}-$pEnd);
			$b3 = ($$hash{sStart}-1) + ($$hash{length}-$$hash{qEnd}) + ($$hash{length}-$pStart);
		}

		$$hash{b5}=$b5;
		$$hash{b3}=$b3;
		$$hash{finalStrand}=$strand;

		# peptide frame:
		$$hash{frame}=$b5 % 3;
	}
}

sub parseTabFile
{
	my ($input,$reads,$level)=@_;

	open(IN,"$input");
	while(<IN>)
	{
		next if (/^#/);
		chomp;
		#4_1108_11983_108489-1   hg19_knownGene_uc001api.2       98.89   90      1       0       1       90      21      110     2e-42    170
		my @tmp=split(/\t/,$_);
		unless (scalar(@tmp)==12) {die scalar(@tmp),"\n$_\n";}

		if (defined($$reads{$tmp[0]}))
		{
			my $replace;
			if (defined($$reads{$tmp[0]}{$level}{'score'}) and $$reads{$tmp[0]}{$level}{'score'}>$tmp[11]+$tmp[3])
			{
				$replace=0;
			}
			else { $replace=1; }

			if ($replace)
			{
				$$reads{$tmp[0]}{$level}{'score'}=$tmp[11]+$tmp[3]; # bitScore + alignment_length
				#$tmp[1] =~ s/hg19_knownGene_//; 
				my @t=split /\_/,$tmp[1]; $tmp[1]=$t[$#t];
				$$reads{$tmp[0]}{$level}{'gene'}=$tmp[1];
				$$reads{$tmp[0]}{$level}{'align_length'}=$tmp[3];
				$$reads{$tmp[0]}{$level}{'mismatch'}=$tmp[4];
				$$reads{$tmp[0]}{$level}{'gap'}=$tmp[5];
				$$reads{$tmp[0]}{$level}{'qStart'}=$tmp[6];
				$$reads{$tmp[0]}{$level}{'qEnd'}=$tmp[7];

				if ( $tmp[8]<$tmp[9] )
				{
					$$reads{$tmp[0]}{$level}{'sStart'}=$tmp[8];
					$$reads{$tmp[0]}{$level}{'sEnd'}=$tmp[9];
					$$reads{$tmp[0]}{$level}{'strand'}='+';
				}
				else
				{
					$$reads{$tmp[0]}{$level}{'sStart'}=$tmp[9];
					$$reads{$tmp[0]}{$level}{'sEnd'}=$tmp[8];
					$$reads{$tmp[0]}{$level}{'strand'}='-';
				}
#				$$reads{$tmp[0]}{'level'}=$level;
#				$$reads{$tmp[0]}{'validated'}=isValidated(\%{$$reads{$tmp[0]}});
			}
		}
		else
		{
			die "not existed read: $tmp[0]\n";
		}
	}
	close IN;
}

sub build_out2align
{
	my ($reads,$out2read,$levels)=@_;

	my (%out2align);
	foreach my $out (keys %{$out2read})
	{
		my $maxS=0;my $maxID='';my $maxLevel=0;
		# check for each level: try to identify the 'best' alignment with the 'top' level
		for (my $i=1; $i<=$levels; $i++)
		{
			$maxS=0;$maxLevel=0; $maxID='';
			# find the best alignment within this level
			foreach my $id (keys %{$$out2read{$out}{readID}})
			{
				if (defined($$reads{$id}{$i}))
				{
					if ( defined($$reads{$id}{$i}{'score'}) and 
						$$reads{$id}{$i}{'score'}>$maxS and
						isValidated(\%{$$reads{$id}{$i}},$$reads{$id}{'length'}))
					{ $maxS=$$reads{$id}{$i}{'score'}; $maxID=$id; }
				}
				#else {die "not exist read: $id (level: $i)\n";}
			}
			# otherwise: go to the next lower level
			#if ($maxS and isValidated(\%{$$reads{$maxID}{$i}},$$reads{$maxID}{'length'})) { $maxLevel=$i; last; }
			if ($maxS) { $maxLevel=$i; last; }
		}

		if ($maxLevel)
		{
			$out2align{$out}{'score'}=$maxS;
			$out2align{$out}{'level'}=$maxLevel;
			$out2align{$out}{'readID'}=$maxID;
			$out2align{$out}{'gene'}=$$reads{$maxID}{$maxLevel}{'gene'};
			$out2align{$out}{'align_length'}=$$reads{$maxID}{$maxLevel}{'align_length'};
			$out2align{$out}{'mismatch'}=$$reads{$maxID}{$maxLevel}{'mismatch'};
			$out2align{$out}{'gap'}=$$reads{$maxID}{$maxLevel}{'gap'};
			$out2align{$out}{'qStart'}=$$reads{$maxID}{$maxLevel}{'qStart'};
			$out2align{$out}{'qEnd'}=$$reads{$maxID}{$maxLevel}{'qEnd'};
			$out2align{$out}{'length'}=$$reads{$maxID}{'length'};

			$out2align{$out}{'sStart'}=$$reads{$maxID}{$maxLevel}{'sStart'};
			$out2align{$out}{'sEnd'}=$$reads{$maxID}{$maxLevel}{'sEnd'};
			$out2align{$out}{'strand'}=$$reads{$maxID}{$maxLevel}{'strand'};
		}
		else { $out2align{$out}{'score'}=0; }
	}

	return %out2align;
}

sub isValidated
{
	my ($align,$readLength)=@_;

	if ($$align{'qStart'}==1 and $$align{'qEnd'}==$readLength ) { return 1; }
	else { return 0; }
}
