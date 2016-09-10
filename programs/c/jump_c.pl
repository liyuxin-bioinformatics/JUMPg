#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

if (scalar(@ARGV)!=1)
{
        die "Usage: perl jump_g.pl jump_g.params\n";
}

# code path
my $code_path="/home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/c/customizedDB";

#initialization
my (%parahash);
parse_params($ARGV[0],\%parahash);

# output folder
#my $prefix='gn_';
my $prefix='';
my $outputDir=getcwd;
chomp($outputDir);
$outputDir.="\/$prefix";
$outputDir.=$parahash{output_folder};
if (! -e $outputDir) { system("mkdir $outputDir"); }
system(qq(cp $ARGV[0] $outputDir));
$parahash{output_folder}=$outputDir;

# deal with different mode
if ($parahash{'mode'} eq 'mutation')
{
	run_mutation(\%parahash,$code_path);
}
elsif ($parahash{'mode'} eq 'junction')
{
	run_junction(\%parahash,$code_path);
}
elsif ($parahash{'mode'} eq '6FT')
{
	run_6FT(\%parahash,$code_path);
} elsif ($parahash{'mode'} eq 'transcript_6FT') {
	run_transcript_6FT(\%parahash,$code_path);
} else { die "Unexpected running mode for jump -g: $parahash{'mode'}!!!\n"; }
#print "Done.\n";
#-------------------------------------------------------------------------------------
sub parse_params
{
        my ($par,$parahash)=@_;
        open(IN,$par);
        my $k=0;
        while (<IN>)
        {
                s/^\s+//; # rm space at front
                next if (/^#/); # skip comment lines
                chomp;

                if (/ = /) # normal parameters
                {
                        s/\s*([;\#].*)$//; # delete comments
                        my ($key,$value) = split(' = ',$_);
                        $$parahash{$key}=$value;
                }
 	}
        close IN;
}

sub isSangerPhred
{
	my ($fastq,$n)=@_;

	#$n |= 10000;
	my ($k,$b,$s); my $line;
	$k=$b=$s=0; 
	open(IN,$fastq) || die "cannot open fastq file:$fastq\n";
	print "Checking Phred score type using the first $n rows of fastq file ($fastq) ...\n";
	while(<IN>)
	{
		$k++;
		$line=<IN>; $k++;
		$line=<IN>; $k++;
		$line=<IN>; $k++;
		chomp($line);
		if ($line =~ m/BBBBBBBBBB$/) { $b++; }
		elsif ($line =~ m/\#\#\#\#\#\#\#\#\#\#$/) { $s++; }
		last if ($k>=$n);
	}
	close IN;
	print "  Sanger label row counts: $s\n";
	print "  Solexa label row counts: $b\n";

	return ($s>$b)?1:0;
}

sub run_transcript_6FT {
	my ($parahash,$code_path)=@_;
	my $rdt_path="$code_path/transcript_6FT";
	print "Creating transcript 6FT database:\n";
	my $input=$$parahash{input_file};
	my $outputDir=$$parahash{output_folder};

	my @t=split /\//,$input;
	my $input_base=$t[$#t];

	my $output="$outputDir/$input_base";
	$output =~ s/\.[^.]+$//;
	#my @out; $out[1]=$output;

	system(qq(perl $rdt_path/trinity_tr6.pl $input $$parahash{min_ORF_AA} > $output\_peptides.fas));
	system(qq(perl $code_path/shortID_convert.pl $output\_peptides.fas $output\_peptides.newIDs));
	system(qq(perl $rdt_path/formatFas.pl $input > $output\_formated.fasta));
}

sub run_6FT
{
	my ($parahash,$code_path)=@_;
	my $rdt_path="$code_path/RDT";

	#print "Running 6FT mode:\n";
	print "Creating 6FT database:\n";

	my $input=$$parahash{input_file};
	my $outputDir=$$parahash{output_folder};
	my $input_base;
	if (1)
	{
		my @t=split /\//,$input;
		$input_base=$t[$#t];
	}
	my $output="$outputDir/$input_base";
	$output =~ s/\.[^.]+$//;;
	my @out; $out[1]=$output;

	#print "Fastq reads QC ...\n";
	my $sangerPhred=isSangerPhred($input,10000);
	if ($sangerPhred)
	{
		print "Phred score type: sanger\n";
		print "Fastq reads QC ...\n";
		system(qq(perl $rdt_path/readQC_IDshort.pl $input sanger $out[1]_QC.fas));
	}
	else
	{
		print "Phred score type: solexa\n";
		print "Fastq reads QC ...\n";
		system(qq(perl $rdt_path/readQC_IDshort.pl $input solexa $out[1]_QC.fas));
	}
	print "Translating reads to amino acid sequences ...\n";
	system(qq(perl $rdt_path/fas_tr6.pl $out[1]_QC.fas $out[1]_QC));
	print "Performing in silico digestion of amino acid sequences ...\n";
	system(qq(perl $rdt_path/tryptic_digest.pl $out[1]_QC_peptide.fas KR P $out[1]_QC_peptide_digested.fas));
	print "Removing duplicated peptides ...\n";
	system(qq(perl $rdt_path/uniq_digested_peptides.pl $out[1]_QC_peptide_digested.fas $out[1]_QC_peptide_digested_uniq  >/dev/null 2>&1));
	system(qq(perl $rdt_path/pepToReads.pl $out[1]_QC_peptide_digested_uniq.readCounts $out[1]_QC_peptide_digested_uniq.pepToReads));
	if (!defined($$parahash{read_coverage_cutoff})) 
	{
		$$parahash{read_coverage_cutoff}=1;
	}
	if ($$parahash{read_coverage_cutoff}>1)
	{
		print "Removing peptides with read coverage less than $$parahash{read_coverage_cutoff} ...\n";
		system(qq(perl $rdt_path/multiSptPep.pl $out[1]_QC_peptide_digested_uniq.readCounts $$parahash{read_coverage_cutoff} $out[1]_QC_peptide_digested_uniq_$$parahash{read_coverage_cutoff}reads.fas));
	}
	else
	{
		system(qq(mv $out[1]_QC_peptide_digested_uniq.fas $out[1]_QC_peptide_digested_uniq_$$parahash{read_coverage_cutoff}reads.fas));
	}
	system(qq(perl $code_path/shortID_convert.pl $out[1]_QC_peptide_digested_uniq_$$parahash{read_coverage_cutoff}reads.fas $out[1]_QC_peptide_digested_uniq_$$parahash{read_coverage_cutoff}reads.newIDs));

	# .jump_c_tmp
	open(OUT,">>",'.jump_c_tmp');
	print OUT "$out[1]_QC_peptide_digested_uniq_$$parahash{read_coverage_cutoff}reads.newIDs.fas\n";
	close OUT;
}

sub run_junction
{
	my ($parahash,$code_path)=@_;
	my $jun_path="$code_path/MFM_junction";

	#print "Running junction mode:\n";
	print "Creating junction peptide database:\n";

	my ($genome,$geneStrand);
	if ( !defined($$parahash{reference_genome}) ) { 
		die "Please spcify reference_genome!!!\n"; 
	} else {
		$genome=$$parahash{reference_genome};
		if (! -e $genome) {
			die "reference_genome file does not exist: $genome!!!\n";
		}
	}

	my $input=$$parahash{input_file};
	my $outputDir=$$parahash{output_folder};
	my $input_base;
	if (1)
	{
		my @t=split /\//,$input;
		$input_base=$t[$#t];
	}
	my $output="$outputDir/$input_base";
	#$output =~ s/\.txt$//;
	#$output =~ s/\.[^.]+$//;

	system(qq(perl $jun_path/STARjunc2tab.pl $input > $output.junc));
	$output = "$output.junc";
	print "Extracting genome information ...\n";
	#system(qq(perl $jun_path/extractFlankingGenomicSeq.pl $genome $input 66 $output.fas));
	#system(qq(perl $jun_path/extractFlankingGenomicSeq.pl $genome $input 90 $output.fas));
	system(qq(perl $jun_path/extractFlankingGenomicSeq.pl $genome $output 90 $output.fas));
	print "Translating nucleotide sequences to amino acid sequences ...\n";
	#system(qq(perl $jun_path/junction_seq_translation.pl $output.fas $geneStrand $output.AA.fas));
	system(qq(perl $jun_path/junction_seq_translation.pl $output.fas $output.AA.fas));
	system(qq(perl $jun_path/check_AA_sequence.pl $output.AA.fas $output.AA.central.fas));
	system(qq(perl $code_path/shortID_convert.pl $output.AA.central.fas $output.AA.central.newIDs));

	# .jump_c_tmp
	open(OUT,">>",'.jump_c_tmp');
	print OUT "$output.AA.central.newIDs.fas\n";
	close OUT;
}

sub run_mutation
{
	my ($parahash,$code_path)=@_;
	my $annovar_path="$code_path/annovar";
	my $mut_path="$code_path/MFM_SNV";

	#print "Running mutation mode:\n";
	print "Creating mutation peptides database:\n";

	my ($annovarLib,$annovarAnn,$anVersion);
	my $annovarCode="$annovar_path/annotate_variation.pl";
	#if ( !defined($$parahash{species}) ) { die "Please spcify species!!!\n"; }
	if ( !defined($$parahash{annovar_annotation_file}) ) { 
		die "Please spcify annovar_annotation_file!!!\n"; 
	}else {
		if (! -e $$parahash{annovar_annotation_file}) {
			die "file does not exist: $$parahash{annovar_annotation_file}!!!\n";
		}

		my @t=split /\//,$$parahash{annovar_annotation_file};
		my $filename=pop @t;
		$annovarLib=join("/",@t);
		$filename =~ s/\.txt$//;
		($anVersion,$annovarAnn)=split /\_/,$filename
	}

=head
	if ( $$parahash{species} eq 'human' ) 
	{
		$annovarLib="$annovar_path/humandb/";
		$annovarAnn='knowngene';
	}
	else { die "Unexpected species ($$parahash{species})!!!\nCurrently only support human.\n"; }
=cut

	my $input=$$parahash{input_file};
	my $outputDir=$$parahash{output_folder};
	my $input_base;
	if (1)
	{
		my @t=split /\//,$input;
		$input_base=$t[$#t];
	}
	my $output="$outputDir/$input_base"; 
	#print "$output\n";
	#system(qq(cd $outputDir));
	print "Running mutation annotation (AnnoVar) ...\n";
	system(qq(perl $annovarCode -buildver $anVersion -geneanno --seq_padding 30 -dbtype $annovarAnn $input $annovarLib  >/dev/null 2>&1));
	system(qq(mv $input.* $outputDir));
	print "Generating amino acid sequences ...\n";
	system(qq(perl $mut_path/seqpad2fas.pl $output.seqpad mut $output.seqpad.mut.fas));
	system(qq(perl $mut_path/seqpad2fas.pl $output.seqpad ref $output.seqpad.ref.fas));
	system(qq(perl $mut_path/rm_dup.pl $output.seqpad.mut.fas $output.seqpad.mut.uniq.fas));
	system(qq(perl $code_path/shortID_convert.pl $output.seqpad.mut.uniq.fas $output.seqpad.mut.uniq.newIDs));

	# .jump_c_tmp
	open(OUT,">>",'.jump_c_tmp');
	print OUT "$output.seqpad.mut.uniq.newIDs.fas\n";
	close OUT;
}
