#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/g

use strict;
use warnings;
use Cwd;

use PrimarySeq;
use IDtxt_parser;

if (scalar(@ARGV)!=1)
{
        die "Usage: perl jump_g_postAnnotations.pl jump_g_pa.params\n";
}

# code path
my $code_path="/home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/pa";
#my $annotationFile2='/home/yli4/annotations/knownGenes_uniPro_012314.txt';

#initialization
my (%parahash);
parse_params($ARGV[0],\%parahash);
my $annotationFile2=$parahash{gene_ID_convert_file};

# output folder
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
	#print "Filtering RNA-seq 6FT ID.txt ($parahash{ID_txt}) by reference proteins ($parahash{ref_proteins})\n";
	#filterIDtxtByRef($parahash{ID_txt},$parahash{ref_proteins},"$outputDir/filtered_ID.txt","$outputDir/ref_protein_ID.txt");
	#$parahash{ID_txt}="$outputDir/filtered_ID.txt";
        run_6FT(\%parahash,$code_path);
#} elsif ($parahash{'mode'} eq 'ref') {
} elsif ($parahash{'mode'} eq 'ref_proteins') {
	run_ref(\%parahash,$code_path);
} elsif ($parahash{'mode'} eq 'transcript_6FT') {
	run_transcript_6FT(\%parahash,$code_path);
} else { 
	die "Unexpected running mode for jump -g: $parahash{'mode'}!!!\n"; 
}

#-------------------------------------------------------------------------------------
sub run_transcript_6FT {
	my ($parahash,$code_path)=@_;
	my $local_path="$code_path/transcript_6FT";
	my $output='transcript_6FT';

	my $blat="$local_path/blat";
	my $trackName=$$parahash{track_name};
	my $IDtxt=$$parahash{ID_txt};
	my $refFlat=$$parahash{annovar_annotation_file};
	$$parahash{RNA_transcript}=$$parahash{core_file_name};

	print "Annotating 6FT peptides ...\n";
	#system(qq(export PERL5LIB=''));
	# annoPep.pl
	if (defined($$parahash{ref_proteins})) {
		system(qq(export PERL5LIB='' && cd $$parahash{output_folder} && perl $local_path/annoPep.pl -idtxt $IDtxt -rna $$parahash{RNA_transcript} -path-to-blat $blat -reference-genome $$parahash{reference_genome} -reference-protein $$parahash{ref_proteins} -refFlat $refFlat));
	} else {
		system(qq(cd $$parahash{output_folder} && perl $local_path/annoPep.pl -idtxt $IDtxt -rna $$parahash{RNA_transcript} -path-to-blat $blat -reference-genome $$parahash{reference_genome} -refFlat $refFlat));
	}
	# BED file
	system(qq(cd $$parahash{output_folder} && perl $code_path/scanCounts_BED.pl output.bed > output.sc.bed ));
}

sub run_ref
{
	my ($parahash,$code_path)=@_;
	my $local_path="$code_path/ref";
	my $output='reference';
	my $trackName=$$parahash{track_name};
	my $IDtxt=$$parahash{ID_txt};

	#my $ucscProtein='/home/yli4/genomes/protein_sequence_hg19_06272013.fas';
	#my $refFlat='hg19_knownGene.txt';
	my $refFlat=$$parahash{annovar_annotation_file};;
	#my $uniPro2ucsc='uniProt2015toUCSChg19_v1.0.txt';
	#my $uniPro2ucsc=$$parahash{gene_ID_convert_file};
	my $uniPro2ucsc=$annotationFile2;

	#print "Running ref mode:\n";
	print "Annotating reference peptides ...\n";
	system(qq(cd $$parahash{output_folder} && perl $local_path/uniPro2ucsc.pl $IDtxt $uniPro2ucsc > $IDtxt.ucsc.txt));
	#system(qq(cd $$parahash{output_folder} && perl $local_path/IDtxt2BED_ucsc.pl $IDtxt.ucsc.txt $refFlat $ucscProtein $trackName > $output\_peptides.bed));
	system(qq(cd $$parahash{output_folder} && perl $local_path/IDtxt2BED_ucsc.pl $IDtxt.ucsc.txt $refFlat $trackName  $output\_peptides.bed >/dev/null 2>&1));
	system(qq(cd $$parahash{output_folder} && perl $code_path/scanCounts_BED.pl $output\_peptides.bed > $output\_peptides.sc.bed));

}

sub filterIDtxtByRef
{
	my ($IDtxt,$refPro,$output,$output2)=@_;

	# parse IDtxt
	my $idp=IDtxt_parser->new();
	my ($scanhash,$peptidehash)=$idp->parse_IDtxt($IDtxt);

	# parse ref prot FASTA
	my $ps=PrimarySeq->new();
	my $seqhash=$ps->parseFas($refPro);

	# for each scan in %scanhash: mark if ref peptides
	foreach my $outfile (keys %$scanhash)
	{
		my $nomod=$idp->nomod_pep($scanhash->{$outfile}->{peptide});
		foreach my $pro (keys %$seqhash)
		{
			if ($seqhash->{$pro} =~ m/$nomod/) {
				$scanhash->{$outfile}->{ref_protein}=1;
				last;
			}
		}
	}

	# only print out non-ref peptide lines
	open(OUT,">$output");
	open(OUT2,">$output2");
	foreach my $outfile (keys %$scanhash)
	{
		# ref peptides
		if (defined($scanhash->{$outfile}->{ref_protein}) 
		and $scanhash->{$outfile}->{ref_protein} == 1) {
			print OUT2 $scanhash->{$outfile}->{IDtxt_line},"\n";
		# novel peptides
		} else {
			print OUT $scanhash->{$outfile}->{IDtxt_line},"\n";
		}
	}
	close OUT;
	close OUT2;
}

sub run_junction
{
	my ($parahash,$code_path)=@_;
	my $local_path="$code_path/junction";
	my $core_file_name="$$parahash{core_file_name}\.junc\.AA";
	my $output='junction';
	my $trackName=$$parahash{track_name};
	my $IDtxt=$$parahash{ID_txt};
	my $junc_AApadding=30;

	my $refFlat='/home/yli4/downloads/svinframe/resources/refFlat.15AUG2014.txt';
	#my $refFlat='/home/yli4/downloads/svinframe/resources/knownGene_refFlat.txt';
	my $refGenome='/home/yli4/genomes/GRCh37-lite.fa';
	my $refgene2protein='/home/yli4/downloads/svinframe/resources/refgene2protein.tab';

	#print "Running junction mode:\n";
	print "Annotating junction peptides ...\n";
	system(qq(cd $$parahash{output_folder} && perl $local_path/consolidate_junction_peptides.pl $IDtxt $core_file_name.fas $junc_AApadding $output));
	system(qq(cd $$parahash{output_folder} && perl $local_path/junction_frameCheck_input.pl $output\_peptides_frame.txt $core_file_name.fas $junc_AApadding $output\_frameCheck_input.tab));
	system(qq(cd $$parahash{output_folder} && perl $local_path/ID2bed.pl valideted_ID.txt $core_file_name.fas $junc_AApadding $trackName $output\_peptides.bed));
	system(qq(cd $$parahash{output_folder} && perl $code_path/scanCounts_BED.pl $output\_peptides.bed > $output\_peptides.sc.bed));
	# only for internal use currently
	#system(qq(cd $$parahash{output_folder} && perl $local_path/sv_inframe.pl -refflat $refFlat -refgene2protein $refgene2protein -fasta $refGenome -single $output\_frameCheck_input.tab -junction-mode >/dev/null 2>&1));
	# only for internal use currently
	#system(qq(cd $$parahash{output_folder} && perl $local_path/classify_junc.pl $refFlat $output\_frameCheck_input.tab.frame.tab junction_types >/dev/null 2>&1));
}

sub run_6FT
{
	my ($parahash,$code_path)=@_;
	my $local_path="$code_path/6FT";
	my $core_file_name="$$parahash{core_file_name}\_QC";
	my $output='RNAseq_6FT';
	my $trackName=$$parahash{track_name};
	my $IDtxt=$$parahash{ID_txt};

	$$parahash{reference_CDS} =~ s/\.[a-zA-Z]+$//;
	$$parahash{reference_genome} =~ s/\.[a-zA-Z]+$//;
	$$parahash{reference_mRNA} =~ s/\.[a-zA-Z]+$//;
	$$parahash{reference_gene_locus} =~ s/\.[a-zA-Z]+$//;

	print "Annotating RNAseq 6FT peptides ...\n";
	print "Retrieving RNA-seq reads of accpeted peptides ...\n";
	system(qq(cd $$parahash{output_folder} && perl $local_path/extract_MS_reads.pl $IDtxt $core_file_name\_peptide_digested_uniq.seqIDs $core_file_name.fas $output));

	print "Running alignment ...\n";
	#system(qq(cd $$parahash{output_folder} && sh $local_path/megablast_run.sh $output\_reads.fas $output >/dev/null 2>&1));
	#system(qq(cd $$parahash{output_folder} && megablast -d $$parahash{reference_CDS} -i $output\_reads.fas  -e 1 -D 3 -v 3 > $output\_reads_vs_CDS.tab 2>&1));
	#system(qq(cd $$parahash{output_folder} && megablast -d $$parahash{reference_mRNA} -i $output\_reads.fas  -e 1 -D 3 -v 3 > $output\_reads_vs_mRNA.tab 2>&1));
	#system(qq(cd $$parahash{output_folder} && megablast -d $$parahash{reference_genome} -i $output\_reads.fas  -e 1 -D 3 -v 3 > $output\_reads_vs_WG.tab 2>&1));
	#system(qq(cd $$parahash{output_folder} && megablast -d $$parahash{reference_gene_locus} -i $output\_reads.fas  -e 1 -D 3 -v 3 > $output\_reads_vs_geneLocus.tab 2>&1));
	system(qq(cd $$parahash{output_folder} && $local_path/blastall -p blastn -d $$parahash{reference_CDS} -i $output\_reads.fas -m 8  -e 1 -o $output\_reads_vs_CDS.tab 2>&1));
	system(qq(cd $$parahash{output_folder} && $local_path/blastall -p blastn -d $$parahash{reference_mRNA} -i $output\_reads.fas  -m 8  -e 1 -o $output\_reads_vs_mRNA.tab 2>&1));
	system(qq(cd $$parahash{output_folder} && $local_path/blastall -p blastn -d $$parahash{reference_genome} -i $output\_reads.fas  -m 8  -e 1 -o $output\_reads_vs_WG.tab 2>&1));
	system(qq(cd $$parahash{output_folder} && $local_path/blastall -p blastn -d $$parahash{reference_gene_locus} -i $output\_reads.fas  -m 8  -e 1 -o $output\_reads_vs_geneLocus.tab 2>&1));

	print "Parsing alignment ...\n";
	system(qq(cd $$parahash{output_folder} && perl $local_path/build_out2aligns.pl $output\_reads.fas $output\_out2read.hash $output\_pep2readPos.txt $output\_reads.alg $output\_reads_vs_CDS.tab $output\_reads_vs_mRNA.tab $output\_reads_vs_contaminants.tab $output\_reads_vs_altExon.tab $output\_reads_vs_WG.tab >/dev/null 2>&1));
	system(qq(cd $$parahash{output_folder} && perl $local_path/IDtxt2BED_rdt_ucscVersion.pl $output\_reads.alg $trackName $output\_peptides.bed));
	system(qq(cd $$parahash{output_folder} && perl $local_path/scanCounts_BED.pl $output\_peptides.bed > $output\_peptides.sc.bed));
	system(qq(cd $$parahash{output_folder} && perl $local_path/annotateIntronic.pl $output\_reads.alg.updated $output\_reads_vs_geneLocus.tab $output\_pep2readPos.txt > $output\_reads.alg.updated2));
	system(qq(cd $$parahash{output_folder} && perl $local_path/pepTab.pl $output\_reads.alg.updated2 > $output.unique.peptides.txt));
	system(qq(cd $$parahash{output_folder} && perl $local_path/publication_table.pl $annotationFile2 $output\_peptides.sc.bed $output.unique.peptides.txt > $output\_peptides.txt));
	validated_IDtxt("$$parahash{output_folder}/$output\_peptides.txt",$IDtxt,"$$parahash{output_folder}/valideted_ID.txt");

	# further annotate peptides:
	# attach genomic sequence
	# highlight single AA change
	system(qq(cd $$parahash{output_folder} && perl $local_path/attachGenomicPos.pl $output\_peptides.sc.bed $output.unique.peptides.txt > $output.unique.peptides.pos.txt));
	system(qq(cd $$parahash{output_folder} && perl $local_path/attachGenomicSeq.pl ~/genomes/hg19_genome_raw.fa $output.unique.peptides.pos.txt $output.unique.peptides.pos.nt.txt));
	system(qq(cd $$parahash{output_folder} && perl $local_path/attachRefAA.pl $output.unique.peptides.pos.nt.txt > $output.unique.peptides.pos.nt.aa.txt));
}

sub validated_IDtxt
{
	my ($pepTab,$IDtxt,$output)=@_;

	# parse pepTab => %pephash{$pep}
	my %pephash;
	open(IN,$pepTab) || die "Cannot open peptide table $pepTab\n";
	while (<IN>){
		next if (/^Peptide/ || /^peptide/); # header
		chomp;
		my $pep=(split /\t/,$_)[0];
		$pephash{$pep}='';
	}
	close IN;

	# filter IDtxt and print
	my $idp=IDtxt_parser->new();
	open(IN,$IDtxt) || die "Cannot open IDtxt $IDtxt\n";
	open(OUT,">",$output);
	while (<IN>){
		next if (/^Database/ || /^Peptide/); # header
		chomp;
		my $pep=(split /\;/,$_)[0];
		$pep=$idp->intpep($pep);
		if (defined($pephash{$pep})) {
			print OUT $_,"\n";
		}
	}
	close IN;
	close OUT;
}

sub run_mutation
{
	my ($parahash,$code_path)=@_;
	my $mut_path="$code_path/mutation";
	my $core_file_name=$$parahash{core_file_name};
	my $IDtxt=$$parahash{ID_txt};
	#my $output=$$parahash{output_folder};
	my $output='mutation';
	my $trackName=$$parahash{track_name};

	#print "Running mutation mode:\n";
	print "Annotating mutation peptides ...\n";

	# mutation annotation
	system(qq(cd $$parahash{output_folder} && perl $mut_path/consolidate_IDtxt2SNV.pl $core_file_name.seqpad.ref.fas $core_file_name.seqpad.mut.fas $IDtxt $output));
	# visulization
	system(qq(cd $$parahash{output_folder} && perl $mut_path/acpMut_IDtxt.pl $output\_pep.txt ID_posCheck.txt > valideted_ID.txt));
	system(qq(cd $$parahash{output_folder} && perl $mut_path/IDtxt2BED_mut_ucscVersion.pl valideted_ID.txt $trackName $output\_peptides.bed));
	system(qq(cd $$parahash{output_folder} && perl $mut_path/scanCounts_BED.pl $output\_peptides.bed > $output\_peptides.sc.bed));
	system(qq(cd $$parahash{output_folder} && perl $mut_path/pgSnp.pl $output\_mut_all.txt > $output\_mut_all.pgSnp));
}

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


