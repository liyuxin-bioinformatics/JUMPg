#!/usr/bin/perl -I /home/yli4/development/JUMPg/JUMPg_v2.3.5/programs/g

## Created: 08/14/2015
my $VERSION='2.2.7';

use strict;
use warnings;
use Cwd;
use File::Basename;
use IDtxt_parser;
use pepXML_parser;

# usage
if(scalar(@ARGV)<2)
{
	help();
}
my $params=shift @ARGV;
my $currentdir=getcwd;

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump proteogenomics                     ****     #
#       ****  Version $VERSION                           ****     #
#       ****  Copyright (C) 2014 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################

EOF


# check parameter version
#my $VERSION = '1.0.0';
#check_parameter_version($params,$VERSION);

# initialization
my (%parahash);
parse_params($params,\%parahash);
set_program_paths(\%parahash);

# set up output folders
my $customizedDB_dir=$currentdir;
$customizedDB_dir.='/customizedDB';
if (!-e "$customizedDB_dir") {system(qq(mkdir $customizedDB_dir));}
my $output_dir=$currentdir;
$output_dir.="/gnm_$parahash{output_directory}";
if (!-e "$output_dir") {system(qq(mkdir $output_dir));}
system(qq(cp $params $output_dir)); # cp parameter file
$output_dir="$output_dir/intermediate"; # all true runs in intermediate/
if (!-e "$output_dir") {system(qq(mkdir $output_dir));}

#=head # for debug
# jump -params
#if (!defined($parahash{bypass_params}) or $parahash{bypass_params}==0)
if (1)
{
	system(qq(cd $output_dir && perl $parahash{'programs'}{params}));
	set_params_files(\%parahash);
}

print <<EOF;
-------------------------------------------------------------------------------------
Building customized peptide database
-------------------------------------------------------------------------------------
EOF

# jump -c
# initialize .jump_c_tmp
if (-e "$customizedDB_dir/.jump_c_tmp") { system(qq(rm $customizedDB_dir/.jump_c_tmp)); }
system(qq(touch $customizedDB_dir/.jump_c_tmp));
# run jump -c
if (1)
{
	my $mode='mutation';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $fasfile=get_mutation_fas_name(\%parahash,$currentdir);
		#if (-e get_mutation_fas_name(\%parahash,$currentdir))
		if (-e $fasfile and -s $fasfile)
		{
			print "Customized database already generated: ",get_mutation_fas_name(\%parahash,$currentdir),"\n";
			print "Bypass $mode database generation\n";
			my $line=get_mutation_fas_name(\%parahash,$currentdir);
			$line="$line\n";
			print2file("$customizedDB_dir/.jump_c_tmp",$line);
		}
		else
		{
			my $params_c=params_c_update(\%parahash,$mode,$currentdir);
			system(qq(cd customizedDB && perl $parahash{'programs'}{c} $params_c));
			unless (-e $fasfile and -s $fasfile) {
				die "Fail to genenrate $mode peptide database!!!\nPlease check your parameters and annotation files.\n";
			}
		}
	}
	$mode='junction';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $fasfile=get_junction_fas_name(\%parahash,$currentdir);
		#if (-e get_junction_fas_name(\%parahash,$currentdir))
		if (-e $fasfile and -s $fasfile)
		{
			print "Customized database already generated: ",get_junction_fas_name(\%parahash,$currentdir),"\n";
			print "Bypass $mode database generation\n";
			my $line=get_junction_fas_name(\%parahash,$currentdir);
			$line="$line\n";
			print2file("$customizedDB_dir/.jump_c_tmp",$line);
		}
		else
		{
			my $params_c=params_c_update(\%parahash,$mode,$currentdir);
			system(qq(cd customizedDB && perl $parahash{'programs'}{c} $params_c));
			unless (-e $fasfile and -s $fasfile) {
				die "Fail to genenrate $mode peptide database!!!\nPlease check your parameters and annotation files.\n";
			}
		}
	}

	$mode='transcript_6FT';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0') {
		my $fasfile=get_transcript_6FT_fas_name(\%parahash,$currentdir);
		if (-e $fasfile and -s $fasfile) {
			print "Customized database already generated: $fasfile\n";
			print "Bypass $mode database generation\n";
			#my $line=$fasfile;
			#$line="$line\n";
		} else {
			my $params_c=params_c_update(\%parahash,$mode,$currentdir);
			system(qq(cd customizedDB && perl $parahash{'programs'}{c} $params_c));
			unless (-e $fasfile and -s $fasfile) {
				die "Fail to genenrate $mode peptide database!!!\n($fasfile is expected but failed to pass check out)\nPlease check your parameters and annotation files.\n";
			}
		}
		my $line="$fasfile\n";
		print2file("$customizedDB_dir/.jump_c_tmp",$line);
	}

	$mode='6FT';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $fasfile=get_6FT_fas_name(\%parahash,$currentdir);
		#if (-e get_6FT_fas_name(\%parahash,$currentdir))
		if (-e $fasfile and -s $fasfile)
		{
			print "Customized database already generated: ",get_6FT_fas_name(\%parahash,$currentdir),"\n";
			print "Bypass $mode database generation\n";
			my $line=get_6FT_fas_name(\%parahash,$currentdir);
			$line="$line\n";
			print2file("$customizedDB_dir/.jump_c_tmp",$line);
		}
		else
		{
			my $params_c=params_c_update(\%parahash,$mode,$currentdir);
			system(qq(cd customizedDB && perl $parahash{'programs'}{c} $params_c));
			unless (-e $fasfile and -s $fasfile) {
				die "Fail to genenrate $mode peptide database!!!\nPlease check your parameters and annotation files.\n";
			}
		}
	}
}

# jump -d
# update -s parameter first
params_update_by_hash("$output_dir/ParameterFiles/$parahash{MS_data_types}/$parahash{'params'}{s}","$output_dir/$parahash{'params'}{s}",\%parahash);
# jump -d
$parahash{'params'}{d}=params_d_update(\%parahash,$output_dir,$currentdir);
system(qq(cd $output_dir/database && perl $parahash{'programs'}{d} $parahash{'params'}{d}));

print <<EOF;

-------------------------------------------------------------------------------------
Running database search
-------------------------------------------------------------------------------------
EOF

# jump -s
#if (0){ # for debug:  to skipp search
if (!defined($parahash{bypass_database_search}) or $parahash{bypass_database_search}==0)
{
	params_s_update("$output_dir/database/.jump_d_tmp","$output_dir/$parahash{'params'}{s}");
	add_fast_run_label("$output_dir/$parahash{'params'}{s}");

	if (inputRawFiles($ARGV[0]))
	{
		foreach my $rawFile (@ARGV) { 
			system(qq(ln -s $currentdir/$rawFile $output_dir)); 
		}
		system(qq(cd $output_dir && $parahash{'programs'}{s} -p $parahash{'params'}{s} @ARGV));
	} else {
		system(qq(cp $ARGV[0] $output_dir));
		system(qq(cd $output_dir && $parahash{'programs'}{rundtas} $parahash{'params'}{s} $ARGV[0]));
	}
} else {
	print "Bypass database search\n";
} # -s
#} # for debug:  to skipp search

print <<EOF;

-------------------------------------------------------------------------------------
Performing PSM filtering
-------------------------------------------------------------------------------------
EOF

# jump -f
if (!defined($parahash{bypass_PSM_filter}) or $parahash{bypass_PSM_filter}==0)
{
	params_update_by_hash("$output_dir/ParameterFiles/$parahash{MS_data_types}/$parahash{'params'}{f}","$output_dir/$parahash{'params'}{f}",\%parahash);
	# run for accepted PSMs
	# set input path
	params_f_update("$output_dir/.jump_s_tmp","$output_dir/$parahash{'params'}{f}",'accepted_PSMs',"$output_dir/$parahash{'params'}{f}");
	system(qq(cd $output_dir && $parahash{'programs'}{f} $parahash{'params'}{f}));
	# run for confident PSMs (FDR = 0)
	#if (1)
	if ($parahash{ref_proteins} ne '0' and inputRawFiles($ARGV[0]) and $parahash{digestion} eq 'full')
	{
		$parahash{'unique_protein_or_peptide'} = 'peptide';
		$parahash{'FDR'} = 0;
		params_update_by_hash("$output_dir/$parahash{'params'}{f}","$output_dir/$parahash{'params'}{f}",\%parahash);
		params_f_update("$output_dir/.jump_s_tmp","$output_dir/$parahash{'params'}{f}",'confident_PSMs',"$output_dir/$parahash{'params'}{f}");
		system(qq(cd $output_dir && $parahash{'programs'}{f} $parahash{'params'}{f}  >/dev/null 2>&1));
		if (-e "$output_dir/.sum_confident_PSMs") { system(qq(rm -rf $output_dir/.sum_confident_PSMs)); }
		system(qq(mv $output_dir/sum_confident_PSMs $output_dir/.sum_confident_PSMs));
	}
} else {
	print "Bypass PSM filtering\n";
} # -f

print <<EOF;

-------------------------------------------------------------------------------------
Filtering spectra by MS2 quality score
-------------------------------------------------------------------------------------
EOF

# jump -qc
if (!defined($parahash{bypass_quality_filter}) or $parahash{bypass_quality_filter}==0)
{
	my $conf_f;
	if ($parahash{ref_proteins} ne '0' and inputRawFiles($ARGV[0]) and $parahash{digestion} eq 'full') {
		$conf_f="$output_dir/.sum_confident_PSMs/ID.txt";
	} else { 
		$conf_f=0; 
	}

	my $acpt_f="$output_dir/sum_accepted_PSMs/ID.txt";
	params_qc_update("$output_dir/.jump_s_tmp",$conf_f,$acpt_f,$parahash{'PSM_recoveray_rate'},'accepted_PSMs',"$output_dir/jump_qc.params");
	system(qq(cd $output_dir && perl $parahash{'programs'}{qc} jump_qc.params));
	system(qq(cp $output_dir/.jump_qc_tmp qc_MSMS_input.txt));
} else {
	print "Bypass MS2 quality scoring\n";
}  # -qc
#=cut # for debug

# summary
if (1)
{
	my ($mut_ids, $jun_ids, $rdt_ids)=set_summary_params(\%parahash,$currentdir);
	system(qq(perl $parahash{'programs'}{summary} $output_dir/sum_accepted_PSMs $mut_ids $jun_ids $rdt_ids $output_dir/JUMPg_results));
}

print <<EOF;

-------------------------------------------------------------------------------------
Annotating identified peptides
-------------------------------------------------------------------------------------
EOF

# postAnnotations
if (!defined($parahash{bypass_peptide_annotation}) or $parahash{bypass_peptide_annotation}==0)
{
	my $mode='mutation';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $bname=basename($parahash{$mode});
		my $params_pa=params_pa_update(\%parahash,$mode,$output_dir,$currentdir,$bname);
		system(qq(cd $output_dir/JUMPg_results && perl $parahash{'programs'}{pa} $params_pa));
	}
	$mode='junction';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $bname=basename($parahash{$mode});
		#$bname =~ s/\.[^.]+$//; 
		#print "junction bname: $bname\n";
		my $params_pa=params_pa_update(\%parahash,$mode,$output_dir,$currentdir,$bname);
		system(qq(cd $output_dir/JUMPg_results && perl $parahash{'programs'}{pa} $params_pa));
	}

	$mode='transcript_6FT';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0') {
		my $bname=basename($parahash{$mode});
		$bname =~ s/\.[^.]+$//;
		my $params_pa=params_pa_update(\%parahash,$mode,$output_dir,$currentdir,$bname);
		system(qq(cd $output_dir/JUMPg_results && perl $parahash{'programs'}{pa} $params_pa));
	}

	$mode='6FT';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $bname=basename($parahash{$mode});
		#print "bname: $bname\n";
		$bname =~ s/\.[^.]+$//;
		#print "bname: $bname\n";
		my $params_pa=params_pa_update(\%parahash,$mode,$output_dir,$currentdir,$bname);
		system(qq(cd $output_dir/JUMPg_results && perl $parahash{'programs'}{pa} $params_pa));
	}
	$mode='ref_proteins';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		my $bname='';
		my $params_pa=params_pa_update(\%parahash,$mode,$output_dir,$currentdir,$bname);
		system(qq(cd $output_dir/JUMPg_results && perl $parahash{'programs'}{pa} $params_pa));
	}
}

print <<EOF;

-------------------------------------------------------------------------------------
Printing publication tables
-------------------------------------------------------------------------------------
EOF

# publication DIR
$output_dir =~ s/intermediate$//;
if (!-e "$output_dir/publications") {system(qq(mkdir $output_dir/publications));}
if (1)
{
	print "Summary of JUMPg peptide identifications:\n";
	my $mode='mutation';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		system(qq(cp $output_dir/intermediate/JUMPg_results/mutation/mutation_peptides.txt $output_dir/publications));
		system(qq(cp $output_dir/intermediate/JUMPg_results/mutation/mutation_events.txt $output_dir/publications));
		system(qq(cp $output_dir/intermediate/JUMPg_results/mutation/mutation_peptides.sc.bed $output_dir/publications/mutation_peptides.bed));
	}
	$mode='junction';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		system(qq(cp $output_dir/intermediate/JUMPg_results/junction/junction_peptides.txt $output_dir/publications));
		system(qq(cp $output_dir/intermediate/JUMPg_results/junction/junction_events.txt $output_dir/publications));
		system(qq(cp $output_dir/intermediate/JUMPg_results/junction/junction_peptides.sc.bed $output_dir/publications/junction_peptides.bed));
	}
	$mode='transcript_6FT';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0') {
		system(qq(cp $output_dir/intermediate/JUMPg_results/transcript_6FT/id_peptide.txt $output_dir/publications/transcript_6FT_peptides.txt));
		system(qq(cp $output_dir/intermediate/JUMPg_results/transcript_6FT/output.sc.bed $output_dir/publications/transcript_6FT_peptides.bed));
	}
	$mode='6FT';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		system(qq(cp $output_dir/intermediate/JUMPg_results/6FT/RNAseq_6FT_peptides.txt $output_dir/publications));
		system(qq(cp $output_dir/intermediate/JUMPg_results/6FT/RNAseq_6FT_peptides.sc.bed $output_dir/publications/RNAseq_6FT_peptides.bed));
	}
	$mode='ref_proteins';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0')
	{
		system(qq(cp $output_dir/intermediate/JUMPg_results/ref_proteins/reference_peptides.sc.bed $output_dir/publications/reference_peptides.bed));
	}
	combined_peptide_pubTab(\%parahash,$output_dir);
}

# simplified publication DIR
if (!-e "$output_dir/simplified_publications") {system(qq(mkdir $output_dir/simplified_publications));}
system(qq(cp $output_dir/publications/identified_peptide_table.txt $output_dir/simplified_publications));

=head
print <<EOF;
-------------------------------------------------------------------------------------
Preparing 
-------------------------------------------------------------------------------------
EOF
=cut
# for multistage
if (!-e "$output_dir/multistage") {system(qq(mkdir $output_dir/multistage));}
#if (-e "./qc_MSMS_input.txt") {system(qq(mv qc_MSMS_input.txt $output_dir/multistage));}
if (-e "./qc_MSMS_input.txt") {
	system(qq(mv qc_MSMS_input.txt $output_dir/multistage));
}
if (1) {
	my $mode='ref_proteins';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0') {
		expProt("$output_dir/intermediate/sum_accepted_PSMs/publications/id_uni_prot.txt",$parahash{'ref_proteins'},"$output_dir/multistage/identified_ref_proteins.fasta");
	}
	$mode='mutation';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0') {
		expMut("$output_dir/publications/mutation_events.txt","$output_dir/multistage/identified_mutations.txt");
	}
	$mode='junction';
	if (defined($parahash{$mode}) and $parahash{$mode} ne '0') {
		expJun("$output_dir/publications/junction_events.txt","$output_dir/multistage/identified_junctions.txt");
	}
}

print "\nJUMPg finished\n\n";

#----------------------------------------------------------------------------------------
sub expJun
{
	my ($events,$output)=@_;
	open(IN,$events) or die "Cannot open $events\n";
	open(OUT,'>',$output);
	my $line=<IN>;
	while (<IN>)
	{
		chomp;
		my ($junction,$strand,$others)=(split /\t/,$_);
		$strand =~ s/strand=//;
		# M.E. format
		# 13:67205417:+,13:67799537:+
		my @t=split /\:/,$junction;
		print OUT "chr$t[0]\t",$t[1]+1,"\t",$t[3]-1,"\t$strand\t$strand\t0\n";
	}
	close IN;
	close OUT;
}

sub expMut
{
	my ($events,$output)=@_;
	open(IN,$events) or die "Cannot open $events\n";
	open(OUT,'>',$output);
	my $line=<IN>;
	while (<IN>)
	{
		chomp;
		my $mutation=(split /\t/,$_)[0];
		# chr15:99669628-99669628.A_T
		$mutation =~ m/^(.*?)\:(\d+)\-(\d+)\.(.*?)\_(.*?)$/;
		print OUT "$1\t$2\t$3\t$4\t$5\n";
	}
	close IN;
	close OUT;
}

sub expProt
{
	my ($id_uni_prot,$dbFasta,$output)=@_;
	# parse id_uni_prot.txt, establish %prot
	my (%prot);
	expprot_parseTab($id_uni_prot,\%prot);
	# parse and filter fasta using %prot
	expprot_filterFas($dbFasta,\%prot,$output);
}

sub expprot_parseTab
{
        my ($tab,$prot)=@_;

        open(IN,$tab) or die "Cannot open $tab\n";
        while (<IN>)
        {
                chomp;
                next if (/^Unique/);
                next if (/^Protein/);

                my @t=split /\t/,$_;
                $$prot{$t[1]}='';
        }
        close IN;
}

sub expprot_filterFas
{
        my ($fasta,$prot,$output)=@_;

        open(IN,$fasta) or die "Cannot open $fasta\n";
        open(OUT,">$output");
        my ($pr,$seq,$title)=('','','');
        while (<IN>)
        {
                chomp;
                if (/>/)
                {
                        if ($pr ne '' and defined($$prot{$pr}))
                        {
                                print OUT "$title\n$seq\n";
                        }

                        $title=$_;
                        s/>//;
                        my @t=split / /,$_;
                        $pr=$t[0];
                        $seq='';
                }
                else
                { $seq.=$_; }
        }
        close IN;
}

sub combined_peptide_pubTab
{
	my ($parahash,$output_dir)=@_;
	my $idp=IDtxt_parser->new();
	my $pepxml=pepXML_parser->new();

	# parse pepXML
	my (%pxmlparahash,%acceptedOutHash);
	$pepxml->pepXML2Hashes(\%pxmlparahash,\%acceptedOutHash,"$output_dir/intermediate/sum_accepted_PSMs/html/accepted_PSM.pepXML");

	my ($ref_scanhash,$ref_peptidehash,$mut_scanhash,$mut_peptidehash,$jun_scanhash,$jun_peptidehash,$rdt_scanhash,$rdt_peptidehashi,$ast_scanhash,$ast_peptidehash);

	# parse IDtxt
	my $mode='ref_proteins';
	if (defined($$parahash{$mode}) and $$parahash{$mode} ne '0')
	{
		($ref_scanhash,$ref_peptidehash)=$idp->parse_IDtxt("$output_dir/intermediate/sum_accepted_PSMs/ID.txt");
		# rm non-ref peptides
		foreach my $pep (keys %$ref_peptidehash)
		{
			my $refpep=0;
			foreach my $pro (keys %{$ref_peptidehash->{$pep}->{proteins}})
			{
				unless ($pro =~ m/^cu\|/) { $refpep=1; last; }
			}
			if ($refpep) {}
			else { delete $ref_peptidehash->{$pep}; }
		}
	}
	$mode='mutation';
	if (defined($$parahash{$mode}) and $$parahash{$mode} ne '0')
	{
		($mut_scanhash,$mut_peptidehash)=$idp->parse_IDtxt("$output_dir/intermediate/JUMPg_results/mutation/valideted_ID.txt");
	}
	$mode='junction';
	if (defined($$parahash{$mode}) and $$parahash{$mode} ne '0')
	{
		($jun_scanhash,$jun_peptidehash)=$idp->parse_IDtxt("$output_dir/intermediate/JUMPg_results/junction/valideted_ID.txt");
	}
	$mode='transcript_6FT';
	if (defined($$parahash{$mode}) and $$parahash{$mode} ne '0') {
		($ast_scanhash,$ast_peptidehash)=$idp->parse_IDtxt("$output_dir/intermediate/JUMPg_results/transcript_6FT/novelscan.txt");
	}
	$mode='6FT';
	if (defined($$parahash{$mode}) and $$parahash{$mode} ne '0')
	{
		#($rdt_scanhash,$rdt_peptidehash)=$idp->parse_IDtxt("$output_dir/intermediate/JUMPg_results/6FT/valideted_ID.txt");
	}

	# print 
	my $output="$output_dir/publications/identified_peptide_table.txt";
	open(OUT,">$output");
	#print OUT "peptide\tscan_counts\tbest_scan\tsearch_score\tdelta_score\tsequence_ID\ttype\n";
	print OUT "Peptide\t\# PSMs\tBest PSM\tTag\tJscore\tdJn\tDatabase entry\tDatabase\n";
	close OUT;

	my %stahash;
	($stahash{ref}{pepN},$stahash{ref}{SC})=print_peptidehash('reference',$ref_peptidehash,\%acceptedOutHash,$output);
	($stahash{mut}{pepN},$stahash{mut}{SC})=print_peptidehash('mutation',$mut_peptidehash,\%acceptedOutHash,$output);
	($stahash{jun}{pepN},$stahash{jun}{SC})=print_peptidehash('junction',$jun_peptidehash,\%acceptedOutHash,$output);
	($stahash{ast}{pepN},$stahash{ast}{SC})=print_peptidehash('transcript_6FT',$ast_peptidehash,\%acceptedOutHash,$output);
	#($stahash{rdt}{pepN},$stahash{rdt}{SC})=print_peptidehash('6FT',$rdt_peptidehash,\%acceptedOutHash,$output);

	# print summary statistics
	print "Type\t\tPSMs\tPeptides\n";
	print "Reference\t$stahash{ref}{SC}\t$stahash{ref}{pepN}\n";
	print "Mutation\t$stahash{mut}{SC}\t$stahash{mut}{pepN}\n";
	print "Junction\t$stahash{jun}{SC}\t$stahash{jun}{pepN}\n";
	print "RNA 6FT\t\t$stahash{ast}{SC}\t$stahash{ast}{pepN}\n";
	#print "RNAseq 6FT\t$stahash{rdt}{SC}\t$stahash{rdt}{pepN}\n";
	
	# print scanhash
	#print_scanhash('transcript_6FT',$ast_scanhash,\%acceptedOutHash,"$output_dir/publications/identified_PSMs_table.txt");
}

sub print_scanhash {
	my ($mode,$scanhash,$acceptedOutHash,$output)=@_;

	open(OUT,">>$output");
	foreach my $out (keys %$scanhash) {
		my $outfile=basename($out);
		$outfile =~ s/\.[spout]+$//;
		my $tag=(defined($$acceptedOutHash{$outfile}{TagSeq}[1]))?$$acceptedOutHash{$outfile}{TagSeq}[1]:'N/A';
		print OUT "$outfile\t$tag\t$scanhash->{$out}->{XCorr}\t$scanhash->{$out}->{dCn}\t$mode\n";
	}
	close OUT;
}

sub print_peptidehash
{
	my ($mode,$peptidehash,$acceptedOutHash,$output)=@_;

	my ($pepN,$SC)=(0,0);
	open(OUT,">>$output");
	foreach my $pep (keys %$peptidehash)
	{
		my $protein=(keys %{$peptidehash->{$pep}->{proteins}})[0];
		my $outfile=basename($peptidehash->{$pep}->{top_score_outfile});
		$outfile =~ s/\.[spout]+$//;
		my $tag=(defined($$acceptedOutHash{$outfile}{TagSeq}[1]))?$$acceptedOutHash{$outfile}{TagSeq}[1]:'N/A';
		print OUT "$pep\t",scalar(keys %{$peptidehash->{$pep}->{outfiles}}),"\t$outfile\t$tag\t$peptidehash->{$pep}->{XCorr}\t$peptidehash->{$pep}->{dCn}\t$protein\t$mode\n";
		$pepN++;
		$SC+=scalar(keys %{$peptidehash->{$pep}->{outfiles}});
	}
	close OUT;

	return ($pepN,$SC);
}

sub params_pa_update
{
	# set 'mode'
	# update 'input_file', 'species'
	# set 'output_folder' as 'mode'
	my ($parahash,$mode,$output_dir,$currentdir,$bname)=@_;
	my $output.="$output_dir/JUMPg_results/jump_g_pa.params";
	open(OUT,">$output");

	# decide core_file_name
	my $core_file_name; #=get_core_file_name($parahash,$mode);
	if ($mode eq 'mutation') 
	{
		$core_file_name="$currentdir/customizedDB/mutation/$bname";
	}
	elsif ($mode eq 'junction') 
	{
		#$bname =~ s/\.[^.]+$//;
		$core_file_name="$currentdir/customizedDB/junction/$bname";
	}
	elsif ($mode eq 'transcript_6FT') {
		$core_file_name="$currentdir/customizedDB/transcript_6FT/$bname\_formated.fasta";
	}
	elsif ($mode eq '6FT') 
	{
		#$bname =~ s/\.[^.]+$//;
		$core_file_name="$currentdir/customizedDB/6FT/$bname";
		print OUT "gene_ID_convert_file = $$parahash{gene_ID_convert_file}\n";
	}
	elsif ($mode =~ m/^ref/)
	{
		$core_file_name="";
		if (-e "$$parahash{gene_ID_convert_file}_updated") {
			print OUT "gene_ID_convert_file = $$parahash{gene_ID_convert_file}_updated\n";
		} else {
			print OUT "gene_ID_convert_file = $$parahash{gene_ID_convert_file}\n";
		}
	}

	# decide ID_txt
	my $ID_txt="$output_dir/JUMPg_results/$mode\_ID.txt";;
	#if ($mode eq 'mutation') {$ID_txt="$output_dir/JUMPg_results/mut_ID.txt";}
	#elsif ($mode eq 'junction') {}
	#elsif ($mode eq '6FT') {}

	# print
	print OUT "mode = $mode\n";
	print OUT "core_file_name = $core_file_name\n";
	print OUT "ID_txt = $ID_txt\n";
	print OUT "output_folder = $mode\n";
	print OUT "track_name = $mode\n";
	#print OUT "species = $$parahash{species}\n";
	print OUT "annovar_annotation_file = $$parahash{annovar_annotation_file}\n";
	print OUT "reference_genome = $$parahash{reference_genome}\n";
	#print OUT "reference_mRNA = $$parahash{reference_mRNA}\n";
	#print OUT "reference_CDS = $$parahash{reference_CDS}\n";
	#print OUT "reference_gene_locus = $$parahash{reference_gene_locus}\n";

	#if ($mode eq '6FT') {
	if (defined($$parahash{reference_protein})) {
		print OUT "ref_proteins = $$parahash{reference_protein}\n";
	}

	close OUT;

	return $output;
}

sub inputRawFiles
{
	my ($input)=@_;
	$input =~ /\.([^.]+)$/;
	my $sf=uc($1);
	my $raw=($sf =~ /RAW/ or $sf =~ /MZXML/)?1:0;
	return $raw;
}

sub print2file
{
	my ($output,@lines)=@_;
	open(OUT,">>$output");
	for (my $i=0;$i<=$#lines;$i++) {print OUT $lines[$i];}
	close OUT;
}

sub get_mutation_fas_name
{
        my ($parahash,$currentdir)=@_;
	my $file=basename($$parahash{mutation});
	return "$currentdir/customizedDB/mutation/$file.seqpad.mut.uniq.newIDs.fas";
}

sub get_junction_fas_name
{
        my ($parahash,$currentdir)=@_;
	my $file=basename($$parahash{junction});
	#$file =~ s/\.[^.]+$//;
	#return "$currentdir/customizedDB/junction/$file.AA.central.newIDs.fas";
	return "$currentdir/customizedDB/junction/$file.junc.AA.central.newIDs.fas";
}

sub get_transcript_6FT_fas_name {
	my ($parahash,$currentdir)=@_;
	my $file=basename($$parahash{transcript_6FT});
	$file =~ s/\.[^.]+$//;
	return "$currentdir/customizedDB/transcript_6FT/$file\_peptides.newIDs.fas";
}

sub get_6FT_fas_name
{
        my ($parahash,$currentdir)=@_;
	my $file=basename($$parahash{'6FT'});
	$file =~ s/\.[^.]+$//;
	my @t;
	$t[0]=($$parahash{read_coverage_cutoff}>2)?$$parahash{read_coverage_cutoff}:2;
	return "$currentdir/customizedDB/6FT/$file\_QC_peptide_digested_uniq_$t[0]reads.newIDs.fas";
}

sub set_summary_params
{
	my ($parahash,$currentdir)=@_;

	my ($mut_ids, $jun_ids, $rdt_ids)=('na','na','na');
	if ($$parahash{mutation} ne '0')
	{
		$mut_ids=get_mutation_fas_name($parahash,$currentdir);
		$mut_ids =~ s/fas$/ids/;
	}
	if ($$parahash{junction} ne '0')
	{
		$jun_ids=get_junction_fas_name($parahash,$currentdir);
		$jun_ids =~ s/fas$/ids/;
	}
	if ($$parahash{transcript_6FT} ne '0')
	{
		$rdt_ids=get_transcript_6FT_fas_name($parahash,$currentdir);
		$rdt_ids =~ s/fas$/ids/;
	}
=head
	if ($$parahash{'6FT'} ne '0')
	{
		$rdt_ids=get_6FT_fas_name($parahash,$currentdir);
		$rdt_ids =~ s/fas$/ids/;
	}
=cut
	return ($mut_ids, $jun_ids, $rdt_ids);
}

sub params_qc_update
{
        my ($input,$conf_f,$acpt_f,$PSM_recoveray_rate,$outname,$output)=@_;

        # read -s result paths
        my %updatedHash;
        open(IN,$input) or die "Cannot open $input!!!";
        my @sResults=<IN>;
        close IN;

        # write -qc params
        open(OUT,'>',$output);
        for (my $i=0; $i<=$#sResults;$i++)
        {
                print OUT "run",$i+1," = $sResults[$i]";
        }
        print OUT "confident_IDtxt = $conf_f\n";
        print OUT "accepted_IDtxt = $acpt_f\n";
        print OUT "PSM_recoveray_rate = $PSM_recoveray_rate\n";
        print OUT "output_folder = $outname\n";
        close OUT;
}


sub set_params_files
{
	my ($parahash)=@_;

	if ($$parahash{search_engine} eq 'JUMP')
	{
		$$parahash{'params'}{s}="jump_sj_$$parahash{MS_data_types}.params";
		$$parahash{'params'}{f}="jump_fj_$$parahash{MS_data_types}.params";
	}
	elsif ($parahash{search_engine} eq 'SEQUEST')
	{
		$$parahash{'params'}{s}="jump_ss_$$parahash{MS_data_types}.params";
		$$parahash{'params'}{f}="jump_fs_$$parahash{MS_data_types}.params";
	}
	
}

sub params_f_update
{
	my ($input,$param,$outname,$output)=@_;

	# read -s result paths
	my %updatedHash;
        open(IN,$input) or die "Cannot open $input!!!";
        my @sResults=<IN>;
        close IN;

        # read -f params
        open(IN,$param) or die "Cannot open $param!!!\n";
        my @lines=<IN>;
        close IN;

        # re-write -f params
        open(OUT,'>',$output);
        print OUT "$outname\: $sResults[0]\n";
        for (my $i=1; $i<=$#sResults;$i++)
        {
                print OUT "$sResults[$i]\n";
        }


        for (my $i=0; $i<=$#lines;$i++)
        {
                $lines[$i] =~ s/^\s+//; # rm space at front
                next if ($lines[$i] =~ /^#/); # skip comment lines
                #chomp;

                if ( $lines[$i]  =~ / = /  )
                {print OUT $lines[$i];}
        }
        close OUT;

}

sub add_fast_run_label
{
	my ($p)=@_;

	open(IN,$p) or die "Cannot open $p!!!\n";
	my @lines=<IN>;
	close IN;

	unshift @lines,"fast_run = 1\n";

	open(OUT,'>',$p);
	for (my $i=0; $i<=$#lines;$i++) {print OUT $lines[$i];}
	close OUT;
}

sub params_s_update
{
        my ($input,$param_s)=@_;

        my %updatedHash;
        open(IN,$input) or die "Cannot open $input!!!";
        $updatedHash{database_name}=<IN>; chomp($updatedHash{database_name});
        $updatedHash{pit_file}=<IN>; chomp($updatedHash{pit_file});
        close IN;

        params_update_by_hash($param_s,$param_s,\%updatedHash);
}


sub params_d_update
{
	# mkdir $output_dir/database
	# set @database[1]=$$parahash{ref_proteins}
	# read customizedDB/.jump_c_tmp => push @database
	# read $output_dir/ParameterFiles/jump_d.params: set input_database, output_prefix, jump.params
	#
	my ($parahash,$output_dir,$currentdir)=@_;

	# mkdir $output_dir/database
	if (!-e "$output_dir/database") { system(qq(mkdir $output_dir/database)); }

	# define @database
	my @database;
	# set @database[1]=$$parahash{ref_proteins}
	if (defined($$parahash{ref_proteins}) and $$parahash{ref_proteins} ne '0')
	{ $database[0]=$$parahash{ref_proteins}; }
	# read customizedDB/.jump_c_tmp => push @database
	if (-e "$currentdir/customizedDB/.jump_c_tmp")
	{
		open(IN,"$currentdir/customizedDB/.jump_c_tmp") || die "Cannot open $currentdir/customizedDB/.jump_c_tmp!!!\n";
		while (<IN>)
		{
			chomp;
			push @database, $_;
		}
		close IN;
	}

	# read $output_dir/ParameterFiles/jump_d.params: set input_database, output_prefix,jump.params
	my $output="$output_dir/database/jump_d.params";
	open(IN,"$output_dir/ParameterFiles/jump_d.params") || die "Cannot open $output_dir/ParameterFiles/jump_d.params!!!\n";
	open(OUT,">$output");
	my $dbmark=0; # mark if DB rows already print out
	while (<IN>)
	{
		if (/^input_database/)
		{
			next if ($dbmark);
			for (my $i=0;$i<=$#database;$i++)
			{
				print OUT "input_database",$i+1," = $database[$i]\n";
			}
			$dbmark=1;
		}
		elsif (/^output_prefix/)
		{
			print OUT "output_prefix = $$parahash{output_directory}\n";
		}
		elsif (/^jump.params/)
		{
			#print OUT "jump.params = $$parahash{'params'}{s}\n";
			print OUT "jump.params = $output_dir/$$parahash{'params'}{s}\n";
		#} elsif (/^include_contaminants/ and defined($$parahash{include_contaminants})) {
		#	print OUT "include_contaminants = $$parahash{include_contaminants}\n";
		}
		#else {print OUT $_;}
	}

	print OUT "include_contaminants = 0\n";
	print OUT "decoy_generation = 1\n";
	print OUT "decoy_generation_method = 1\n";
	print OUT "bypass_db_generation = 0\n";


	close OUT;
	close IN;

	return $output;
}

sub params_update_by_hash
{
        my ($input,$output,$updatedHash)=@_;

        open(IN,$input) or die "Cannot open $input!!!\n";
        my @lines=<IN>;
        close IN;

        open(OUT,'>',$output);
        for (my $i=0; $i<=$#lines;$i++)
        {
                my @t=split(' = ',$lines[$i]);
                if (defined($$updatedHash{$t[0]}))
                {
                        print OUT "$t[0] = $$updatedHash{$t[0]}\n";
                }
                else {print OUT $lines[$i];}
        }
        close OUT;
}

sub params_c_update
{
	# read ./genm_userDefinedSuffix/ParameterFiles/jump_c.params
	# set 'mode'
	# update 'input_file', 'species'
	# set 'output_folder' as 'mode'
	my ($parahash,$mode,$currentdir)=@_;
	#system(qq(cp ./genm_userDefinedSuffix/ParameterFiles/jump_c.params customizedDB));
	my $output=$currentdir;
	$output.='/customizedDB/jump_c.params';
	open(OUT,">$output");
	print OUT "mode = $mode\n";
	print OUT "input_file = $$parahash{$mode}\n";
	print OUT "output_folder = $mode\n";
	#print OUT "species = $$parahash{species}\n";
	print OUT "annovar_annotation_file = $$parahash{annovar_annotation_file}\n";
	print OUT "reference_genome = $$parahash{reference_genome}\n";
	if (defined($$parahash{read_coverage_cutoff})) { 
		print OUT "read_coverage_cutoff = $$parahash{read_coverage_cutoff}\n"; 
	} else {
		print OUT "read_coverage_cutoff = 2\n";
	}
	if (defined($$parahash{min_ORF_AA})) {
		print OUT "min_ORF_AA = $$parahash{min_ORF_AA}\n";
	} else {
		print OUT "min_ORF_AA = 66\n";
	}
	close OUT;

	return $output;
}

sub set_program_paths
{
	my ($parahash)=@_;

	#my $path='/home/yli4/development/JUMPg/v2.0/programs/';
	my $path="/home/yli4/development/JUMPg/JUMPg_v2.3.5/programs";
	#my $path=;
	#my $path=aa";
	$$parahash{'programs'}{params}="$path/params/jump_params.pl";
	$$parahash{'programs'}{c}="$path/c/jump_c.pl";
	$$parahash{'programs'}{d}="$path/d/builddb.pl";
	$$parahash{'programs'}{s}="$path/s/jump.pl";
	$$parahash{'programs'}{f}="$path/f/jump_f.pl";
	$$parahash{'programs'}{qc}="$path/qc/spectrumQC.pl";
	$$parahash{'programs'}{summary}="$path/summary/splitIDtxt.pl";
	$$parahash{'programs'}{rundtas}="$path/rundtas/rundtas.pl";
	$$parahash{'programs'}{pa}="$path/pa/jump_g_postAnnotations.pl";
}

sub check_parameter_version
{
	my ($par,$VERSION)=@_;

	open(IN,$par);
	my $line=<IN>;
	$line =~ m/Version\:\s+(\d+)\.(\d+)\.(\d+)\,/;
	my $v="$1\.$2\.$3";
	close IN;

	if ($VERSION ne $v)
	{
		die "\nInconsistent parameter file version (detected as $v).\n\nPlease use version $VERSION parameter file.\n\n";
	}
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
			$key =~ s/^input_//;
			# temporary solution
			if ($key eq 'RNAseq_6FT') { $key='6FT'; }

                        $$parahash{$key}=$value;
                }
        }
        close IN;

	$$parahash{MS_data_types}='HH';
	$$parahash{search_engine}='JUMP';
	if ($$parahash{MS_data_types} =~ m/TMT/ ) { $$parahash{TMT_data}=1; }
	if (!defined($$parahash{ref_proteins})) { $$parahash{ref_proteins}=0; }
	if (!defined($$parahash{mutation})) { $$parahash{mutation}=0; }
	if (!defined($$parahash{junction})) { $$parahash{junction}=0; }
	if (!defined($$parahash{'6FT'})) { $$parahash{'6FT'}=0; }

}

sub help
{

        print  "\nUsage: perl JUMPg.pl JUMPg.params *.raw/mzXML\n\n";
	print  "Or\n";
        print  "\nUsage: perl JUMPg.pl JUMPg.params qc_MSMS_input.txt\n\n";
	exit;
}


