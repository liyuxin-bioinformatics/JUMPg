#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Basename;

our $VERSION = 			'$Revision: 529 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2013-08-23 08:53:41 -0700 (Fri, 23 Aug 2013) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $remove, $checkfile, $protocol, $operation, $otherinfo, $onetranscript, $nastring, $genericdbfile, $gff3dbfile, $bedfile, $vcfdbfile, $csvout, $argument);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'remove'=>\$remove, 'checkfile!'=>\$checkfile, 
	'protocol=s'=>\$protocol, 'operation=s'=>\$operation, 'otherinfo'=>\$otherinfo, 'onetranscript'=>\$onetranscript, 'nastring=s'=>\$nastring,
	'genericdbfile=s'=>\$genericdbfile, 'gff3dbfile=s'=>\$gff3dbfile, 'bedfile=s'=>\$bedfile, 'vcfdbfile=s'=>\$vcfdbfile,
	'csvout'=>\$csvout, 'argument=s'=>\$argument
	) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($queryfile, $dbloc) = @ARGV;


my @unlink;					#a list of files to be deleted
my @header;					#annotation output header
my %varanno;					#varstring as key of 1st hash, anno_db as key of 2nd hash, anno as value of 2nd hash
my (@protocol, @operation, @dbtype1, @argument);		#@dbtype1 is the translated file name for @protocol
my (@genericdbfile, @gff3dbfile, @vcfdbfile, @bedfile);
my (@genericdbfile1, @gff3dbfile1, @vcfdbfile1, @bedfile1);

processArgument ();
@dbtype1 = proxyDBType(@protocol);
$checkfile and checkFileExistence (@dbtype1);

for my $i (0 .. @protocol-1) {
	print STDERR "-----------------------------------------------------------------\n";
	print STDERR "NOTICE: Processing operation=$operation[$i] protocol=$protocol[$i]\n";
	if ($operation[$i] eq 'g') {
		geneOperation ($protocol[$i], $dbtype1[$i], $argument[$i]||undef);
	} elsif ($operation[$i] eq 'r') {
		regionOperation ($protocol[$i], $dbtype1[$i], $argument[$i]||undef);
	} elsif ($operation[$i] eq 'f') {
		filterOperation ($protocol[$i], $dbtype1[$i], $argument[$i]||undef);
	}
}

printOrigOutput ();


$remove and unlink(@unlink);



sub printOrigOutput {
	my $final_out;
	if ($csvout) {
		$final_out="$outfile.${buildver}_multianno.csv";
	} else {
		$final_out="$outfile.${buildver}_multianno.txt";
	}
	open OUT,">",$final_out or die "Cannot write to $final_out: $!\n";
	print STDERR "NOTICE: output file is written to $final_out\n";
	
	my @expanded_header;
	for my $i (0 .. @header-1) {
		if ($header[$i] eq 'ljb_all') {
			push @expanded_header, qw/LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++/;
		} elsif ($header[$i] eq 'ljb2_all') {
			push @expanded_header, qw/LJB2_SIFT LJB2_PolyPhen2_HDIV LJB2_PP2_HDIV_Pred LJB2_PolyPhen2_HVAR LJB2_PolyPhen2_HVAR_Pred LJB2_LRT LJB2_LRT_Pred LJB2_MutationTaster LJB2_MutationTaster_Pred LJB_MutationAssessor LJB_MutationAssessor_Pred LJB2_FATHMM LJB2_GERP++ LJB2_PhyloP LJB2_SiPhy/;
		} elsif ($header[$i] eq 'popfreq_all') {
			push @expanded_header, qw/PopFreqMax 1000G2012APR_ALL 1000G2012APR_AFR 1000G2012APR_AMR 1000G2012APR_ASN 1000G2012APR_EUR ESP6500si_AA ESP6500si_EA CG46 NCI60 SNP137 COSMIC65 DISEASE/;
		} else {
			push @expanded_header, $header[$i];
		}
	}
	
	if ($csvout) {
		print OUT join(",", qw/Chr Start End Ref Alt/, @expanded_header), $otherinfo?",Otherinfo":"", "\n";
	} else {
		print OUT join("\t", qw/Chr Start End Ref Alt/, @expanded_header), $otherinfo?"\tOtherinfo":"", "\n";
	}
	
	open (FH, $queryfile) or die "Error: cannot read from inputfile $queryfile: $!\n";
	while (<FH>) {
		s/[\r\n]+$//;
		m/^(\S+\s+\S+\s+\S+\s+\S+\s+\S+)\s*(.*)/ or next;
		my ($varstring, $info) = ($1, $2);
		my @varstring = split (/\s+/, $varstring);
		$varstring =~ s/\s+/\t/g;
		
		my @oneline;
		for my $item (@header) {
			my $expanded_field;
			if ($item eq 'ljb_all') {
				$expanded_field = 11;
			} elsif ($item eq 'ljb2_all') {
				$expanded_field = 15;
			} elsif ($item eq 'popfreq_all') {
				$expanded_field = 13;
			}
			if ($csvout) {
				if ($expanded_field) {
					push @oneline,($varanno{$varstring}{$item} || join (",", ($nastring) x $expanded_field));	#value is already comma-delimited
				} else {
					push @oneline,(defined $varanno{$varstring}{$item}) ? qq/"$varanno{$varstring}{$item}"/ : $nastring;	#value of $varanno{$varstring}{$item} could be zero
				}
			} else {
				if ($expanded_field) {
					if ($varanno{$varstring}{$item}) {
						for my $nextscore (split (/,/, $varanno{$varstring}{$item})) {
							push @oneline, $nextscore;
						}
					} else {
						for (1 .. $expanded_field) {
							push @oneline, $nastring;
						}
					}
				} else {
					push @oneline, (defined $varanno{$varstring}{$item})? $varanno{$varstring}{$item} : $nastring;	#value of $varanno{$varstring}{$item} could be zero
				}
			}
		}
		if ($csvout) {
			print OUT join (",", @varstring, @oneline), $otherinfo?qq/,"$info"/:"", "\n";
		} else {
			print OUT join ("\t", @varstring, @oneline), $otherinfo?"\t$info":"", "\n";
		}
	}
}

sub processArgument {
	$outfile ||= $queryfile;
	
	if (not defined $buildver) {
		$buildver = 'hg18';
		print STDERR "NOTICE: the --buildver argument is set as 'hg18' by default\n";
	}
	
	not defined $checkfile and $checkfile = 1;
	
	defined $nastring or $nastring = '';
	
	if (not $protocol) {
		$operation and pod2usage ("Error in argument: you must specify --protocol if you specify --operation");
		if ($buildver eq 'hg18') {
			$protocol = 'gene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb_all';
			$operation = 'g,r,r,f,f,f,f';
			print STDERR "NOTICE: the --protocol argument is set as 'gene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all' by default\n";
		} elsif ($buildver eq 'hg19') {
			$protocol = 'gene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb_all';
			$operation = 'g,r,r,f,f,f,f';
			print STDERR "NOTICE: the --protocol argument is set as 'gene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all' by default\n";
		} else {
			pod2usage ("Error in argument: please specify --protocol argument for the --buildver $buildver");
		}
	}
	
	if ($protocol =~ m/\bgeneric\b/) {
		$genericdbfile or pod2usage ("Error in argument: please specify -genericdbfile argument when 'generic' operation is specified");
	}
	
	#additional preparation work
	@protocol = split (/,/, $protocol);
	@operation = split (/,/, $operation);
	@argument = split (/,/, $argument||'', -1);
	@protocol == @operation or pod2usage ("Error in argument: different number of elements are specified in --protocol and --operation argument");
	@argument and @protocol == @argument || pod2usage ("Error in argument: different number of elements are specified in --protocol and --argument argument");
	for my $op (@operation) {
		$op =~ m/^g|r|f$/ or pod2usage ("Error in argument: the --operation argument must be comma-separated list of 'g', 'r', 'f'");
	}
		
	my %uniq_protocol;
	for (@protocol) {
		$uniq_protocol{$_}++;
	}
	if ($uniq_protocol{generic}) {
		$genericdbfile or pod2usage ("Error in argument: please specify -genericdbfile argument when 'generic' protocol is specified");
		@genericdbfile = split (/,/, $genericdbfile);
		@genericdbfile == $uniq_protocol{generic} or pod2usage ("Error in argument: you specified $uniq_protocol{generic} 'generic' in 'protocol' argument, but only ${\(scalar @genericdbfile)} filenames in 'genericdbfile' argument");
	}
	if ($uniq_protocol{gff3}) {
		$gff3dbfile or pod2usage ("Error in argument: please specify -gff3dbfile argument when 'gff3' protocol is specified");
		@gff3dbfile = split (/,/, $gff3dbfile);
		@gff3dbfile == $uniq_protocol{gff3} or pod2usage ("Error in argument: you specified $uniq_protocol{gff3} 'gff3' in 'protocol' argument, but only ${\(scalar @gff3dbfile)} filenames in 'gff3dbfile' argument");
	}
	if ($uniq_protocol{vcf}) {
		$vcfdbfile or pod2usage ("Error in argument: please specify -vcfdbfile argument when 'vcf' protocol is specified");
		@vcfdbfile = split (/,/, $vcfdbfile);
		@vcfdbfile == $uniq_protocol{vcf} or pod2usage ("Error in argument: you specified $uniq_protocol{vcf} 'vcf' in 'protocol' argument, but only ${\(scalar @vcfdbfile)} filenames in 'vcfdbfile' argument");
	}
	if ($uniq_protocol{bed}) {
		$bedfile or pod2usage ("Error in argument: please specify -bedfile argument when 'bed' protocol is specified");
		@bedfile = split (/,/, $bedfile);
		@bedfile == $uniq_protocol{bed} or pod2usage ("Error in argument: you specified $uniq_protocol{bed} 'bed' in 'protocol' argument, but only ${\(scalar @bedfile)} filenames in 'bedfile' argument");
	}

	#prepare PATH environmental variable
	my $path = File::Basename::dirname ($0);
	$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in
}

#gene annotations
sub geneOperation {
	my ($protocol, $dbtype1, $argument) = @_;
	#generate gene anno
	my $sc;
	$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $protocol -outfile $outfile.$protocol -exonsort $queryfile $dbloc";
	$argument and $sc .= " $argument";
	system ($sc) and die "Error running system command: <$sc>\n";
	
	#read in gene anno
	my $anno_outfile="$outfile.$protocol.variant_function";
	my $e_anno_outfile="$outfile.$protocol.exonic_variant_function";
	
	open (FUNCTION, "<",$anno_outfile) or die "Error: cannot read from $anno_outfile: $!\n";	
	open (EFUNCTION,"<",$e_anno_outfile) or die "Error: cannot read from $e_anno_outfile: $!\n";
	
	push @header,"Func.$protocol", "Gene.$protocol", "ExonicFunc.$protocol", "AAChange.$protocol"; #header
	while (<FUNCTION>) 
	{
		s/[\r\n]+$//;
		m/^(\S+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+).*/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
		my ($function, $gene, $varstring) = ($1, $2, $3);
		$varstring =~ s/\s+/\t/g;
		$varanno{$varstring}{"Func.$protocol"}=$function;
		$varanno{$varstring}{"Gene.$protocol"}=$gene;
	}
	close FUNCTION;
	while (<EFUNCTION>)
	{
		m/^line\d+\t([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
		my ($efunc, $aachange, $varstring) = ($1, $2, $3);
		$varstring =~ s/\s+/\t/g;
		my @aachange = split (/:|,/, $aachange);

		$varanno{$varstring}{"ExonicFunc.$protocol"}=$efunc;
		if (not $onetranscript) {
			$varanno{$varstring}{"AAChange.$protocol"}=$aachange;
		} else {
			if (@aachange >= 5) 
			{
			    $varanno{$varstring}{"AAChange.$protocol"}="$aachange[1]:$aachange[3]:$aachange[4]"; #only output aachange in first transcript
			} else {
			    $varanno{$varstring}{"AAChange.$protocol"}=$aachange;		#aachange could be "UNKNOWN"
			}
		}
	}
	close EFUNCTION;
	push @unlink, $anno_outfile, $e_anno_outfile, "$outfile.$protocol.log";
}

sub regionOperation {
	my ($dbtype, $dbtype1, $argument) = @_;
	my ($userfile, $header);
	my $sc = "annotate_variation.pl -regionanno -dbtype $dbtype -buildver $buildver -outfile $outfile $queryfile $dbloc";
	$argument and $sc .= " $argument";
	
	if ($dbtype eq 'bed') {
		$userfile = shift @bedfile1;
		$sc .= " -bedfile $userfile";
	} elsif ($dbtype eq 'gff3') {
		$userfile = shift @gff3dbfile1;
		$sc .= " -gff3dbfile $userfile";
	}
	
	print STDERR "\nNOTICE: Running with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";

	if ($dbtype eq 'bed') {
		if (@bedfile-@bedfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@bedfile-@bedfile1);
		}
	} elsif ($dbtype eq 'gff3') {
		if (@gff3dbfile-@gff3dbfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@gff3dbfile-@gff3dbfile1);
		}
	} else {
		$header = $dbtype;
	}
	push @header,$header;
	
	open (FH, "$outfile.${buildver}_$dbtype1") or die "Error: cannot open file\n";
	while (<FH>) 
	{
		m/^([^\t]+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+).*/ or die "Error: invalid record found in annovar outputfile: \n";
		my ($db,$anno,$varstring)=($1,$2,$3);
		$varstring =~ s/\s+/\t/g;
		#preprocess $anno
		$varanno{$varstring}{$header}=$anno;
	}
	close (FH);
	
	push @unlink, "$outfile.${buildver}_$dbtype1", "$outfile.log";
}

sub filterOperation {
	my ($dbtype, $dbtype1, $argument) = @_;
	my ($userfile, $header);
	my $sc = "annotate_variation.pl -filter -dbtype $dbtype -buildver $buildver -outfile $outfile $queryfile $dbloc";
	$argument and $sc .= " $argument";
	
	if ($dbtype =~ m/^ljb\d*_all/ or $dbtype =~ m/^ljb\d*_mod/ or $dbtype =~ m/^popfreq/) {
		$sc .= " -otherinfo";
	} elsif ($dbtype eq 'generic') {
		$userfile = shift @genericdbfile1;
		$sc .= " -genericdbfile $userfile";
	} elsif ($dbtype eq 'vcf') {
		$userfile = shift @vcfdbfile1;
		$sc .= " -vcfdbfile $userfile";
	}
	
	if ($dbtype eq 'avsift') {
		$sc .= " -sift_threshold 0";		#for historical reasons
	}
	
	print STDERR "\nNOTICE: Running system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";

	if ($dbtype eq 'generic') {
		if (@genericdbfile-@genericdbfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@genericdbfile-@genericdbfile1);
		}
	} elsif ($dbtype eq 'vcf') {
		if (@vcfdbfile-@vcfdbfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@vcfdbfile-@vcfdbfile1);
		}
	} else {
		$header = $dbtype;
	}
	push @header,$header;
	
	open (FH, "$outfile.${buildver}_${dbtype1}_dropped") or die "Error: cannot open file\n";
	while (<FH>) 
	{
		m/^([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+).*/ or die "Error: invalid record found in annovar outputfile\n";
		my ($db,$anno,$varstring)=($1,$2,$3);
		$varstring =~ s/\s+/\t/g;
		$varanno{$varstring}{$header}=$anno;
	}
	close (FH);
	
	push @unlink, "$outfile.${buildver}_${dbtype1}_dropped", "$outfile.${buildver}_${dbtype1}_filtered", "$outfile.log";
}

sub proxyDBType {
	my @db_names = @_;
	#    my @db_for_check;
	
	#    map {push @db_for_check,"${buildver}_$_.txt" } @db_names;
	#to overcome case-sensitive matching problem, check file existence case-INsensitively
	#    my @all_db=glob File::Spec->catfile ($dbloc,"*");
	#    my $all_db_string = join (',',@all_db);
	
	#    for (@db_for_check)
	#    {
	#	$all_db_string =~ m/$_/i or die "Error: the required database file $_ does not exist. Please download it via -downdb argument by annotate_variation.pl.\n";
	#    }
	my %dbalias=(					#for backward compatibility (historical reasons)
		'gene'=>'refGene', 
		'refgene'=>'refGene', 
		'knowngene'=>'knownGene', 
		'ensgene'=>'ensGene', 
		'band'=>'cytoBand', 
		'cytoband'=>'cytoBand', 
		'tfbs'=>'tfbsConsSites', 
		'mirna'=>'wgRna',
		'mirnatarget'=>'targetScanS', 
		'segdup'=>'genomicSuperDups', 
		'omimgene'=>'omimGene', 
		'gwascatalog'=>'gwasCatalog',
		'1000g_ceu'=>'CEU.sites.2009_04', 
		'1000g_yri'=>'YRI.sites.2009_04', 
		'1000g_jptchb'=>'JPTCHB.sites.2009_04',
		'1000g2010_ceu'=>'CEU.sites.2010_03', 
		'1000g2010_yri'=>'YRI.sites.2010_03', 
		'1000g2010_jptchb'=>'JPTCHB.sites.2010_03',
		'1000g2010jul_ceu'=>'CEU.sites.2010_07', 
		'1000g2010jul_yri'=>'YRI.sites.2010_07', 
		'1000g2010jul_jptchb'=>'JPTCHB.sites.2010_07',
		'1000g2010nov_all'=>'ALL.sites.2010_11', 
		'1000g2011may_all'=>'ALL.sites.2011_05'
	);
    
	for (@db_names) {
		$_ = $dbalias{$_} || $_;
		if (m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
			my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
			$_ = uc($3) . ".sites." . $1 . '_' . $monthhash{$2};
			$monthhash{$2} or die "Error: the '$_' operation does not contain a valid month abbreviation\n";
		}
	}
	for (@db_names) {
		m/^mce(\d+way)$/ and $_ = "phastConsElements$1";		#for backward compatibility in ANNOVAR
	}
	return @db_names;
}

sub checkFileExistence {    
	my @db_names = @_;
	@genericdbfile1 = @genericdbfile;
	@vcfdbfile1 = @vcfdbfile;
	@bedfile1 = @bedfile;
	@gff3dbfile1 = @gff3dbfile;
	for my $i (0 .. @db_names-1) {
		my $dbfile;
		if ($db_names[$i] eq 'generic') {
			$dbfile = File::Spec->catfile ($dbloc, shift @genericdbfile1);
		} elsif ($db_names[$i] eq 'vcf') {
			$dbfile = File::Spec->catfile ($dbloc, shift @vcfdbfile1);
		} elsif ($db_names[$i] eq 'bed') {
			$dbfile = File::Spec->catfile ($dbloc, shift @bedfile1);
		} elsif ($db_names[$i] eq 'gff3') {
			$dbfile = File::Spec->catfile ($dbloc, shift @gff3dbfile1);
		} else {
			$dbfile = File::Spec->catfile ($dbloc, "${buildver}_$db_names[$i].txt");
		}
		-f $dbfile or die "Error: the required database file $dbfile does not exist.\n";
	}
	@genericdbfile1 = @genericdbfile;
	@vcfdbfile1 = @vcfdbfile;
	@bedfile1 = @bedfile;
	@gff3dbfile1 = @gff3dbfile;
}

=head1 SYNOPSIS

 table_annovar.pl [arguments] <query-file> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --protocol <string>		comma-delimited string specifying database protocol
            --operation <string>	comma-delimited string specifying type of operation
            --outfile <string>		output file name prefix
            --buildver <string>		genome build version (default: hg18)
            --remove			remove all temporary files
            --(no)checkfile		check if database file exists (default: ON)
            --genericdbfile <files>	specify comma-delimited generic db files
            --gff3dbfile <files>	specify comma-delimited GFF3 files
            --bedfile <files>		specify comma-delimited BED files
            --vcfdbfile <files>		specify comma-delimited VCF files
            --otherinfo			print out otherinfo (infomration after fifth column in queryfile)
            --onetranscript		print out only one transcript for exonic variants (default: all transcripts)
            --nastring <string>		string to display when a score is not available (default: null)
            --csvout			generate comma-delimited CSV file (default: tab-delimited txt file)
            --argument <string>		comma-delimited strings as optional argument for each operation
            

 Function: automatically run a pipeline on a list of variants and summarize 
 their functional effects in a comma-delimited file, to be opened by Excel for 
 manual filtering
 
 Example: #recessive disease
          table_annovar.pl infile humandb/ -protocol refGene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,ljb_all -operation g,r,r,f,f,f,f -nastring NA
          table_annovar.pl infile humandb/ -buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all -operation g,r,r,f,f,f,f -nastring .
                  
 Version: $LastChangedDate: 2013-08-23 08:53:41 -0700 (Fri, 23 Aug 2013) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--protocol>

comma-delimited string specifying annotation protocol. These strings typically 
represent database names in ANNOVAR.

=item B<--operation>

comma-delimited string specifying type of operation. These strings can be g 
(gene), r (region) or f (filter).

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--(no)checkfile>

the program will check if all required database files exist before execution of annotation

=item B<--genericdbfile>

specify the genericdb files used in -dbtype generic. Note that multiple comma-
delimited file names can be supplied.

=item B<--gff3dbfile>

specify the GFF3 dbfiles files used in -dbtype gff3. Note that multiple comma-
delimited file names can be supplied.

=item B<--bedfile>

specify the GFF3 dbfiles files used in -dbtype bed. Note that multiple comma-
delimited file names can be supplied.

=item B<--vcfdbfile>

specify the VCF dbfiles files used in -dbtype vcf. Note that multiple comma-
delimited file names can be supplied.

=item B<--otherinfo>

print out otherinfo in the output file. "otherinfo" refers to all the 
infomration after fifth column in the input queryfile.

=item B<--onetranscript>

print out only one random transcript for exonic variants. By default, all 
transcripts are printed in the output.

=item B<--nastring>

string to display when a score is not available. By default, empty string is 
printed in the output file.

=item B<--csvout>

generate comma-delimited CSV file. By default, tab-delimited text file is generated.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.

The table_annovar.pl program is designed to replace summarize_annovar.pl in 
earlier version of ANNOVAR. Basically, it takes an input file, and run a series 
of annotations on the input file, and generate a tab-delimited output file, 
where each column represent a specific type of annotation. Therefore, the new 
table_annovar.pl allows better customization for users who want to annotate 
specific columns.

ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.


=cut