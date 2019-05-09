----------------------
Introduction
----------------------

JUMPg is a proteogenomics software pipeline for analyzing large mass spectrometry (MS) and functional genomics datasets. The pipeline includes customized database building, tag-based database search, peptide-spectrum match filtering, and data visualization. The program is written in perl and is designed for high performance parallel computing systems. 
 
For more information, please read: Li et al., JUMPg: an integrative proteogenomics pipeline identifying unannotated proteins in human brain and cancer cells. J. Proteome Res., 15 (2016), pp. 2309–2320

----------------------
Contents of this file
----------------------

* Introduction
* Release notes
* Software Requirements
* Hardware Requirements
* Installation
* Command Line Arguments
* Maintainers

----------------------
Release notes
----------------------

Version 2.3.6:
* Fixed a bug that cause program to crash when the program is used in Ubuntu system.
* Use the latest JUMP search engine (thanks to Xusheng and Ji-Hoon).
* During the 2nd stage search, some PSMs may be missed during database search because of inaccurate file calculation. This error is fixed.
* Previous version failed to generate the BED file for mutation peptides. This error is fixed.
* The PSM filtering function previously relied on a local version of clone.pm, but may lead to crash when a local clone.pm file was found with a differnt version. Now the local dependence was removed and this error should be fixed.

----------------------
Software Requirements
---------------------- 

The program is written by a combination of Perl and R. It should run on every system with Perl5 and R 3.1.0. The minimum required Perl version should be Perl 5.8 or better.

Perl modules are needed:
- Parallel::ForkManager
- Statistics::Distributions
- Clone
- Excel/Writer/XLSX.pm (optional but recommended)

To install perl modules, please refer to:
http://www.cpan.org/modules/INSTALL.html

R module:
- MASS

----------------------
Hardware Requirements
---------------------- 

Starting from JUMPg_v2.3.1, the program can be run on either high performance computing systme or a single server. 

The program has been successfully tested on the following system:

+ Cluster mode (key parameters: 'cluster = 1' & 'Job_Management_System = SGE'):
  - Batch-queuing system: SGE, version 6.1u5, qsub and qstat
  - 128 GB memory and 64 cores on each node

+ Single server mode (key parameters: 'cluster = 0' & 'processors_used = 8'): 
  - 32 GB memory
  - 2 GHz CPU processors with 8 cores
 
----------------------
Installation
---------------------- 

After download the source code, you can put it in any working directory (e.g. /home/usr/JUMPg). 
IMPORTANT: The folder containing all the source code needs to be accessible by each node.


INSTALLATION:

__step 0__: unzip the source code and test data packages in the current directory by running the following commands:
 - unzip JUMPg_v2.3.1.zip (already unzipped if downloaded from github)
 - unzip exampleData_v6.2.zip

__Step 1__: set up module and program paths by running the following commands:
 - cd programs
 - perl path_setup.pl
 - cd ..
 
__Step 2__: download AnnoVar annotation files by running the following commands (take human hg19 as an example):
 - mkdir annotations
 - cd annotations
 - perl ../programs/c/customizedDB/annovar/annotate_variation.pl -downdb -buildver hg19 -downdb knownGene human/
 - cd ..
 
__Step 3__: prepare reference genome (FASTA format)  # this file is required for splice junction and RNA 6FT analysis
 - For model organisms, the reference genome should be available through UCSC genome browser database (http://genome.ucsc.edu/). 
 
   Example command: 
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
   Concatenate all the fasta files using the following command
   cat *.fa > hg19_genome_raw.fa

__Step 4__: (Optional but recommended) update uniProt-UCSC ID convertion table based on pairwise sequence alignment. 
> For human, an updated uniProt-UCSC ID convertion table (named 'ID_convert_uniProt2UCSChg19.txt') is available in the exampleData_v6.2/ folder. This table is used to correct the UCSC downloaded table, using the following command:
  - perl tools/updateIdTab.pl annotations/human/hg19_kgXref.txt exampleData_v6.2/resources/ID_convert_uniProt2UCSChg19.txt


----------------------
Input File Preparation
----------------------

(You can skip this section if you only want to run through the example.)

Besides MS data, the program takes up to three types of genomics/NGS data as inputs for the customized peptide database building. 

1) For building the mutation peptide database: an AnnoVar input format file that contains the genomic locations and allele (ref/alt) information file is required.
2) For building the junction peptide database: a splice junction file reported by STAR (or with equivalent format) is required.
3) For building the RNA 6FT database: RNA nucleotide sequences in FASTQ format is required (can be either de novo assemebled transcripts or reference mRNA sequences).

Users can find example input files in exampleData_v6.2.tar.gz (within the rna_database/ folder).

If users process the data from RNA-seq raw reads, to assist the input file preparation process, a wrapper (called 'RNAseq_preprocess.pl' located in the tools/ folder) is provided that contains RNAseq alignment (by STAR), mutation detection (by GATK), and novel splice junction filtering. 

To run the wrapper, the following 3rd party software is required:
1) STAR: https://github.com/alexdobin/STAR/releases
Note that the reference genome should be index by the STAR 'genomeGenerate' module 
2) Picard: http://broadinstitute.github.io/picard/
3) GATK: https://www.broadinstitute.org/gatk/
Note that GATK requires the reference genome to be sorted by kayrotypic order
4) Samtools: http://samtools.sourceforge.net/
Only used for indexing BAM files

Users can check each requirement by filling the parameter file of the wrapper ('RNAseq_preprocess.params' located in the tools/ folder).

Once finished with the steps above, the wrapper can be executed using the following command:

perl tools/RNAseq_preprocess_v1.0.0.pl RNAseq_preprocess.params <test_PE1.fq> <test_PE2.fq>
<test_PE1.fq> <test_PE2.fq> Pair 1 and 2 of RNAseq FASTAQ files

The wrapper will output:
1) Mutation file in AnnoVar input format: filtered_output.annovar
2) Novel splice junction file in STAR format: SJ.out.filtered2uniq.tab

Users can take these two files as input for the JUMPg program.

----------------------
Run the example
----------------------

Since the program supports multistage database search analysis, we have designed a two-stage test dataset. For the 1st stage, the MS/MS spectra are searched against a peptide database pooled from uniProt, mutations and splice junctions; the matching results are filtered to ~1% false discovery rate. For the 2nd stage, the remaining high quality spectra are searched against the 6FT database.

1) 1st stage analysis:

-Step 1: cp 'jump_g_v2.2.stage1.params' (included in exampleData_v6.2/parameters folder) to your working directory and edit the following parameters:

The following assumes that your working directory is at "/home/usr/JUMPg"

  a) input_ref_proteins = /home/usr/JUMPg/exampleData_v6.2/uniProt/reference_uniProt_human.fasta
  b) input_mutation = /home/usr/JUMPg/exampleData_v6.2/rna_database/mutations.txt
  c) input_junction = /home/usr/JUMPg/exampleData_v6.2/rna_database/junctions.txt
  d) annovar_annotation_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_knownGene.txt
  e) gene_ID_convert_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_kgXref.txt
  f) reference_genome = /home/usr/genomes/hg19_genome_raw.fa

-Step 2: cp spectra file 'exampleData_v6.2/spectra/MS_test.mzXML' (included in exampleData_v6.2.tar.gz) to your working diredirectory and run the program using the following command:

perl /home/usr/JUMPg_v2.3.1/programs/JUMPg_v2.3.1.pl jump_g_v2.2.stage1.params MS_test.mzXML

-Output: the program will generate a folder named 'gnm_round1_test1'. The final results are all included in the 'publications' folder that includes 6 files:
  a) identified_peptide_table.txt: a text table, showing the identified peptide sequences, PSM counts, tags and scores of the best PSM, and database entries.
  b) mutation_peptides.txt: a text table showing the identified mutation peptides
  c) mutation_peptides.bed: mutation peptides with genomic locations in BED format, which can be co-displayed with other genomic information on the UCSC genome browser
  d) junction_peptides.txt: a text table showing the identified novel junction peptides
  e) junction_peptides.bed: novel junction peptides with genomic locations in BED format for visualization on the UCSC genome browser
  f) reference_peptides.bed: reference peptides with genomic locations in BED format for visualization on the UCSC genome browser

For multistage analysis, the program also generates the unmatched high quality MS/MS spectra, of which the path is recorded in this file: gnm_round1_test1/multistage/qc_MSMS_input.txt. This file will be taken as input for 2nd stage analysis.

2) 2nd stage analysis:

-Step 1: cp 'jump_g_v2.2.stage2.params' (included in exampleData_v6.2/parameters) to your working directory and edit the following parameters:
  a) input_transcript_6FT = /home/usr/JUMPg/exampleData_v6.2/rna_database/assembled_transcripts.fasta
  b) annovar_annotation_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_knownGene.txt
  c) gene_ID_convert_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_kgXref.txt
  d) reference_genome = /home/usr/genomes/hg19_genome_raw.fa
  f) reference_protein = /home/usr/JUMPg/exampleData_v6.2/uniProt/reference_uniProt_human.fasta

-Step 2: copy qc_MSMS_input.txt (that records the path to unmatched high quality MS/MS spectra) to current directory:

cp gnm_round1_test1/multistage/qc_MSMS_input.txt .

-Step 3: run the program by the command:
perl /home/usr/JUMPg_v2.3.1/programs/JUMPg_v2.3.1.pl jump_g_v2.2.stage2.params qc_MSMS_input.txt

-Output: similar to the 1st stage result, the program will generate a folder named 'gnm_round2_test1' containing results in its 'publications' folder.


----------------------
Maintainers
----------------------

* To submit bug reports and feature suggestions, please contact:

  Yuxin Li (yuxin.li@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

