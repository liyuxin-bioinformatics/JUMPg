#!/usr/bin/perl

use strict;
use warnings;

if(scalar(@ARGV)<3)
{
        help();
}

my $params=shift @ARGV;
my $fq1=shift @ARGV;
my $fq2=shift @ARGV;

my (%parahash);
parse_params($params,\%parahash);

run_commands(\%parahash,$fq1,$fq2);

#--------------------------------------------------------------------------------------------------------------------------------------
sub run_commands {
	my ($parahash,$fq1,$fq2)=@_;

	# STAR
	system(qq($parahash{STAR} --runThreadN $parahash{runThreadN} --genomeDir $parahash{genome_dir} --readFilesIn $fq1 $fq2 --outFileNamePrefix $parahash{outFileNamePrefix} --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate));

	# Picard
	# AddOrReplaceReadGroups
	system(qq(java  -jar $parahash{Picard} AddOrReplaceReadGroups I=$parahash{outFileNamePrefix}Aligned.sortedByCoord.out.bam O=rg_added_sorted.bam SO=coordinate RGID=$parahash{outFileNamePrefix} RGPL=Illumina RGLB=v2 RGPU=HiSeq2000 RGSM=$parahash{outFileNamePrefix}));
	# MarkDuplicates
	system(qq(java  -jar $parahash{Picard} MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics));
	# ReorderSam (by kayrotypic order)
	system(qq(java  -jar $parahash{Picard} ReorderSam I= dedupped.bam O= kayrotypic.bam REFERENCE=$parahash{kayrotypic_ordered_genome}));

	# samtools: index bam
	system(qq(samtools index kayrotypic.bam));

	# GATK
	# SplitNCigarReads
	system(qq(java -Djava.io.tmpdir=tmp -jar $parahash{GATK} -T SplitNCigarReads -R $parahash{kayrotypic_ordered_genome} -I kayrotypic.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS));
	# HaplotypeCaller
	system(qq(java -jar $parahash{GATK} -T HaplotypeCaller -R $parahash{kayrotypic_ordered_genome} -I split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o output.vcf));
	# VariantFiltration
	system(qq(java -jar $parahash{GATK} -T VariantFiltration -R $parahash{kayrotypic_ordered_genome} -V output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o filtered_output.vcf));

	# vcf2annovar.pl
	vcf2annovar('filtered_output.vcf','filtered_output.annovar');

	# filterSTARjunc.pl
	filterSTARjunc("$parahash{outFileNamePrefix}SJ.out.tab",$parahash{junction_min_uniq_reads},"$parahash{outFileNamePrefix}SJ.out.filtered2uniq.tab");
}

sub filterSTARjunc {
	my $star_junc=shift @_;
	my $minUniReads=shift @_;
	my $output=shift @_;

	open(IN,$star_junc) || die "Cannot open STAR junction file $star_junc\n";
	open(OUT,">$output");

	while(<IN>)
	{
        	chomp;
	        my ($chr,$start,$end,$strand,$motifCode,$annotated,$uniReads,$multiReads,$maxOverhang)=split /\t/,$_;

        	next if ($annotated);   # only novel junction
	        next if ($uniReads<$minUniReads); # min unique mapped reads
        	next if ($motifCode ==0 ||      # non-canonical motif
                	$strand == 1 && $motifCode % 2 == 0 || # + strand, '-' motif
	                $strand == 2 && $motifCode % 2 == 1);  # - strand, '+' motif

        	print OUT $_,"\n";
	}
	close IN;
	close OUT;
}

sub vcf2annovar {
	my ($vcf,$output)=@_;
	open(IN,$vcf) || die "Cannot open vcf file $vcf\n";
	open(OUT,">$output");

	while(<IN>)
	{
	        chomp;
        	s/^\s*//;
	        next if (/^#/);

        	my ($chr,$pos,$id,$ref,$alts,$score,$filter,$infor)=split /\t/,$_;
	        next unless ($filter eq 'PASS'); # only use 'filtered' mutations

        	my @altss=split /\,/,$alts;
	        foreach my $alt (@altss)
        	{
                	my ($start,$end);
	                if ( length($ref)==1 and length($alt)==1 ){
        	                $start=$end=$pos;
                	}
	                elsif ( length($ref)==1 and length($alt)>1 ){
        	                $start=$end=$pos+1;
                	        $ref='-';
                        	$alt=shiftHeadLetter($alt);
	                }
        	        elsif ( length($ref)>1 and length($alt)==1 ){
                	        $start=$pos+1;
                        	$end=$pos+length($ref)-1;
	                        $alt='-';
                	        $ref=shiftHeadLetter($ref); #
        	        }
	                else { next; }
        	        	print OUT "$chr\t$start\t$end\t$ref\t$alt\n";
        		}
	}
	close IN;
	close OUT
}

sub shiftHeadLetter
{
        my ($s)=@_;
        $s=reverse($s);
        chop($s);
        return reverse($s);
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
#                        $key =~ s/^input_//;
                       $$parahash{$key}=$value;
                }
        }
        close IN;
}

sub help {
	print  "\nUsage: perl RNAseq_preprocess_v1.0.0.pl RNAseq_preprocess.params RNAseq_reads_PE1.fq RNAseq_reads_PE2.fq\n\n";
	exit;
}
