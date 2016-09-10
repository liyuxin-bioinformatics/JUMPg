~/bin/blast-2.2.18/bin/megablast -d ~/genomes/CDS_sequence_hg19_06272013 -i $1  -e 1 -D 3 -v 3 > $2_reads_vs_CDS.tab
~/bin/blast-2.2.18/bin/megablast -d ~/genomes/mRNA_sequence_hg19_06272013 -i $1  -e 1 -D 3 -v 3 > $2_reads_vs_mRNA.tab
~/bin/blast-2.2.18/bin/megablast -d ~/genomes/rnaseq_contaminants_07312013 -i $1  -e 1 -D 3 -v 3 > $2_reads_vs_contaminants.tab
~/bin/blast-2.2.18/bin/megablast -d ~/genomes/hg19_genome_raw -i $1  -e 1 -D 3 -v 3 > $2_reads_vs_WG.tab
~/bin/blast-2.2.18/bin/megablast -d ~/genomes/refFlat-sharp-alt_exon -i $1  -e 1 -D 3 -v 3 > $2_reads_vs_altExon.tab
~/bin/blast-2.2.18/bin/megablast -d ~/genomes/hg19_geneLocus_seq -i $1  -e 1 -D 3 -v 3 > $2_reads_vs_geneLocus.tab
#perl build_out2aligns.pl $1 out2read.hash $2_reads.alg $2_reads_vs_CDS.tab $2_reads_vs_mRNA.tab $2_reads_vs_contaminants.tab $2_reads_vs_altExon.tab $2_reads_vs_WG.tab

