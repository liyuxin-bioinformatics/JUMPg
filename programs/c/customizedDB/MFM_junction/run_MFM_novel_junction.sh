perl /home/yli4/bin/customizedDB/MFM_junction/select_novel_junctions.pl $1.txt 3 3 $1_novel.txt
perl /home/yli4/bin/customizedDB/MFM_junction/extractFlankingGenomicSeq.pl /home/yli4/genomes/hg19_genome_raw.fa $1_novel.txt 66 $1_novel.fas
perl /home/yli4/bin/customizedDB/MFM_junction/junction_seq_translation.pl $1_novel.fas /home/yli4/annotations/knownGenes_geneSymble_strand_120613.txt $1_novel_AA.fas
perl /home/yli4/bin/customizedDB/MFM_junction/check_AA_sequence.pl $1_novel_AA.fas $1_novel_AA_central.fas
perl /home/yli4/bin/customizedDB/MFM_junction/shortID_convert.pl $1_novel_AA_central.fas $1_novel_AA_central_newIDs
