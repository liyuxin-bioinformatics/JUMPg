perl /home/yli4/bin/customizedDB/RDT/readQC_IDshort.pl $1.fq $2 $1-QC.fas
perl /home/yli4/bin/customizedDB/RDT/fas_tr6.pl $1-QC.fas $1-QC
perl /home/yli4/bin/customizedDB/RDT/tryptic_digest.pl $1-QC_peptide.fas KR P $1-QC_peptide_digested.fas
perl /home/yli4/bin/customizedDB/RDT/uniq_digested_peptides.pl $1-QC_peptide_digested.fas $1-QC_peptide_digested_uniq
perl /home/yli4/bin/customizedDB/RDT/pepToReads.pl $1-QC_peptide_digested_uniq.readCounts $1-QC_peptide_digested_uniq.pepToReads
perl /home/yli4/bin/customizedDB/MFM_junction/shortID_convert.pl $1-QC_peptide_digested_uniq.fas $1-QC_peptide_digested_uniq.newIDs

wc -l $1.fq
wc -l $1-QC.fas
wc -l $1-QC_peptide.fas
wc -l $1-QC_peptide_digested.fas
wc -l $1-QC_peptide_digested_uniq.fas
wc -l $1-QC_peptide_digested_uniq.pepToReads

