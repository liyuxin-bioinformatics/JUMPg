#awk '{print $4"\t"$5"\t"$5"\t"$10"\t"$11}' $1 > $1.5col
perl /home/yli4/bin/customizedDB/MFM_SNV/extract_5cols.pl $1 > $1.5col
perl ~/annovar/annotate_variation.pl  -buildver hg19 -geneanno --seq_padding 30 -dbtype knowngene $1.5col /home/yli4/annovar/humandb/
perl /home/yli4/bin/customizedDB/MFM_SNV/seqpad2fas.pl $1.5col.seqpad mut $1.5col.seqpad.mut.fas
perl /home/yli4/bin/customizedDB/MFM_SNV/seqpad2fas.pl $1.5col.seqpad ref $1.5col.seqpad.ref.fas
#perl /home/yli4/bin/customizedDB/MFM_SNV/digest_seqpad_peptide.pl $1.5col.seqpad.mut.fas 30 $1.5col.seqpad.mut.digested.fas
#perl /home/yli4/bin/customizedDB/MFM_SNV/digest_seqpad_peptide.pl $1.5col.seqpad.ref.fas 30 $1.5col.seqpad.ref.digested.fas
#perl /home/yli4/bin/customizedDB/MFM_SNV/uniq_digested_peptides.pl $1.5col.seqpad.mut.digested.fas $1.5col.seqpad.mut.digested.uniq
perl /home/yli4/bin/customizedDB/shortID_convert.pl $1.5col.seqpad.mut.fas $1.5col.seqpad.mut.newIDs

