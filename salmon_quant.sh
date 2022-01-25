#!/bin/bash
dir="/Volumes/SEAGATE/Post_viva_rnaseq/Trimmed_files"
files=`ls /Volumes/SEAGATE/Post_viva_rnaseq/Trimmed_files | grep -E "JBE*|7.*"`
for i in $files
do
salmon quant -i gencode.vM23.transcripts_index_17_04_2020 -l ISR -1 $dir/$i/*val_1.fq.gz -2 $dir/$i/*val_2.fq.gz --useVBOpt --numBootstraps 30 --seqBias --gcBias --posBias -o $i --threads 2
done
