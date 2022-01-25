#!/bin/bash
file_dir=/Volumes/SEAGATE/Post_viva_rnaseq
curr_dir=/Volumes/SEAGATE/Post_viva_rnaseq/Trimmed_files
cd $curr_dir
ls /Volumes/SEAGATE/Post_viva_rnaseq | sort | sed s/...fastq.gz/""/g | uniq > files.txt
for file in `cat files.txt`
do
SAMPLENAME=$(echo $file | awk -F "_" '{print $3}')
echo $SAMPLENAME
mkdir -p $SAMPLENAME
cd $SAMPLENAME
trim_galore -q 30 --fastqc --cores 8 --stringency 1 --max_n 1 --length 25 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" -a2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" --paired ${file_dir}/${file}_1.fastq.gz ${file_dir}/${file}_2.fastq.gz
cd $curr_dir
done
