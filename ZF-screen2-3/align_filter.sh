#!/bin/bash
mkdir aligned
mkdir analysis

# filter
filtlong --min_length 2000 --min_mean_q 10  fastq/reads_zfscreen2.fastq > fastq/filtered_reads_zfscreen2.fastq
filtlong --min_length 2000 --min_mean_q 10  fastq/reads_zfscreen3.fastq > fastq/filtered_reads_zfscreen3.fastq

# align
minimap2 --secondary no -t 11 -ax map-ont references/zfscreen2_references.fa fastq/filtered_reads_zfscreen2.fastq | samtools sort | samtools view -Sb > aligned/zfscreen2.bam
samtools index aligned/zfscreen2.bam
minimap2 --secondary no -t 11 -ax map-ont references/zfscreen3_references.fa fastq/filtered_reads_zfscreen3.fastq | samtools sort | samtools view -Sb > aligned/zfscreen3.bam
samtools index aligned/zfscreen3.bam

# extract
echo -e "ID\tFlag\tMatch\tStart\tLength\tEditDistance" > analysis/filtered_zfscreen2.tsv
samtools view -h -F 4 -F 256 -F 2048 -@ 4 -L references/rec_zfscreen2.bed aligned/zfscreen2.bam | samtools view -@ 4 -L references/ts_zfscreen2.bed - | awk -F '\t' '{ OFS="\t" } { print $1,$2,$3,$4,length($10),$12}' >> analysis/filtered_zfscreen1.tsv

echo -e "ID\tFlag\tMatch\tStart\tLength\tEditDistance" > analysis/filtered_zfscreen3.tsv
samtools view -h -F 4 -F 256 -F 2048 -@ 4 -L references/rec_zfscreen3.bed aligned/zfscreen3.bam | samtools view -@ 4 -L references/ts_zfscreen3.bed - | awk -F '\t' '{ OFS="\t" } { print $1,$2,$3,$4,length($10),$12}' >> analysis/filtered_zfscreen3.tsv

