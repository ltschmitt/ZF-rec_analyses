#!/bin/bash
mkdir aligned
mkdir analysis

# filter
filtlong --min_length 2000 --min_mean_q 10  brec1_znf.fastq > filtered_brec1_znf.fastq

# align
minimap2 -a -N 50 -k 13 -t 11 -x map-ont references/references.fa fastq/filtered_brec1_znf.fastq | samtools sort | samtools view -Sb > aligned/n50_brec1_znf.bam
samtools index aligned/n50_brec1_znf.bam

echo -e "ID\tFlag\tMatch\tStart\tLength\tEditDistance" > analysis/filtered_n50_brec1_znf.tsv
samtools view -h -F 4 -F 256 -F 2048 -@ 4 -L references/rec.bed aligned/n50_brec1_znf.bam | samtools view -@ 4 -L references/ts.bed - | awk -F '\t' '{ OFS="\t" } { print $1,$2,$3,$4,length($10),$12}' >> analysis/filtered_n50_brec1_znf.tsv

