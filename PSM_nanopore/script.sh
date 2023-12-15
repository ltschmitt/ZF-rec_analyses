mkdir basecalled

guppy_basecaller --input_path fast5 --save_path basecalled -c dna_r9.4.1_450bps_hac.cfg --device cuda:0 --min_qscore 10

cat basecalled/pass/*.fastq | gzip -c > reads.fastq.gz

zcat reads.fastq.gz | parallel --pipe --max-lines 40 --cat 'sed -i "3~4d;4~4d;1~4s/^.\([^ ]*\).*/>\1/" {} ; exonerate --ryo "{%ti\t%Pl\t%Pqb\t%Pqs\t%Pts\n}" --verbose 0 --showalignment no --showvulgar no -E -S no --percent 80 --geneticcode FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG --model protein2dna:bestfit vika_aa.fa {}' ::: basecalled/pass/*.fastq | gzip -c > reads_vs_vika_aa.exonerate.txt.gz
