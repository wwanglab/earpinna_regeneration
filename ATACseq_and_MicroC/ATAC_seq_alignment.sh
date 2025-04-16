#!/bin/bash
  
bowtie2_index="/data/Genome/XXX/bowtie2_index/index"
fastq_list=` ls "/data/XXX/sequencing/00.CleanData/" | grep _1.clean.fq.gz`
echo $fastq_list

for i in $fastq_list
do
	echo ${i%_1.clean.fq}
	echo "samples were runing in bowtie2"
	file_name=${i%_1.clean.fq.gz}
	bowtie2 -p 16 -q -I 10 -X 1000 --dovetail --no-unal --very-sensitive-local --no-mixed --no-discordant -x $bowtie2_index -1 "/data/XXX/sequencing/00.CleanData/"$file_name"_1.clean.fq.gz" -2 "/data/XXX/sequencing/00.CleanData/"$file_name"_2.clean.fq.gz" 2>"/data/atac-seq_analysis/alignment/"$file_name".alignment.summary" | samtools view -F 4 -u - | samtools sort -@ 6 -o "/data/atac-seq_analysis/bam/"$file_name".bam"

	samtools view -b -f 2 -q 30 -o "/data/atac-seq_analysis/bam/"$file_name".f2.q30.bam" "/data/atac-seq_analysis/bam/"$file_name".bam"

	#preseq
	preseq lc_extrap -P -B -D -v -o "/data/atac-seq_analysis/preseq/"$file_name".dat" "/data/atac-seq_analysis/bam/"$file_name".f2.q30.bam" 2>"/data/atac-seq_analysis/preseq/"$file_name".log" 

	#picard
	picard MarkDuplicates PG=null VERBOSITY=ERROR QUIET=true CREATE_INDEX=false REMOVE_DUPLICATES=ture INPUT="/data/atac-seq_analysis/bam/"$file_name".f2.q30.bam" OUTPUT="/data/atac-seq_analysis/bam/"$file_name".f2.q30.dedup.bam" M="/data/atac-seq_analysis/bam/"$file_name".pe.markduplicates.log"

done
