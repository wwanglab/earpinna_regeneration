#!/bin/bash
  
bowtie2_index="/data/Genome/XXX/bowtie2_index/index"
fastq_list=` ls "/data/XXX/sequencing/00.CleanData/" | grep _1.clean.fq.gz`
echo $fastq_list

for i in $fastq_list
do
        echo ${i%_1.clean.fq}
	echo " samples were runing in bowtie2"
        file_name=${i%_1.clean.fq.gz}
	
        #call peak
        macs2 callpeak -t "/data/atac-seq_analysis/bam/"$file_name".f2.q30.bam" -f BAMPE -n $file_name -g mm -B --SPMR --keep-dup all --outdir "/data/atac-seq_analysis/callpeak"

        #tss enrichment
        #1. chrom size
        samtools view -H "/data/atac-seq_analysis/bam/"$file_name".f2.q30.bam"|perl -ne 'if(/SN:(\S+)\s+LN:(\d+)/){print "$1\t$2\n"}' > "/data/atac-seq_analysis/bam/"$file_name".chromSize"

        #bdg to bw
        bedSort "/data/atac-seq_analysis/callpeak/"$file_name"_treat_pileup.bdg" stdout|bedClip -truncate stdin "/data/atac-seq_analysis/bam/"$file_name".chromSize" stdout|perl -ane 'print if($F[1]<$F[2])' > tmp.bdg
        bedGraphToBigWig tmp.bdg "/data/atac-seq_analysis/bam/"$file_name".chromSize" "/data/atac-seq_analysis/bw/"$file_name"_treat_pileup.bw"


done
