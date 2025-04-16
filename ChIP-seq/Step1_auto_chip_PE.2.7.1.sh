#!/bin/bash

# -i: bowtie2 index folder path
# -f: sequencing data folder path
# -o: output folder
# -e: effectiveGenomeSize for bamCoverage

while getopts 'i:f:o:e:' optname; do
	case $optname in
		i) genome_index="$OPTARG";;
		f) fastqfolder="$OPTARG";;
		o) output="$OPTARG";;
		e) effectiveGenomeSize="$OPTARG";;
		?) echo "ERROR, the input option is worse";;
	esac
done

echo "###########input_options############"
echo "genome_index_path:"$genome_index
echo "fastq_folder_path:"$fastqfolder
echo "output_path:"$output
echo "effective_Genomesize:"$effectiveGenomeSize

fastq_list=`ls $fastqfolder |grep 1.clean.fq.gz`

echo "############input_fastq_files##############"
echo $fastq_list

for filenames in $fastq_list
do
	# extract replicate names from files name. only *1.clean.fq.gz could be selected
	replicate_name=${filenames%_1.clean.fq.gz}
	echo "<----------$replicate_name file is analysizing ------------->"
	echo "<----------$replicate_name mapping_result ----------------->" >> $output"mapping_log.txt"
	time bowtie2 -p 16 -x $genome_index -1 $fastqfolder$replicate_name"_1.clean.fq.gz" -2 $fastqfolder$replicate_name"_2.clean.fq.gz" -S $output$replicate_name".unsorted.sam" 2>> $output"mapping_log.txt"
	samtools view -@ 8 -h -F 268 -q 10 -bS $output$replicate_name".unsorted.sam" -o $output$replicate_name".unsorted.bam"	
	samtools sort -@ 8 $output$replicate_name".unsorted.bam" -o $output$replicate_name".sorted.bam"
	samtools rmdup $output$replicate_name".sorted.bam" $output$replicate_name".sorted.filtered.bam"
	samtools index $output$replicate_name".sorted.filtered.bam"
	bamCoverage -p 8 --bam $output$replicate_name".sorted.filtered.bam" -o $output$replicate_name".sorted.filtered.bam.bw" --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize $effectiveGenomeSize
	if [ -f $output$replicate_name".sorted.filtered.bam.bw" ]
	then 
		rm $output$replicate_name".unsorted.sam" $output$replicate_name".unsorted.bam" $output$replicate_name".sorted.bam"
	else 
		echo "something wrong and this program didn't delete *.unsorted.sam and *.unsorted.bam files"
	fi
done
