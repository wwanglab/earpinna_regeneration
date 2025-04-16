#!/bin/bash

# -i: star index folder path
# -d: sequencing data folder path
# -g: path of gtf file

# make tree option for input star index path and sequencing data path
while getopts 'i:d:g:' optname; do
	case $optname in
		i) index_path="$OPTARG";;
		d) data_path="$OPTARG";;
		g) gtf_path="$OPTARG";;
		?) echo "ERROR, the input option is worse";;
	esac
done
echo $index_path
echo $data_path
echo $gtf_path
# make a new folder termed "output" for output files
if [ ! -d "./output/" ];then
	mkdir ./output
else
	echo "Warning:output folder already existed"
fi

# make a list of all fastq files
# the file name must be *.1.clean.fq.gz format, if you want to change other format, you should change this line
fastq_list=`ls $data_path |grep 1.clean.fq.gz`
echo $fastq_list


for filenames in $fastq_list
do
	# extract replicate names from files name. only *.1.clean.fq.gz could be selected
	replicate_name=${filenames%_1.clean.fq.gz}
	echo "<----------$replicate_name is alignmenting with STAR------------->"
	STAR --runMode alignReads \
		--outSAMtype BAM SortedByCoordinate \
		--readFilesCommand gunzip -c \
		--genomeDir $index_path \
		--quantMode GeneCounts \
		--sjdbGTFfile $gtf_path \
		--outFileNamePrefix ./output/$replicate_name \
		--readFilesIn $data_path$replicate_name"_1.clean.fq.gz" $data_path$replicate_name"_2.clean.fq.gz" \
		--runThreadN 8
done
