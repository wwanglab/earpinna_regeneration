#/bin/bash

# -b: folders containing sorted bamfile
# -i: input .sorted.filtered.bam file path
# -o: output folder
# -g: -g parameters for macs2, for mice is "mm", for zebrafish is "1.4e9", for human is "hs"

while getopts 'b:i:o:g:' optname; do
        case $optname in
                b) bamfile_folders="$OPTARG";;
		i) input_file="$OPTARG";;
		o) output_folders="$OPTARG";;
		g) genome_size="$OPTARG";;
                ?) echo "ERROR, the input option is worse";;
        esac
done

file_list=`ls $bamfile_folders | grep .sorted.filtered.bam$`
#input_file=`ls $input_folders | grep .sorted.filtered.bam$`

echo $file_list
echo $input_file

for filenames in $file_list
do
	prefix=${filenames%.sorted.filtered.bam}
	echo "<--------bamfile is:$bamfile_folders$filenames ------------->"
	echo "<--------inputfile is:$input_file ----------->"

	macs2 callpeak -t $bamfile_folders$filenames -c $input_file -g $genome_size -n $prefix".narrow" -f BAMPE -B -q 0.01 --outdir $output_folders 2> $output_folders$prefix".narrow.log"
#	macs2 callpeak -t $bamfile_folders$filenames -c $input_file -g $genome_size -n $prefix".narrow" -f BAMPE -B -q 0.01 --keep-dup all --outdir $output_folders 2> $output_folders$prefix".narrow.log"
	
#	macs2 callpeak -t $bamfile_folders$filenames -c $input_file -g $genome_size -n $prefix".broad" --broad --broad-cutoff 0.01 --outdir $output_folders 2> $output_folders$prefix".broad.log"
#	macs2 callpeak -t $bamfile_folders$filenames -c $input_file -g $genome_size -n $prefix".nomodel" --broad --broad-cutoff 0.01 --outdir $output_folders --nomodel --shift 72 2> $output_folders$prefix".nomodel.log"

done
