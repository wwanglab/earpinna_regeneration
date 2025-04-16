#alignment
bwa mem -5SP -T0 -t16 /data/Genome/rabbit_UM_NZW_1.0/UM_NZW_1.0_genomic.fna /data/micro-C_analysis/orig.data_610G/rb_D10_1.clean.fq.gz /data/micro-C_analysis/orig.data_610G/rb_D10_2.clean.fq.gz -o Rb_bulk_micro-C.sam

#find ligation junctions
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path /data/Genome/rabbit_UM_NZW_1.0/genometable.UM.NZW.txt Rb_bulk_micro-C.sam > Rb_bulk_micro-C.parsed.pairsam

#Sort pairsam file
pairtools sort --nproc 8 --tmpdir /data/micro-C_analysis/analysis/temp/ Rb_bulk_micro-C.parsed.pairsam > Rb_bulk_micro-C.parsed.sorted.pairsam

#Remove PCR duplicates
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output Rb_bulk_micro-C.parsed.dedup.sorted.pairsam Rb_bulk_micro-C.parsed.sorted.pairsam

#Generate .pairs and bam files
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs Rb_bulk_micro-C.mapped.pairs --output-sam Rb_bulk_micro-C.unsorted.bam Rb_bulk_micro-C.parsed.dedup.sorted.pairsam

#Generate the final bam file
samtools sort -@16 -T /data/micro-C_analysis/analysis/temp/temp.bam -o Rb_bulk_micro-C.mapped.PT.bam Rb_bulk_micro-C.unsorted.bam
samtools index Rb_bulk_micro-C.mapped.PT.bam

python3 /data/Micro-C/get_qc.py -p stats.txt > QC.txt

preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output out.preseq Rb_bulk_micro-C.mapped.PT.bam

java -Xmx48000m -Djava.awt.headless=true -jar /data/juicertools.jar pre --threads 16 Rb_bulk_micro-C.mapped.pairs Rb_bulk_micro-C.contact_map.hic /data/Genome/rabbit_UM_NZW_1.0/genometable.UM.NZW.txt

