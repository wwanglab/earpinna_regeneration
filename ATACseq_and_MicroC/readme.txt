######### ATAC-seq
ATAC_seq_alignment.sh: bowtie2 alignment
ATAC_seq_callpeak.sh: macs2 call peak


######### micro-C
micro_C.sh: micro-C data analysis including bwa alignment, generating valid pairs file and HiC contact maps, learned from https://micro-c.readthedocs.io/en/latest/
get_qc.py: culculate the quality of micro-C library, downloaded from https://github.com/dovetail-genomics/Micro-C.git 
