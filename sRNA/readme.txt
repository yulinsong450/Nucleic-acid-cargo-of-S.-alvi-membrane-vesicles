The custom Bash scripts were used to align small RNA reads to reference templates and sort alignments.

The scripts were run in a conda environment (python==3.12.8) with the following packages installed:
cutadapt v4.9
bowtie v1.3.1
samtools v1.21

To reproduce the analysis, before running the scripts: 
1) change a suitable directory for "${mainfolder}",
2) copy the raw .fastq.gz files of sRNA sequencing (NCBI BioProject: PRJNA1404345) into "${mainfolder}YS_G_raw/",
3) copy the reference sequence file B2_GFPw_pDS.fasta into "${mainfolder}bowtie1_ref/".

Then, run "bash sRNA_Bash_scripts.sh"

Bash outputs were used to plot in R. Reads that mapped to GFP (G_mapped_bam_GFP_sorted/) were further visualized via IGV.