Raw data of Illumina whole genome sequencing are available on NCBI (BioProject: PRJNA1404345).

Short reads were aligned to reference sequences via breseq.
E.g.:
breseq -r wkB2.fasta -r pDS-GFP.gb Cells_G1.fastq.gz -o ./Cells_G1

breseq coverages ("fit mean") that mapped to chromosome (CP007446.1) and pDS-GFP were used to calculate the plasmid-chromosome ratio.