#!/usr/bin/env bash
# Author: Yulin Song <yulinsong450@gmail.com>
# 2026-01
# Function: align small RNA reads to reference templates and sort the alignments
# IMPORTANT: to reproduce the analysis, before running the scripts, make sure that the raw .fastq.gz files of sRNA sequencing have been copied into "${mainfolder}YS_G_raw/" and the reference sequence file B2_GFPw_pDS.fasta has been copies into "${mainfolder}bowtie1_ref/".

set -e

mainfolder="/stor/work/Ochman/yulin/MV_sRNA_GFP/" # Change to any suitable working directory.


#############################
### 1. cutadapt trimming
#############################

mkdir "${mainfolder}YS_G_trimmed/"
cd "${mainfolder}YS_G_raw/"

# Loop through .fastq.gz files and write a trimming report
for file in *.fastq.gz; do 
printf "\n\n${file}\n" >> ../YS_G_trimming_report.txt && \
cutadapt --discard-untrimmed --report=minimal -j 20 -q 15 -u 1 -a TGGAATTCTCGGGTGCCAAGG -m 15 -o ../YS_G_trimmed/"${file}" "${file}" >> ../YS_G_trimming_report.txt
done
echo "finish cutadapt trimming"



#############################
### 2. bowtie aligning
#############################

bowtiefolder="${mainfolder}YS_G_trimmed_bowtie1_v0/"
mkdir -p "${bowtiefolder}G_unmapped/" "${bowtiefolder}G_mapped/" "${bowtiefolder}G_mapped_bam" "${bowtiefolder}G_mapped_bam_sorted/"

log="${bowtiefolder}YS_G_trimmed_bowtie1_v0_log.txt"
printf "YS_G_trimmed_bowtie1_v0_align_to_B2_GFPw_pDS\n\n" > "$log"

# Build a bowtie index for the reference sequences
cd "${mainfolder}bowtie1_ref/"
bowtie-build B2_GFPw_pDS.fasta B2_GFPw_pDS

# Loop through trimmed files
for file in ${mainfolder}YS_G_trimmed/*.fastq; do
filename=$(basename "$file" .fastq) && \

# Bowtie
printf "\n##########\n\n${filename}\n" >> "${log}" && \
bowtie -v 0 --no-unal --threads 30 -q -S --un "${bowtiefolder}G_unmapped/${filename}_unmapped.fastq" -x B2_GFPw_pDS "${file}" > "${bowtiefolder}G_mapped/${filename}_mapped.sam" 2>> "${log}" && \

# Convert SAM to BAM
samtools view -bS "${bowtiefolder}G_mapped/${filename}_mapped.sam" > "${bowtiefolder}G_mapped_bam/${filename}_mapped.bam" && \

# Sort BAM
samtools sort "${bowtiefolder}G_mapped_bam/${filename}_mapped.bam" -o "${bowtiefolder}G_mapped_bam_sorted/${filename}_mapped.bam" && \

# Index sorted BAM
samtools index "${bowtiefolder}G_mapped_bam_sorted/${filename}_mapped.bam"

done

echo "finish bowtie aligning"



#############################
### 3. count reads
#############################

cd "${bowtiefolder}"

echo "sample,template,count" > "YS_G_count_reads.csv"

# Loop through .sam files
for sam_file in G_mapped/*.sam; do

    # Extract sample name from filename
    filename=$(basename "$sam_file" .sam)

    # Loop over the three reference templates
    for template in NZ_CP007446.1 GFP_wPT pDS-noGFPw; do

        # Count reads aligned to different templates
        count=$(awk -v ref="$template" '
            $1 !~ /^@/ && $3 == ref {c++}
            END {print c+0}
        ' "$sam_file")

        # Write to CSV
        echo "$filename,$template,$count" >> "YS_G_count_reads.csv"

    done
done

echo "finish counting reads"



#############################
### 4. count read lengths
#############################

echo "Sample,Category,Length,Count" > YS_G_aligned_read_length_distribution.csv

# Loop through each sample
for sample in ${bowtiefolder}G_mapped/*_mapped.sam; do
  sample_name=$(basename "$sample" .sam)
  echo "Processing $sample_name"
  
  # Loop through each category
  for category in NZ_CP007446.1 GFP_wPT pDS-noGFPw; do
    
    # Extract read lengths and their counts
    awk -v sample="$sample_name" -v category="$category" '
      $1 !~ /^@/ && $3 == category {lengths[length($10)]++}
      END {
        for (len in lengths) {
          print sample "," category "," len "," lengths[len]
        }
      }
    ' "$sample" >> YS_G_aligned_read_length_distribution.csv
  done
done

echo "finish counting read lengths"



#############################
### 5. extract GFP reads
#############################

mkdir "${bowtiefolder}G_mapped_bam_GFP/" "${bowtiefolder}G_mapped_bam_GFP_sorted/"

# Loop
for file in ${bowtiefolder}G_mapped_bam_sorted/*.bam; do
filename=$(basename "${file}" .bam) && \

# Extract reads mapped to GFP_wPT
samtools view -b "${file}" GFP_wPT > "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP.bam" && \

# Split sense and antisense reads
samtools view -b -F 0x10 "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP.bam" > "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP_sense.bam" && \
samtools view -b -f 0x10 "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP.bam" > "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP_antisense.bam" && \

# Sort BAM
samtools sort "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP.bam" -o "${bowtiefolder}G_mapped_bam_GFP_sorted/${filename}_GFP.bam" && \
samtools sort "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP_sense.bam" -o "${bowtiefolder}G_mapped_bam_GFP_sorted/${filename}_GFP_sense.bam" && \
samtools sort "${bowtiefolder}G_mapped_bam_GFP/${filename}_GFP_antisense.bam" -o "${bowtiefolder}G_mapped_bam_GFP_sorted/${filename}_GFP_antisense.bam" && \

# Index sorted BAM
samtools index "${bowtiefolder}G_mapped_bam_GFP_sorted/${filename}_GFP.bam"
samtools index "${bowtiefolder}G_mapped_bam_GFP_sorted/${filename}_GFP_sense.bam"
samtools index "${bowtiefolder}G_mapped_bam_GFP_sorted/${filename}_GFP_antisense.bam"

done

echo "finish extracting GFP reads"


echo "all finished!"