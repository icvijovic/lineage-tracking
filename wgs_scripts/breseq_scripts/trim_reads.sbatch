#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 0-10:00
#SBATCH --mem 32000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o trim.out      #File to which standard out will be written
#SBATCH -e trim.err      #File to which standard err will be written

TRIMMOMATIC=/n/regal/desai_lab/icvijovic/Trimmomatic

java -jar ${TRIMMOMATIC}/trimmomatic-0.35.jar PE -phred33 -threads 1 \
${sample_fastq_directory}/${sample_name}.R1.fastq.gz ${sample_fastq_directory}/${sample_name}.R2.fastq.gz \
data/trimmed_fastq_files/${sample_name}.R1.fastq data/trimmed_fastq_files/${sample_name}.R1.unpaired.fastq \
data/trimmed_fastq_files/${sample_name}.R2.fastq data/trimmed_fastq_files/${sample_name}.R2.unpaired.fastq \
ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 \
2> data/trimmed_fastq_files/${sample_name}.trim.paired.txt

touch data/trimmed_fastq_files/${sample_name}.finished