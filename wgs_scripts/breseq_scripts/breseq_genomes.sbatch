#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 0-10:00
#SBATCH --mem 32000
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o slurm_out_and_err/breseq.out      #File to which standard out will be written
#SBATCH -e slurm_out_and_err/breseq.err      #File to which standard err will be written

export REFERENCE_FILE='../../data/wgs/genome_feature_files/reference.gff3'

export FASTQ_FILES=""

for f in $(ls data/trimmed_fastq_files/${sample_name}*.fastq)
do
   export FASTQ_FILES="${FASTQ_FILES} ${f}"
done

# echo ${FASTQ_FILES}

breseq -p -j 8 --brief-html-output -o data/breseq_output/${sample_name} -r ${REFERENCE_FILE} ${FASTQ_FILES} > data/breseq_output/${sample_name}.out 2> data/breseq_output/${sample_name}.err
