#!/bin/bash

export population=$1

python population_parameters.py samples_directories ${population} > temp_parameters_${population}

mkdir -p data/trimmed_fastq_files

while read -r i j; do

    export sample_name=$i
    export sample_fastq_directory=$j

    mkdir -p data/breseq_output/${sample_name}
    if [ ! -e "data/trimmed_fastq_files/${sample_name}.finished" ]
    then

    sbatch --job-name=trim-${i} trim_reads.sbatch
    # bash trim_reads.sbatch
    fi
done < temp_parameters_${population}

rm temp_parameters_${population}