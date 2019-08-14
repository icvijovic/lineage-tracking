#!/bin/bash

export population=$1

python population_parameters.py samples_directories ${population} > temp_parameters_${population}

mkdir -p data/breseq_output

while read -r i j; do

    export sample_name=$i

    mkdir -p data/breseq_output/${sample_name}
    if [ ! -e "data/breseq_output/${sample_name}/output/evidence/evidence.gd" ]
    then

    rm -rf data/breseq_output/${sample_name}
    mkdir -p data/breseq_output/${sample_name}

    sbatch --job-name=breseq-${i} breseq_genomes.sbatch
    # bash breseq_genomes.sbatch
    fi
done < temp_parameters_${population}

rm temp_parameters_${population}