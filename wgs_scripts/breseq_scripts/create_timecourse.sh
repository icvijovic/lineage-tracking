#!/bin/bash
export population=$1
export bam_files=""
export sample_names=""

for i in $(python population_parameters.py samples ${population}); do
    export sample_name=$i
    export sample_names="${sample_names} ${sample_name}"
    export reference_file="data/breseq_output/${sample_name}/data/reference.fasta"
    export bam_files="${bam_files} data/breseq_output/${sample_name}/data/reference.bam"
done
#echo ${bam_files}
#echo ${reference_file}

bp_per_run=100000

python population_parameters.py chromosomes ${bp_per_run} > temp_chromosome_info

while read -r a b; do
	export chr=${a}
    export chr_length=$b
    for (( i=$bp_per_run; i<=$chr_length ; i+=$bp_per_run));
    do
    	# echo $i
    	export start_position=`expr $i - $bp_per_run`
    	export end_position=$i
    	# echo ${population}_${chr}_${start_position}
        # bash create_timecourse.sbatch
        sbatch --job-name=timecourse_${population}_${chr}_${start_position} create_timecourse.sbatch
    done

done < temp_chromosome_info

rm temp_chromosome_info
