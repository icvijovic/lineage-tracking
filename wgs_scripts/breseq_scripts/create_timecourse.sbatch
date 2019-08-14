#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 16:00:00
#SBATCH --mem 16000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o slurm_out_and_err/create_full_timecourse_${population}.out      #File to which standard out will be written
#SBATCH -e slurm_out_and_err/create_full_timecourse_${population}.err      #File to which standard err will be written

samtools mpileup -q10 -f ${reference_file} -r ${chr}:${start_position}-${end_position} ${bam_files} > data/timecourse_files/${population}_${chr}_${start_position}.pileup

cat data/timecourse_files/${population}_${chr}_${start_position}.pileup \
 | python create_timecourse.py ${start_position} ${end_position} ${population} ${sample_names} > data/timecourse_files/${population}_${chr}_${start_position}_timecourse.txt

# remove temporary file
rm data/timecourse_files/${population}_${chr}_${start_position}.pileup
