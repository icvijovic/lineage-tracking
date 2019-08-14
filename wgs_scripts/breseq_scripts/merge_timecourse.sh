#!/bin/bash
export population=$1
export sample_names=""

mkdir -p data/merged_timecourse_files/

# write everything to a merged file
cat data/timecourse_files/${population}_*timecourse.txt | bzip2 -c > data/merged_timecourse_files/${population}_merged_timecourse.bz2