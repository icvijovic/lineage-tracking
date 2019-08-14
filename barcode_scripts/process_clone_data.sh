#!/bin/bash
#

echo "Compiling clone lists based on flags stored in ../data/flags/ ..."
bash compile_clone_list.sh C1
bash compile_clone_list.sh D1

echo "Calculating clone frequencies..."
python calculate_clone_frequencies.py

echo "Jointly inferring evolution and barcoding condition fitness for all clones in population C1 (YPD)..."
python measure_evolution_fitness.py C1
python measure_barcoding_fitness.py C1

echo "Jointly inferring evolution and barcoding condition fitness for all clones in population D1 (YPA)..."
python measure_evolution_fitness.py D1
python measure_barcoding_fitness.py D1

echo "Calcualting p-values for the number of mutations acquired by each clone in each epoch..."
python calculate_segment_p_values.py