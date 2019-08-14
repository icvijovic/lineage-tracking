import sys, os, glob, csv, re
import string, math, numpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import matplotlib.patches as patches
import argparse
from scipy.special import erf
from copy import deepcopy
import matplotlib.gridspec as gridspec
import itertools

import matplotlib.animation as manimation

## add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import file_parser as file_parser
	
from lineage import *
import inference_params

from tree_utils import *
from plot_utils_clone import *

clone_directory = config.clone_data_directory

def read_clone_data(population, read_fitness = True, assign_colors = True):

	if population == "D1":
		fitness_norm = 0.12
	else:
		fitness_norm = 0.05

	clone_list_files = file_parser.natural_sort(glob.glob(clone_directory+population+'*clone_list.tsv'))

	data = []
	sys.stderr.write(" reading freq files...\n")
	clone_freq_file = file_parser.natural_sort(glob.glob(clone_directory+population+'-clone_frequencies.tsv'))[0]
	times, clone_data = file_parser.parse_frequency_file(clone_freq_file,prepend_zeros = False)
	clone_dict = clone_data['']

	sys.stderr.write(" reading clone list\n")
	clone_list = []

	population_tree = {}
	for filename in clone_list_files:
		f = open(filename,'r')
		for row in csv.reader(f,delimiter = '\t'):
			clone_list.append(row[0])
			insert_into_tree(row[0],population_tree)
		f.close()

	#check tree doesn't have major flaws
	for key in population_tree.keys():
		ancestral_barcode_options = ["_".join(key.split("_")[0:l]) for l in range (1,len(key.split("_")))]
		parents = [p_ID for p_ID in clone_list if p_ID in ancestral_barcode_options]
		if len(parents) > 0:
			print "Error!"

	lineage_dict = deepcopy(clone_dict)

	for lineage_id, lineage in lineage_dict.iteritems():

		ancestral_barcode_options = ["_".join(lineage_id.split("_")[0:l]) for l in range (1,len(lineage_id.split("_")))]
		parents = [p_ID for p_ID in clone_list if p_ID in ancestral_barcode_options]
		children = [p_ID for p_ID in clone_list if (p_ID.startswith(lineage_id) and p_ID > lineage_id)]
		for child_id in children:
			lineage_dict[lineage_id].freqs += clone_dict[child_id].freqs


	total_freqs = numpy.zeros(len(times))
	for basal_lineage in population_tree.keys():
		total_freqs += lineage_dict[basal_lineage].freqs
	wt_freqs = numpy.ones(len(total_freqs)) - total_freqs
	clone_dict.update({"":Lineage("",wt_freqs)})
	clone_dict[""].color = 'w'
	mutant_tree = population_tree

	if read_fitness:
		evolution_fitness_file = config.clone_data_directory+population+'-evolution_fitnesses.tsv'
		barcoding_fitness_file = config.clone_data_directory+population+'-barcoding_fitnesses.tsv'

		for clone_fitness_file in [evolution_fitness_file,barcoding_fitness_file]:
			if len(clone_fitness_file) > 0:
				sys.stderr.write(" reading fitness inferences...\n")
				# print clone_fitness_file
				f = open(clone_fitness_file,'r')
				f.readline()
				for row in csv.reader(f,delimiter='\t'):
					bcd = row[0]
					if bcd == "ancestor":
						bcd = ""
					fitnesses = numpy.asarray(row[1:],dtype = float)
					target_lineage = clone_dict[bcd]
					if 'barcoding' in clone_fitness_file:
						target_lineage.barcoding_fitness = fitnesses[0]
						target_lineage.barcoding_fitness_CI = (fitnesses[1],fitnesses[2])
					else:
						SCALE = inference_params.scale_fitness_per_interval['evolution']
						target_lineage.evolution_fitness = fitnesses[0]/SCALE
						target_lineage.evolution_fitness_CI = (fitnesses[1]/SCALE,fitnesses[2]/SCALE)

				f.close()
			else:
				sys.stderr.write("warning! no fitness file found. all clone fitnesses are zero.")


	if assign_colors:

		big_color_divisons = len([key for key in population_tree.keys() if len(population_tree[key].keys())>0])

		max_freqs = compile_max_freqs(population_tree.keys(),lineage_dict)
		sorted_ids_by_freqs = [x for y,x in reversed(sorted(zip(max_freqs,population_tree.keys())))]

		def sort_indices_for_muller(list_to_sort):

			if len(list_to_sort) < 2:
				return list_to_sort
			else:
				middle_item = list_to_sort[0]
				list_above = list_to_sort[1::2]
				list_below = list_to_sort[2::2]

				return sort_indices_for_muller(list_below) + [middle_item] + sort_indices_for_muller(list_above)
		
		sorted_lineages = sort_indices_for_muller(sorted_ids_by_freqs)
		sorted_colored_lineages = [x for x in sorted_lineages if len(population_tree[x].keys()) > 0]

		if 'C1' in population:
			big_color_divisons = 8
			colors = {0:7.5, 1:0, 2:4.8, 3:6,
					  4:3, 5:1, 6:3.5, 7:4,
					  8:7, 9:2, 10:6.5}

		elif 'D1' in population:
			big_color_divisons = 11
			colors = {0:1, 1:4, 2:2, 3:5, 4:5.5,
					  5:0.5, 6:0, 7:6, 8:2, 9:8,
					  10:9, 11:10}


		color_indices = {	sorted_colored_lineages[i] : colors[i]
					 		for i in range(len(sorted_colored_lineages))}

		large_color_cycler = ColorCycler(big_color_divisons)


		# assign clone colors
		for ID in clone_list:
			fitness_alpha = clone_dict[ID].evolution_fitness/fitness_norm*0.6 +0.4
			if fitness_alpha > 1.:
				fitness_alpha = 1.
			# print fitness_alpha

			if ID in population_tree.keys():
				if len(population_tree[ID].keys()) > 0:
					lineage_dict[ID].color = large_color_cycler.get_color_by_index(color_indices[ID],alpha = fitness_alpha)
					clone_dict[ID].color = lineage_dict[ID].color
				else:
					lineage_dict[ID].color = (.5,.5,.5,fitness_alpha)
					clone_dict[ID].color = lineage_dict[ID].color
			else:
				parents = ancestor_list(ID,population_tree)
				parents = list(parents)
				last_parent_ID = parents[-1]

				R, G, B, A = clone_dict[last_parent_ID].color
				lineage_dict[ID].color = (R,G,B,fitness_alpha)
				clone_dict[ID].color = lineage_dict[ID].color

	return clone_list, times, clone_dict, lineage_dict, population_tree



