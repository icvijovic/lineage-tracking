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
import numpy.ma as ma

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import lineage.inference_params as inference_params

import lineage.file_parser as file_parser
from lineage.read_clone_data import *
			

parser = argparse.ArgumentParser()
parser.add_argument("pop_name", help='name of population')

args = parser.parse_args()

population = args.pop_name

ENVIRONMENT = 'barcoding'

coverage_directory = config.barcode_data_root_directory + population +'/'

count_file = glob.glob(coverage_directory+population+'*read_coverage.txt')[0]
count_dict = file_parser.parse_count_file(count_file)
counts = numpy.asarray([count_dict[key] for key in sorted(count_dict.keys())])

kappas = file_parser.read_kappas_from_file(config.error_model_directory + population + '-kappas.tsv')

clone_list, times, clone_dict, lineage_dict, population_tree  = read_clone_data(population, read_fitness = False, assign_colors = False)

max_barcode = config.max_barcode[population]

bc_fitness_dict = {}
bc_fitness_ci_dict = {}

converged = False
it = 0
max_iterations = 100
negative_llhs = [10**40]

for clone_ID in clone_dict.keys():
	bc_fitness_dict.update({clone_ID : 0. })
	bc_fitness_ci_dict.update({clone_ID : np.zeros(2) })

population_mean_fitness = numpy.zeros(2*max_barcode-2)

barcoding_kappas = kappas[10::11][:max_barcode-1]

bc_count_ratios = counts[11::11]*1./(counts[10::11][:-1])
bc_count_ratios = bc_count_ratios[:max_barcode-1]


num_fitnesses = inference_params.num_fitnesses
fitness_grid = inference_params.fitness_grid[ENVIRONMENT]*inference_params.scale_fitness_per_interval[ENVIRONMENT]

bc_fitness_dict = {}
bc_fitness_ci_dict = {}

total_intervals = inference_params.BARCODING_INTERVALS_PER_EPOCH*(max_barcode -1)

# initialize clone fitnesses to zero

current_clone_fitnesses = numpy.zeros( (len(clone_dict.keys())) )
current_clone_fitness_CI_lower = numpy.zeros( (len(clone_dict.keys())) )
current_clone_fitness_CI_upper = numpy.zeros( (len(clone_dict.keys())) )
# initialize up matrices fo clone counts before and after barcoding
clone_counts_before = numpy.ones( (len(clone_dict.keys()),total_intervals) )
clone_counts_after = numpy.ones( (len(clone_dict.keys()),total_intervals) )
clone_freqs_before = numpy.ones( (len(clone_dict.keys()),total_intervals) )
clone_freqs_after = numpy.ones( (len(clone_dict.keys()),total_intervals) ) 

for i in range(0,len(clone_dict.keys())):
	lineage = clone_dict[clone_dict.keys()[i]]

	clone_counts_before[i] = (lineage.freqs*counts)[10::11][:-1][:max_barcode-1]
	clone_counts_after[i] = (lineage.freqs*counts)[11::11][:max_barcode-1]

	clone_freqs_before[i] = (lineage.freqs)[10::11][:-1][:max_barcode-1]
	clone_freqs_after[i] = (lineage.freqs)[11::11][:max_barcode-1]

# initialize grid containing expectations of frequencies given fitness vectors
expectation_grid = numpy.ones( (len(clone_dict.keys()),total_intervals,len(fitness_grid)) )

#initialize mask for low frequency counts
mask = clone_counts_before > inference_params.threshold_lineage_size

#now for the coordinate descent

#iterate until convergence or until max_iterations have been completed 
it = -1
while it < max_iterations:
	it += 1	
	#loop through clones
	for this_clone_index in range(0,len(clone_dict.keys())):
		if sum(mask[this_clone_index]) > 0:
			other_clones = numpy.ones(len(clone_dict.keys()),dtype = bool)
			other_clones[this_clone_index] = False
			
			other_lineages_unnormalized_frequencies = np.einsum('ij,i->ij',clone_freqs_before[other_clones],numpy.exp(current_clone_fitnesses[other_clones]))
			other_lineages_total_unnormalized_frequency = numpy.sum(clone_freqs_before[other_clones].T * numpy.exp(current_clone_fitnesses[other_clones]),axis = 1)		
			this_lineage_unnormalized_frequency = numpy.outer(clone_freqs_before[this_clone_index],numpy.exp(fitness_grid))
			total_unnormalized_frequency = numpy.outer(other_lineages_total_unnormalized_frequency,numpy.ones(len(fitness_grid))) + this_lineage_unnormalized_frequency

			expectation_grid[other_clones] = np.repeat(other_lineages_unnormalized_frequencies[:, :, np.newaxis], len(fitness_grid), axis=2)
			expectation_grid[this_clone_index] = this_lineage_unnormalized_frequency

			expectation_grid = expectation_grid/total_unnormalized_frequency
			expectation_grid = np.einsum('ijk,j->ijk', expectation_grid, counts[11::11][:max_barcode-1])

			expectation_grid[expectation_grid < 1] = 1
			clone_counts_after[clone_counts_after<1] = 1

			#calculate negative log likelihoods
			llhs = (np.sqrt(expectation_grid) - numpy.sqrt(clone_counts_after[:,:,np.newaxis]))**2
			llhs = np.einsum('ijk,j->ijk',llhs,1./barcoding_kappas)
			llhs += 0.75* numpy.log(clone_counts_after[:,:,np.newaxis])
			llhs += -0.25* numpy.log(expectation_grid) 
			llhs += 0.5 * numpy.log(4*numpy.pi*barcoding_kappas[np.newaxis,:,np.newaxis])


			llhs = np.einsum('ijk,ij->ijk',llhs,mask)

			llhs = np.einsum('ijk->k',llhs)

			current_clone_fitnesses[this_clone_index] = fitness_grid[llhs == min(llhs)][0]
			fitness_CI_range = fitness_grid[llhs < 2+ min(llhs)]
			current_clone_fitness_CI_lower[this_clone_index] = fitness_CI_range[0]
			current_clone_fitness_CI_upper[this_clone_index] = fitness_CI_range[-1]

			if clone_dict.keys()[this_clone_index] == "":
				# pin ancestor fitness at 0
				current_clone_fitnesses -= current_clone_fitnesses[this_clone_index]
				current_clone_fitness_CI_lower -= current_clone_fitnesses[this_clone_index]
				current_clone_fitness_CI_upper -= current_clone_fitnesses[this_clone_index]

			if len(fitness_grid[llhs == min(llhs)]) > 1:
				print len(fitness_grid[llhs == min(llhs)]), "solutions work equally well"
				print sum(mask[this_clone_index])
		else:
			# if there are no good barcoding intervals,
			# look for ancestor, and copy their fitness onto this child
			ancestors = list(ancestor_list(clone_dict.keys()[this_clone_index], population_tree))
			if len(ancestors)>0:
				parent_ID = ancestors[-1]
			else:
				parent_ID = ''
			parent_index = [i for i, s in enumerate(clone_dict.keys()) if parent_ID == s][0]
			current_clone_fitnesses[this_clone_index] = current_clone_fitnesses[parent_index]
			current_clone_fitness_CI_lower[this_clone_index] = fitness_grid[0]
			current_clone_fitness_CI_upper[this_clone_index] = fitness_grid[-1]		

	negative_llhs.append(min(llhs))
	if abs(negative_llhs[-1] - negative_llhs[-2])<0.001:
		print 'Converged after', it, 'iterations.'
		print 'Negative log-likelihood for each iteration:'
		print negative_llhs
		break
for this_clone_index in range(0,len(clone_dict.keys())):
	ID = clone_dict.keys()[this_clone_index]
	bc_fitness_dict[ID] = current_clone_fitnesses[this_clone_index]
	bc_fitness_ci_dict[ID] = [current_clone_fitness_CI_lower[this_clone_index],current_clone_fitness_CI_upper[this_clone_index]]

with open(config.clone_data_directory+'%s-%s_fitnesses.tsv' % (population,ENVIRONMENT), 'w') as csvfile:
	out_writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	out_writer.writerow(['BC'] + ['Barcoding Fitness (per barcoding procedure)', 'CI_lower', 'CI_upper'])

	clone_list.append("")
	for ID in clone_list:
		lineage = clone_dict[ID]
		if lineage.ID == '':
			row = ['ancestor']
		else:
			row = [lineage.ID]
		row.extend([bc_fitness_dict[lineage.ID]])
		row.extend(bc_fitness_ci_dict[lineage.ID])
		out_writer.writerow(row)


