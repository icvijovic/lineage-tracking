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
from scipy import stats

# add custom modules to path
sys.path.insert(0,'../modules/') 


# import custom modules
import config
import local_matplotlibrc
import lineage.file_parser as file_parser

from lineage.read_clone_data import *
from lineage.plot_utils_clone import * # contains muller and bar chart class
from lineage.tree_utils import *
import lineage.inference_params

popindex = -1

def find_last_parent(ID,population_tree):
		if ID == "":
			return ''
		else:
			parents = ancestor_list(ID,population_tree)
			i = 0
			for p_ID in parents:
				i+= 1
			if i > 0:
				return p_ID
			else:
				return ''


for population in config.populations:
	popindex += 1

	max_barcode = config.max_barcode[population]


	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

	population_tree = {"":population_tree}
	mutant_tree = population_tree[""]
	clone_list.append("")

	mean_average_fitness = numpy.zeros(inference_params.INTERVALS_PER_EPOCH*max_barcode)
	for ID in clone_list:
		clone = clone_dict[ID]
		mean_average_fitness += 100.* clone.evolution_fitness * clone.freqs[:inference_params.INTERVALS_PER_EPOCH*max_barcode]
		mean_average_fitness += clone.barcoding_fitness * clone.freqs[:inference_params.INTERVALS_PER_EPOCH*max_barcode]
	mean_average_fitness /= 2.

	clone_dict[''].color = 'k'

	no_descendants = numpy.zeros((len(clone_list),max_barcode-1))
	weights = numpy.zeros((len(clone_list),max_barcode-1))
	relative_fitnesses = numpy.zeros((len(clone_list),max_barcode-1))

	total_mutations = len(clone_list) - 1
	average_fitness_effects = numpy.zeros(total_mutations)

	for i in range(0,len(clone_list)):
		clone = clone_dict[clone_list[i]]
		arising_time = len(clone_list[i].split("_")) - 1
		ID = clone.ID

		if clone_list[i] != "":
			descendants = child_list(clone_list[i],population_tree)

			parent_ID = find_last_parent(ID, population_tree)
			ev_diff = clone_dict[ID].evolution_fitness - clone_dict[parent_ID].evolution_fitness
			bc_diff = clone_dict[ID].barcoding_fitness - clone_dict[parent_ID].barcoding_fitness
			
			average_fitness_effects[i] = (100.*ev_diff+bc_diff)/2.
		else:
			descendants = mutant_tree.keys()

		for epoch in range(2,max_barcode-1):
			no_descendants[i][epoch] = sum(1 for item in descendants if len(item.split("_")) == epoch + 1)
			weights[i][epoch] = numpy.sum(clone.freqs[11*epoch:11*(epoch+1)])

			average_fitness = (clone.evolution_fitness*100. + clone.barcoding_fitness)/2.
			relative_fitnesses[i][epoch]= average_fitness - sum(mean_average_fitness[11*epoch:11*(epoch+1)])/11.
	

	old_weights = deepcopy(weights)
	old_shape = weights.shape

	frequency_weights = weights.flatten()
	# frequency_weights /= sum(frequency_weights)
	
	if (frequency_weights[frequency_weights<0])<-10**-4:
		print "Warning! Found a lineage with negative frequency."
	frequency_weights[frequency_weights<0] = 0.

	mutation_weights = numpy.ones((total_mutations,len(frequency_weights)))
	mutation_weights *= frequency_weights

	for model in ['frequency','frequency_plus_fitness']:
		print population, model
		for mut_it in range(0,total_mutations):

			if model == 'frequency_plus_fitness':
				fit_weights = relative_fitnesses.flatten()
				fit_weights += average_fitness_effects[mut_it]
				percent_threshold = 5*10**-5
				fit_weights[fit_weights < 100*percent_threshold] = 0.
				mutation_weights[mut_it]*= fit_weights

			mutation_weights[mut_it] /= sum(mutation_weights[mut_it])

		# resample all mutations
		number_of_trials = 10**4
		mutations_in_trial = numpy.zeros((number_of_trials,len(clone_list)*(max_barcode-1)))


		for n in range(0,number_of_trials):
			for mut_it in range(0,total_mutations):
				w = mutation_weights[mut_it]
				mutant_background = numpy.random.choice(numpy.arange(0,len(w)),size = 1,p = w)
				mutations_in_trial[n][mutant_background] += 1

		no_descendants = no_descendants.flatten()
		p_values_more = numpy.zeros(len(clone_list)*(max_barcode-1))
		for i in range(0,len(p_values_more)):
			p_values_more[i] = sum(1 for trial in range(0, number_of_trials) if mutations_in_trial[trial][i] >= no_descendants[i])
			p_values_more[i] /= 1.*number_of_trials


		p_values_less = numpy.zeros(len(p_values_more))
		for i in range(0,len(p_values_less)):
			p_values_less[i] = sum(1 for trial in range(0, number_of_trials) if mutations_in_trial[trial][i] <= no_descendants[i])
			p_values_less[i] /= 1.*number_of_trials

		with open(config.clone_data_directory+'%s-p_values_segment_%s.tsv' % (population,model), 'w') as csvfile:
			out_writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	   		out_writer.writerow(['BC'] + ['P-gtreq'] + ['P-lesseq'])

	 		for i in range(0,len(clone_list)):
	 			ID = clone_list[i]
	 			row = [ID]
	 			first = i*(max_barcode-1)
	 			last = (i+1)*(max_barcode-1)
	 			row.extend(p_values_more[first:last])
	 			row.extend(p_values_less[first:last])
	 			out_writer.writerow(row)

	 	# now calculate the same for each clone
		order = 'C'
		mutations_in_trial = mutations_in_trial.reshape((number_of_trials,len(clone_list),max_barcode-1),order = order)
	 	no_descendants = no_descendants.reshape((len(clone_list),max_barcode-1))
	 	total_mutations_per_clone = numpy.sum(mutations_in_trial,axis = 2)
		experiment_mutations_per_clone = numpy.sum(no_descendants,axis = 1)

		p_values_more = numpy.zeros(len(clone_list))
		for i in range(0,len(clone_list)):
			p_values_more[i] = sum(1 for trial in range(0, number_of_trials) if total_mutations_per_clone[trial][i] >= experiment_mutations_per_clone[i])
			p_values_more[i] /= 1.*number_of_trials


		p_values_less = numpy.zeros(len(clone_list))
		for i in range(0,len(clone_list)):
			p_values_less[i] = sum(1 for trial in range(0, number_of_trials) if total_mutations_per_clone[trial][i] <= experiment_mutations_per_clone[i])
			p_values_less[i] /= 1.*number_of_trials

		with open(config.clone_data_directory+'%s-p_values_%s.tsv' % (population,model), 'w') as csvfile:
			out_writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	   		out_writer.writerow(['BC'] + ['P-gtreq'] + ['P-lesseq'])

	 		for i in range(0,len(clone_list)):
	 			ID = clone_list[i]
	 			row = [ID]
	 			row.extend([p_values_more[i]])
	 			row.extend([p_values_less[i]])
	 			out_writer.writerow(row)
