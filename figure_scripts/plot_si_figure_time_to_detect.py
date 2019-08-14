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

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import lineage.file_parser as file_parser

from lineage.read_clone_data import *
from lineage.plot_utils_clone import * # contains muller and bar chart class
from lineage.tree_utils import *

# import avg_vs_establishment_per_epoch

fontsize = 8

pop_colors = config.pop_colors


# set up figures and axes
fig,ax = plt.subplots(1,1,figsize = (2,1.5))



all_axes = [ax]
for axes in all_axes:
	axes.patch.set_alpha(0.)
	for axis_position in ['top','right']:
		axes.axes.spines[axis_position].set_visible(False)
	axes.axes.xaxis.set_ticks_position('bottom')
	axes.axes.yaxis.set_ticks_position('left')
	# axes.set_yticklabels([])
	# axes.set_xticklabels([])


popindex = -1
for population in reversed(config.populations):
	popindex += 1
	max_barcode = config.max_barcode[population]

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

	mutant_tree = population_tree


	nums_of_epochs = numpy.zeros(max_barcode)

	overall_time_to_detect = {}

	barcode_fitness_files = file_parser.natural_sort(glob.glob(config.lineage_fitness_estimate_directory+population+'-BC*_fitnesses.csv'))
	sys.stderr.write(" reading fitness inferences...\n")
	for filename in barcode_fitness_files:
		f = open(filename,'r')
		f.readline()
		for row in csv.reader(f,delimiter='\t'):
			ID = row[0]

			if ID in clone_list:
				epoch = len(ID.split("_")) - 1

				fitnesses = numpy.asarray(row[1:],dtype = float)[::3]

				first_nonzero = numpy.nonzero(fitnesses)[0][0]
				delta_epochs = first_nonzero - epoch

				if ID in overall_time_to_detect.keys():
					#check if it was detected earlier in a different environment
					if delta_epochs < overall_time_to_detect[ID]:
						overall_time_to_detect[ID] = delta_epochs
					else:
						pass
				else:
					overall_time_to_detect.update({ID:delta_epochs})

		f.close()	

	nums = numpy.arange(max_barcode)
	print "n", "any", "ev", "bc"

	any_condition = numpy.zeros(max_barcode)

	for n in nums:
		any_condition[n] = sum([1 for x in overall_time_to_detect.values() if x ==n])

	percent_detected_any = numpy.cumsum(any_condition)/numpy.cumsum(any_condition)[-1]
	print percent_detected_any

	ax.plot(numpy.arange(max_barcode),percent_detected_any,'-',lw = 2, color = pop_colors[population],label = config.pop_labels[population],clip_on = False,alpha = 0.8)

	ax.set_xlabel('Epochs from founding until\nnull first rejected',labelpad = 1,fontsize = fontsize)
	ax.set_ylabel('Percent detected',labelpad = 1,fontsize = fontsize)

	ax.set_ylim([0,1])

ax.legend(scatterpoints = 1,labelspacing = .8,handlelength = 1,loc = 'best',fontsize = fontsize)
fig.savefig(config.figure_directory+'si/time_to_detect.pdf',bbox_inches = 'tight')



	

# # 		