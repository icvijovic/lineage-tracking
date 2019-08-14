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


from scipy.stats import linregress

# import avg_vs_establishment_per_epoch

fontsize = 8

# set up figures and axes
fig,ax = plt.subplots(figsize = (2,2))

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
	print population

	max_barcode = config.max_barcode[population]

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

	clone_list.append("")

	effect_dict = {}

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

	fitness_effects = []
	max_bc_CI = 0
	for ID in clone_list:
		if ID == "":
			pass
		else:
			parent_ID = find_last_parent(ID, population_tree)
			ev_diff = clone_dict[ID].evolution_fitness - clone_dict[parent_ID].evolution_fitness
			bc_diff = clone_dict[ID].barcoding_fitness - clone_dict[parent_ID].barcoding_fitness

			effect_dict.update({ID:[100.*ev_diff,bc_diff]})
			max_bc_CI = max(max_bc_CI, abs(clone_dict[ID].barcoding_fitness_CI[1] - clone_dict[ID].barcoding_fitness_CI[0]))

	ev_effects = []
	bc_effects = []
	for ID in effect_dict.keys():
		bc_CI = abs(clone_dict[ID].barcoding_fitness_CI[1] - clone_dict[ID].barcoding_fitness_CI[0])
		if bc_CI < max_bc_CI:
			symbol = 'o'
			alpha = 0.8
			ev_effects.append(effect_dict[ID][0])
			bc_effects.append(effect_dict[ID][1])

			ax.scatter(effect_dict[ID][0],effect_dict[ID][1],marker = symbol, color = config.pop_colors[population],
															 s = 16, lw = 0., alpha = alpha, clip_on=False)
	ev_effects = np.asarray(ev_effects)
	bc_effects = np.asarray(bc_effects)
	slope,intercept,rsquare,pvalue,stderr = linregress(ev_effects,bc_effects)
	print 'slope = ',slope
	print 'intercept = ', intercept
	print 'rsquare = ', rsquare
	print 'pvalue = ', pvalue

	if 'D' in population:
		print '\n  repeating regression with only mutations for which evolution env. effect < 2.5\n'
		slope,intercept,rsquare,pvalue,stderr = linregress(ev_effects[ev_effects<2.5],bc_effects[ev_effects<2.5])
		print 'slope = ',slope
		print 'intercept = ', intercept
		print 'rsquare = ', rsquare
		print 'pvalue = ', pvalue

	ax.plot([],[],lw = 8, color= config.pop_colors[population], label = config.pop_labels[population])
ax.legend(loc = 3,handlelength =  0.8)
ax.set_ylim(ax.get_xlim())
ax.set_xlabel('Evolution fitness effect (per epoch)',fontsize = 8)
ax.set_ylabel('Barcoding fitness effect (per interval)', fontsize = 8)

OUTFILE = config.figure_directory+'si/jdfe.pdf'
fig.savefig(OUTFILE,bbox_inches = 'tight')



	

		