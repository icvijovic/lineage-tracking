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
from lineage.plot_utils import * # contains muller and bar chart class
from lineage.tree_utils import *


parser = argparse.ArgumentParser()
parser.add_argument("-dpi", default = 600, type = int, help = 'resolution (dots per inch)')
parser.add_argument("-width", default = 7, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 3, type = float, help = 'figure height (inches)')
parser.add_argument("-fontsize", default = 7, type = int, help = 'fontsize (pt)')

args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize

matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize


fig_mean,ax_mean = plt.subplots(2,1,figsize = (3,4))

pop_index = 0

for population in config.populations:
	
	timepoints, data, counts = file_parser.get_data(population,config.barcode_data_root_directory,
													fitness_files = None,as_matrix = True)
	times = timepoints[0]

	max_barcode = config.max_barcode[population]

	averages = numpy.zeros(10)
	nums = numpy.zeros(10)

	CI_upper = numpy.zeros(10)
	CI_lower = numpy.zeros(10)

	bins = numpy.logspace(0,3,15,base =10)
	for dset in range(0,max_barcode):
		this_dataset_counts = data[dset].T[0]*float(counts[dset*11])
		nums[dset]=len(this_dataset_counts)
		averages[dset] = numpy.sum(this_dataset_counts)/nums[dset]

		this_dataset_counts = sorted(this_dataset_counts)
		
		uniq_counts = numpy.unique(this_dataset_counts)
		for x in uniq_counts:
			if sum([1. for item in this_dataset_counts if item <= x])/nums[dset] >= 0.025:
				break
		CI_lower[dset] = x
		for x in reversed(uniq_counts):
			if sum([1 for item in this_dataset_counts if item >= x]) >= 0.025 * nums[dset]:
				break
		CI_upper[dset] = x

		print CI_lower[dset], averages[dset], CI_upper[dset]

	pop_index += 1
	offset = 1 -0.5 + pop_index*0.25

	ax_mean[0].scatter(numpy.arange(0,10)+offset,CI_lower+1,marker = '_',color = config.pop_colors[population],s = 10,clip_on = False)
	ax_mean[0].scatter(numpy.arange(0,10)+offset,CI_upper,marker = '_',color = config.pop_colors[population],s = 10,clip_on = False)
	ax_mean[0].scatter(numpy.arange(0,10)+offset,averages,marker = 'o',color = config.pop_colors[population],s = 10,clip_on = False)
	for dset in range(0,10):
		ax_mean[0].plot([dset+offset,dset+offset],[CI_lower[dset]+1,CI_upper[dset]],lw = 1,color = config.pop_colors[population],clip_on = False)
	ax_mean[1].bar(numpy.arange(0,10)+offset,nums,color = config.pop_colors[population],lw = 0,width = 0.25)
	ax_mean[1].set_xlabel('Epoch',fontsize = fontsize)
	ax_mean[0].set_ylabel('Mean reads per lineage',fontsize = fontsize)
	ax_mean[1].set_ylabel('No. new lineages',fontsize = fontsize)

	ax_mean[0].set_ylim([1,200])
	ax_mean[0].set_yscale('log')
	ax_mean[0].set_xlim([0.5,10.5])
	ax_mean[1].set_xlim([0.5,10.5])

	ax_mean[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
	mf.set_powerlimits((0,0))
	ax_mean[1].yaxis.set_major_formatter(mf)
	ax_mean[1].set_xticks([1,2,3,4,5,6,7,8,9,10])
	ax_mean[0].set_xticks([1,2,3,4,5,6,7,8,9,10])

	ax_mean[1].plot([],[],lw  =8, label = config.pop_labels[population],color = config.pop_colors[population])
	ax_mean[1].legend(loc = 2,fontsize = fontsize - 1)


	fig_mean.savefig(config.figure_directory+'si/mean_counts.pdf',bbox_inches = 'tight')
