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
from lineage.read_clone_data import *

parser = argparse.ArgumentParser()
parser.add_argument("-dpi", default = 600, type = int, help = 'resolution (dots per inch)')
parser.add_argument("-width", default = 3.4, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 7.65, type = float, help = 'figure height (inches)')
parser.add_argument("-fontsize", default = 7, type = int, help = 'fontsize (pt)')
parser.add_argument("--stretched",dest='stretched', action='store_true',default=True,
                    help='display barcoding intervals as 100 generations long')


args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize
stretched = args.stretched


matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize

fig_muller,ax_muller = plt.subplots(3,1,figsize = (fig_width,fig_height))
i = -1

# display barcoding as length_of_barcoding_generations longs
if stretched:
	length_of_barcoding = 100
else:
	length_of_barcoding = 10

#time to add to 10 generations
additional_length_of_barcoding = length_of_barcoding - 10

lineage_color_dict = {}

for population in ['D1','D1R1','D1R2']:

	if 'R' in population:
			max_barcode = 8
	else:
		max_barcode= config.max_barcode[population]

	sys.stderr.write(" reading freq files...\n")
	directory = config.barcode_data_root_directory + population +'/'
	lineage_freq_files = file_parser.natural_sort(glob.glob(directory+population+'*BC*_frequencies.txt'))

	if population == 'D1':
		clone_list, clone_times, clone_dict, clone_lineage_dict, clone_tree = read_clone_data(population)
		for epoch in range(1,max_barcode):
			clone_times[epoch*11:] += additional_length_of_barcoding

	timepoints = []
	lineage_data = []
	for lineage_freq_file in lineage_freq_files:
		times, dataset = file_parser.parse_frequency_file(lineage_freq_file)
		timepoints.append(times)
		lineage_data.append(dataset)
	times = timepoints[0]
	for epoch in range(1,max_barcode):
		times[epoch*11:] += additional_length_of_barcoding


	freq_cutoff = 10**-3
	lineage_list = []
	population_tree = {}
	lineage_dict = {}
	for epoch in range(0,len(lineage_data)):
		for parental_bcd in lineage_data[epoch].keys():
			for bcd, lineage in lineage_data[epoch][parental_bcd].iteritems():
				if lineage.ID in clone_dict or ('R' in population and len(ID.split("_")) > 7 and any(lineage.freqs > freq_cutoff)):
					lineage_list.append(lineage.ID)
					insert_into_tree(lineage.ID,population_tree)
					lineage_dict.update({lineage.ID:lineage})

	i += 1

	sys.stderr.write("initializing figure and axes objects....\n")

	if population == 'D1':
		# times[times > 870] += 100
		max_time = times[-1]
	muller_diagram = MullerDiagram(fig_muller,ax_muller[i],times)

	max_freqs = compile_max_freqs(population_tree.keys(),lineage_dict)
	if population == 'D1':
		sorted_ids = [x for y,x in reversed(sorted(zip(max_freqs,population_tree.keys())))]

	assign_muller_lower(sorted_ids,lineage_dict,population_tree,
		baseline = numpy.zeros(len(times)),
			range_available = numpy.ones(len(times)))

	number_of_lineages_with_children = len([key for key in population_tree.keys() if len(population_tree[key].keys())>0])
	start_index = 8 # 0, 5, 7
	small_start_index = 10
	small_alpha = 0.5
	large_color_cycler = ColorCycler(number_of_lineages_with_children,start_index = start_index)
	small_color_cycler = ColorCycler(40,start_index = small_start_index)

	no_barcodes = 1
	for ID in lineage_list:
		if ID in clone_list:
			lineage_dict[ID].color = clone_dict[ID].color
			muller_diagram.place(lineage_dict[ID])
		elif 'R' in population and len(ID.split("_")) > 7:
			if ID in lineage_color_dict.keys():
				lineage_dict[ID].color = lineage_color_dict[ID]
			else:
				lineage_dict[ID].color = small_color_cycler.get_new_color(alpha = max(small_alpha,1.0-0.1*len(ID.split("_"))))
				lineage_color_dict.update({ID:lineage_dict[ID].color})

			muller_diagram.place(lineage_dict[ID])

	for no_barcodes in range(1,max_barcode):

		rectangle_right = (110 + additional_length_of_barcoding)*no_barcodes

		if stretched:
			bar_alpha = 0.3
		else:
			bar_alpha = 0.7
		if no_barcodes > 0:
			ax_muller[i].add_patch(patches.Rectangle(
	        	(rectangle_right - length_of_barcoding, -0.1), length_of_barcoding, 1.1,
	        	facecolor= '0.3',
	        	alpha = bar_alpha,
	        	edgecolor = "none"))

	xticks = []
	xticklabels = []
	if max_barcode < 9:
		max_barcode = 10
	xtick_dt = 110+additional_length_of_barcoding
	for epoch in range(0,max_barcode+1):
		if i == 0:
			xticks.extend([xtick_dt*epoch,xtick_dt*epoch+100])
			xticklabels.extend(['%d.0'%(epoch+1),'%d.100'%(epoch+1)])
		else:
			if epoch < 8:
				xticks.extend([xtick_dt*epoch,xtick_dt*epoch+100])
				xticklabels.extend(['%d.0'%(epoch+1),'%d.100'%(epoch+1)])
			if epoch == 8:
				xticks.extend([xtick_dt*epoch])
				xticklabels.extend(['8.200'])

		
	ax_muller[i].set_xticks(xticks)
	ax_muller[i].set_xticklabels(xticklabels)
	ax_muller[i].set_xlim([0,max_time])

	if i == 2:
		ax_muller[i].set_xlabel('Generation')
	ax_muller[i].set_ylabel('Barcodes')
	label = {'D1':'YPA(R0)','D1R1':'YPA(R1)','D1R2':'YPA(R2)'}
	ax_muller[i].text(-0.15,1.05,label[population],fontsize = 7,weight = 'bold',transform = ax_muller[i].transAxes)

	it_tick = 0
	for tick in ax_muller[i].get_xticklabels():
   		tick.set_rotation(90)
   		if stretched:
   			tick.set_horizontalalignment('center')
   		else:
	   		if it_tick % 2 == 0:
	   			tick.set_horizontalalignment('left')
	   		else:
	   			tick.set_horizontalalignment('right')
	   		it_tick+=1


	muller_diagram.save(config.figure_directory+'si/D1_replicate_muller.png',dpi = dpi)

