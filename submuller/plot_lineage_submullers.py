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

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc

import lineage.file_parser as file_parser
from lineage.plot_utils import * 

fitness_cmap = matplotlib.cm.get_cmap('PuOr')


matplotlib.rcParams['font.size'] = 16

def assign_muller_lower(epoch = 0, parent_id = '',low = None,scale = None):
	""" Assigns muller diagram positions to all lineages for quick plotting later. """
	for bcd, lineage in data[epoch][parent_id].iteritems():	
		if low is None:
			low = numpy.zeros(len(lineage.freqs))
		if scale is None:
			scale = numpy.ones(len(lineage.freqs))	
		scale[scale<10**-6] = 10**-6.

		lineage.muller_lower = deepcopy(low)

		if epoch + 1 <= len(data)-1:
			if lineage.ID in data[epoch+1]:
				assign_muller_lower(epoch+1,lineage.ID, low = deepcopy(low),scale = scale)
		low += (lineage.freqs)/(scale)


parser = argparse.ArgumentParser()
parser.add_argument("pop_name", help='name of population')
parser.add_argument("parental_barcode",help = 'parental barcode')
parser.add_argument('--fitness', dest='color_by_fitness', action='store_true')
parser.add_argument('--color', dest='color_by_fitness', action='store_false')
parser.add_argument('-highlight',help = 'lineage to highlight', required = False)
parser.add_argument('-stopat',help = 'number of barcodes after focal to plot',required = False,type = int)
parser.set_defaults(color_by_fitness='both',stopat=-1)

args = parser.parse_args()

population = args.pop_name
target_parental_barcode = args.parental_barcode

color_by_fitness = args.color_by_fitness
highlight = [args.highlight]
stopat = args.stopat

if color_by_fitness == 'both':
	options = [False, True]
else:
	options = [color_by_fitness]

# read data
fitness_files = file_parser.natural_sort(glob.glob(config.lineage_fitness_estimate_directory+population+'*BC*_evolution_fitnesses.csv'))
fitness_files.extend(file_parser.natural_sort(glob.glob(config.lineage_fitness_estimate_directory+population+'*BC*_barcoding_fitnesses.csv')))
timepoints,data,counts = file_parser.get_data(population,config.barcode_data_root_directory,fitness_files = fitness_files)
max_barcode = config.max_barcode[population]
times = timepoints[0]

flag_files = file_parser.natural_sort(glob.glob(config.lineage_flag_directory+population+'*BC*_flags.tsv'))
for filename in flag_files:
	f = open(filename,'r')
	f.readline()
	for row in csv.reader(f,delimiter = '\t'):
		split_ID = row[0].split("_")
		flags = row[1].split(";")

		parent_ID = "_".join(split_ID[:-1])
		last_barcode = split_ID[-1]
		n_bc = len(split_ID) - 1

		data[n_bc][parent_ID][last_barcode].flags = flags

# give barcoding a duration of 100 generations
for epoch in range(1,max_barcode):
	times[epoch*11:]+= 90

dataset = len(target_parental_barcode.split("_"))
if stopat == -1:
	stopat = max_barcode
else:
	stopat = dataset + stopat + 1

grand_parental_barcode = "_".join(target_parental_barcode.split("_")[:-1])
last_barcode = target_parental_barcode.split("_")[-1]
parental_lineage = data[dataset-1][grand_parental_barcode][last_barcode]

#plot submuller
parental_freq_vector = parental_lineage.freqs
scale = parental_freq_vector

assign_muller_lower(epoch = dataset, parent_id = target_parental_barcode,low = numpy.zeros(len(parental_freq_vector)), scale = parental_freq_vector)

color_cycler = ColorCycler()

sys.stderr.write("initializing figure and axes objects....\n")
#set up figure and axis objects
fig_muller,ax = plt.subplots(figsize = (max_barcode*2,12))
outer_grid = gridspec.GridSpec(3,1, height_ratios = [1,4,1])
inner_grid_1 = gridspec.GridSpecFromSubplotSpec(1,2*max_barcode-1,subplot_spec = outer_grid[0],wspace = 0.4)
ax_muller = plt.subplot(outer_grid[1])
ax_freq = plt.subplot(outer_grid[2])
ax_fitness = [plt.subplot(inner_grid_1[i]) for i in range(0,max_barcode*2-1)]
ax_muller = plt.subplot(outer_grid[1])

	
bar_chart = [BarChart(fig_muller,ax_fitness[i],numpy.linspace(-5,5,31)) for i in range(0,len(ax_fitness))]
muller_diagram = MullerDiagram(fig_muller,ax_muller,times)

ax_freq.plot(times,parental_lineage.freqs,color = 'k')
ax_freq.set_yscale('log')
ax_freq.set_xlim(ax_muller.get_xlim())

for color_by_fitness in options:
	ax_muller.clear()
	ax_muller.set_ylim([0.,1.])
	for barcode_level in range(dataset-1,stopat):
		
		# create new bar chart axess
		for axis in ax_fitness:
			axis.clear()
		
		for i in range(0,len(ax_fitness)):
			bar_chart = [BarChart(fig_muller,ax_fitness[i],numpy.linspace(-5,5,31)) for i in range(0,len(ax_fitness))]

		if barcode_level == dataset-1:
			for epoch in range(0,max_barcode-1):
				rectangle_right = 100+200*epoch
				ax_muller.add_patch(patches.Rectangle(
			    	(rectangle_right, -0.1), 100, 1.1,
			    	facecolor= 'k',alpha = 0.2))

		if barcode_level < dataset:
			descendant_barcodes = [grand_parental_barcode]
			lineages_to_plot = [data[barcode_level][grand_parental_barcode][last_barcode]]
		else:
			descendant_barcodes = [p for p in data[barcode_level].keys() if p.startswith(target_parental_barcode)]
			lineages_to_plot = [data[barcode_level][d][l] for d in descendant_barcodes for l in data[barcode_level][d].keys()]
		
		for lineage in lineages_to_plot:

			large_enough = max((lineage.freqs+10**-10)/(scale+10**-6)) > 10**-3
			lineage_ever_selected = isinstance(lineage.relative_fitness,dict)

			if large_enough:

				# determine color and opacity of lineage
				if barcode_level <= dataset or color_by_fitness:
					alpha_value = 1.0
				else:
					alpha_value = 0.5

				if lineage.color == None:
					lineage.color = color_cycler.get_new_color(alpha = alpha_value)
				
				# determine whether to outline lineage
				if lineage.ID in highlight:
					outline = True
				else:
					outline = False

				if ('b' in lineage.flags) and ('x' not in lineage.flags):
					pass

				# place lineage on barchart and muller
				for in_epoch in range(barcode_level,max_barcode):
					if color_by_fitness:
						if lineage_ever_selected:
							fitness_color = fitness_cmap(0.5+lineage.relative_fitness['evolution'][in_epoch]/5)
							barcoding_fitness_color = fitness_cmap(0.5+lineage.relative_fitness['barcoding'][in_epoch]/5)
							
							if (abs(lineage.relative_fitness["evolution"][in_epoch])>10**-4):
								bar_chart[in_epoch*2].add_lineage(lineage,in_epoch,environment = "evolution",
																	color = fitness_color,
																	no_error_bar = False,
																	scale = scale)
							
							if (abs(lineage.relative_fitness["barcoding"][in_epoch])>10**-4):
								bar_chart[in_epoch*2+1].add_lineage(lineage,in_epoch,environment = "barcoding",
																	color = barcoding_fitness_color,
																	no_error_bar = False,
																	scale = scale)
						else:
							fitness_color = fitness_cmap(0.5)
							barcoding_fitness_color = fitness_cmap(0.5)
						muller_kwargs = dict(ignore_epochs = in_epoch,alpha = 1,scale = scale,
											 outline = outline,single_epoch = True)
						muller_diagram.place(lineage,color=fitness_color,environment="evolution",**muller_kwargs)
						muller_diagram.place(lineage,color=barcoding_fitness_color,environment="barcoding",**muller_kwargs)
					else:
						muller_kwargs = dict(ignore_epochs = in_epoch, scale = scale, outline = outline, single_epoch = True)
						muller_diagram.place(lineage,environment="evolution",**muller_kwargs)
						muller_diagram.place(lineage,environment="barcoding",**muller_kwargs)
						
						if lineage_ever_selected:

							if (abs(lineage.relative_fitness['evolution'][in_epoch])>10**-4):
								bar_chart[in_epoch*2].add_lineage(lineage,in_epoch, environment = "evolution", 
																		no_error_bar = False,scale = scale)
							if (abs(lineage.relative_fitness['barcoding'][in_epoch])>10**-4):
								bar_chart[in_epoch*2+1].add_lineage(lineage,in_epoch,environment = "barcoding",
																		no_error_bar = False,scale = scale)
			


		ax_muller.set_xlim([0,200*max_barcode-100])

		if color_by_fitness:
			FILENAME = 'figures/%s_muller_%s_%dBC_fitness_color.png' % (population,target_parental_barcode,barcode_level)
		else:
			FILENAME = 'figures/%s_muller_%s_%dBC_ID_color.png' % (population,target_parental_barcode,barcode_level)
		muller_diagram.save(FILENAME)

