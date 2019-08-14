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


parser = argparse.ArgumentParser()
parser.add_argument("-dpi", default = 600, type = int, help = 'resolution (dots per inch)')
parser.add_argument("-width", default = 7.2, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 5.5, type = float, help = 'figure height (inches)')
parser.add_argument("-fontsize", default = 6, type = int, help = 'fontsize (pt)')
parser.add_argument("-outfile", default = 'Figure2_final.pdf', type = str, help = 'figure name')
parser.add_argument("--stretched",dest='stretched', action='store_true',default=True,
                    help='display barcoding intervals as 100 generations long')

args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize
OUTFILE = args.outfile
stretched = args.stretched

# display barcoding as length_of_barcoding_generations longs
if stretched:
	length_of_barcoding = 100
else:
	length_of_barcoding = 10

#time to add to 10 generations
additional_length_of_barcoding = length_of_barcoding - 10

matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize + 1

submuller_lineages= config.highlighted_barcodes

figures_directory = config.figure_directory

if stretched:
	OUTFILE = OUTFILE.split(".")[0]+ "_stretched." + OUTFILE.split(".")[1]
OUTFILE = figures_directory + 'main/' + OUTFILE

populations = config.populations

fig_muller,ax = plt.subplots(figsize = (fig_width,fig_height))
outer_grid = gridspec.GridSpec(3,len(populations), height_ratios = [0,40,30],hspace = 0.45,wspace = 0.15)

popindex = -1
for population in config.populations:
	print population
	popindex += 1

	max_barcode = config.max_barcode[population]

	# set up all axes objects

	inner_grid_1 = gridspec.GridSpecFromSubplotSpec(1,max_barcode,subplot_spec = outer_grid[0,popindex],wspace = 0.3,hspace = 0.1)
	inner_grid_2 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec = outer_grid[2,popindex],wspace = 0.1,hspace = 0.2)

	ax_muller = plt.subplot(outer_grid[1,popindex])

	ax_muller.set_xlabel('Generation',labelpad = 3*(1-popindex))

	ax_submuller = [plt.subplot(inner_grid_2[i]) for i in range(0,4)]
	ax_submuller[0].set_xticklabels([])
	ax_submuller[1].set_xticklabels([])
	ax_submuller[1].set_yticklabels([])
	ax_submuller[3].set_yticklabels([])
	ax_submuller[2].set_xlabel('Generation',labelpad = 0)
	ax_submuller[3].set_xlabel('Generation',labelpad = 0)

	if popindex == 0:
		# ax_muller.text(-0.12,1.13,'Figure 2', fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)

		ax_muller.text(-0.1,1.05,'a   '+config.pop_labels[population], fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)
		# ax_muller.text(0.9,0.93,'YPD', fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)

		ax_muller.text(-0.1,-0.20,'c', fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)
		ax_muller.set_ylabel('Lineage')
		ax_submuller[0].set_ylabel('Sublineage')
		ax_submuller[2].set_ylabel('Sublineage')
	else:
		ax_muller.text(-0.10,1.05,'b   '+config.pop_labels[population], fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)
		ax_muller.text(-0.10,-0.20,'d', fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)
		# ax_muller.text(0.9,0.93,'YPA', fontsize = 8, weight = 'bold',transform= ax_muller.transAxes)

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)
	for epoch in range(1,max_barcode):
		times[epoch*11:]+= additional_length_of_barcoding

	mutant_tree = population_tree

	muller_diagram = MullerDiagram(fig_muller,ax_muller,times)

	sys.stderr.write("\n file parsing and preprocessing done \n\n")

	max_freqs = compile_max_freqs(population_tree.keys(),lineage_dict)
	sorted_ids = [x for y,x in reversed(sorted(zip(max_freqs,population_tree.keys())))]

	assign_muller_lower(sorted_ids,lineage_dict,population_tree,
		baseline = numpy.zeros(len(times)),
			range_available = numpy.ones(len(times)))

	with_stars = True

	for ID in clone_list:
		children = child_list(ID, population_tree)

		if ID in population_tree.keys():
			outline = False
		else:
			outline = True
		muller_diagram.place(lineage_dict[ID],outline = outline)
	
		if with_stars:
			# place star at likely time of mutation
			barcode_epoch = len(ID.split("_"))
			pos = (barcode_epoch-1)*11 + 5
			t = (110+additional_length_of_barcoding)*(barcode_epoch-1)+50

			y = lineage_dict[ID].muller_lower[pos] + lineage_dict[ID].freqs[pos]/2.

			star_kwargs = {'marker':(5,1),'s':8,'lw':0.1,'zorder':20,'edgecolors' :'0.2'}

			parents = ancestor_list(ID, population_tree)
			for p_ID in parents:
				ax_muller.scatter(t,y,color = clone_dict[p_ID].color,**star_kwargs)

			ax_muller.scatter(t,y,color = clone_dict[ID].color,**star_kwargs)

	for epoch in range(1,max_barcode):
		rectangle_right = (110 + additional_length_of_barcoding)*epoch
		if stretched:
			bar_alpha = 0.3
		else:
			bar_alpha = 0.7
		ax_muller.add_patch(patches.Rectangle(
        	(rectangle_right - length_of_barcoding, -0.1), 10 + additional_length_of_barcoding, 1.1,
        	facecolor= '0.3',
        	alpha = bar_alpha,
        	edgecolor = "none"))

	lineageindex = -1
	for parent_ID in submuller_lineages[population]:
		lineageindex += 1
		offset = deepcopy(lineage_dict[parent_ID].muller_lower)
		scale = numpy.ones(len(lineage_dict[parent_ID].freqs))
		scale[lineage_dict[parent_ID].freqs > 0.] = lineage_dict[parent_ID].freqs[lineage_dict[parent_ID].freqs > 0.]
		scale = deepcopy(scale)
		
		submuller_index_map = {i:i for i in range(len(ax_submuller))}
		submuller_index = submuller_index_map[lineageindex]
		submuller_diagram = MullerDiagram(fig_muller,ax_submuller[submuller_index],times)

		for ID in clone_list:

			if ID.startswith(parent_ID):

				lineage_dict[ID].muller_lower = (lineage_dict[ID].muller_lower - offset)/scale
				# ax_submuller.plot(times, (lineage_dict[ID].muller_lower - offset)/scale, color = lineage_dict[ID].color,zorder = 20)
				# ax_submuller.plot(times, (lineage_dict[ID].muller_lower - offset)/scale + lineage_dict[ID].freqs/scale, color = lineage_dict[ID].color,zorder = 20)
				submuller_diagram.place(lineage_dict[ID],scale = scale,outline = True)	

				if with_stars:
					# place star at likely time of mutation
					barcode_epoch = len(ID.split("_"))

					if barcode_epoch > len(parent_ID.split("_")): # no star for focal mutation

						pos = 11*(barcode_epoch-1)+5
						t = (110+additional_length_of_barcoding)*(barcode_epoch-1)+50

						y = lineage_dict[ID].muller_lower[pos] + lineage_dict[ID].freqs[pos]/2./scale[pos]

						star_kwargs = {'marker':(5,1),'s':8,'lw':0.1,'zorder':20,'edgecolors': '0.2'}

						parents = ancestor_list(ID, population_tree)
						for p_ID in parents:
							ax_submuller[submuller_index].scatter(t,y,color = clone_dict[p_ID].color,**star_kwargs)

						ax_submuller[submuller_index].scatter(t,y,color = clone_dict[ID].color,**star_kwargs)
					lineage_dict[ID].muller_lower *= scale
					lineage_dict[ID].muller_lower += offset
					xtick_dt = 110+additional_length_of_barcoding
					ax_submuller[submuller_index].set_xticks([xtick_dt*i for i in range(0,max_barcode+1)][::2])
					ax_submuller[submuller_index].set_xticklabels(['%d.0'%(2*i+1)for i in range(0,(max_barcode+1)/2)])
					ax_submuller[submuller_index].set_xlim([0,xtick_dt*max_barcode-length_of_barcoding])	

					
					if submuller_index < 2:
						ax_submuller[submuller_index].set_xticklabels([])
	xticks = []
	xticklabels = []
	xtick_dt = 110+additional_length_of_barcoding
	for i in range(0,max_barcode+1):
		if i == 0:
			xticks.extend([xtick_dt*i,xtick_dt*i+100])
		else:
			xticks.extend([xtick_dt*i,xtick_dt*i+100])
		xticklabels.extend(['%d.0'%(i+1),'%d.100'%(i+1)])
	ax_muller.set_xticks(xticks)
	ax_muller.set_xticklabels(xticklabels)
	ax_muller.set_xlim([0,xtick_dt*max_barcode-length_of_barcoding])

	i = 0
	for tick in ax_muller.get_xticklabels():
   		tick.set_rotation(90)
   		if stretched:
   			tick.set_horizontalalignment('center')
   		else:
	   		if i % 2 == 0:
	   			tick.set_horizontalalignment('left')
	   		else:
	   			tick.set_horizontalalignment('right')
	   		i+=1
   	muller_diagram.save(OUTFILE,format = 'pdf')			

