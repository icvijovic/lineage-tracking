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

## add custom modules to path
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
parser.add_argument("-width", default = 2, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 2, type = float, help = 'figure height (inches)')
parser.add_argument("-fontsize", default = 6, type = int, help = 'fontsize (pt)')
parser.add_argument("-outfile", default = 'temp.png', type = str, help = 'figure name')


args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize
OUTFILE = args.outfile

matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize + 1


def place_children(ID,new_list,population_tree):
	for child_ID in child_list(ID,population_tree):
		new_list.append(child_ID)
		place_children(child_ID,new_list,population_tree)

for population in config.populations:

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

	sorted_clone_list = []
	for ID in clone_list:
		if ID in sorted_clone_list:
			pass
		else:
			sorted_clone_list.append(ID)
			place_children(ID,sorted_clone_list,population_tree)


	# hide axes
	
	# columns = ['ID', 'Fitness (percent)', 'Confidence interval (percent)']
	index = 0
	max_lines_per_page = 42
	num_lists = len(sorted_clone_list)/max_lines_per_page + 1
	if len(sorted_clone_list)<=max_lines_per_page:
		sorted_clone_lists = [sorted_clone_list]
	else:
		sorted_clone_lists = [sorted_clone_list[i:i+max_lines_per_page] for i in range(0, len(sorted_clone_list), max_lines_per_page)]

	for sorted_clone_list in sorted_clone_lists:
		index += 1
		fig,ax=plt.subplots(figsize = (6,0.02))
		plt.axis('off')

		table_text = []
		colors = []
		for ID in sorted_clone_list:
			fit = clone_dict[ID].evolution_fitness
			CI_lower, CI_upper = clone_dict[ID].evolution_fitness_CI
			bc_fit = clone_dict[ID].barcoding_fitness
			bc_CI_lower, bc_CI_upper = clone_dict[ID].barcoding_fitness_CI

			ID_text = ID.split("_")
			if len(ID_text) <= 5:
				ID_text = "_".join(ID_text)
			else:
				ID_text = "\n".join(["_".join(ID_text[0:5]),"_"+"_".join(ID_text[5:])])
			row_text = [ID_text,r'%.2f'%(100*fit), r'(%.2f,%.2f)'%(100*CI_lower,100*CI_upper)]
			row_text.extend([r'%.2f'%(bc_fit), r'(%.2f,%.2f)'%(bc_CI_lower,bc_CI_upper), ''])
			w = (1.,1.,1.,1.)
			colors.append([w,w,w,w,w,clone_dict[ID].color])
			table_text.append(row_text)

		first_cell_width = 0.75*1.2
		# if population == 'D1':
		# 	first_cell_width *= 1.25
		columnNames = ['Clone Barcodes','Evolution\nFitness\n(percent)',r'95% CI'+'\n(per cycle)','Barcoding\nFitness\n(percent)',r'95% CI'+'\n(percent)','Color']
		table = plt.table(cellText=table_text,cellColours = colors, cellLoc = 'left',colWidths=[first_cell_width,0.08,0.1,0.08,0.1,0.05],colLabels = columnNames)

		for key, cell in table.get_celld().items():
			cell.set_linewidth(0.1)
			
		def set_pad_for_column(col, pad=0.1):
		    cells = [key for key in table._cells if key[1] == col]
		    for cell in cells:
		        table._cells[cell].PAD = pad

		set_pad_for_column(col = 0,pad = 0.01)
		set_pad_for_column(col = 1, pad = 0.1)
		set_pad_for_column(col = 3, pad = 0.1)

		cellDict = table.get_celld()
		for i in range(0,len(columnNames)):
		    cellDict[(0,i)].set_height(30)

		table.auto_set_font_size(False) 
		table.set_fontsize(5)

		bbox = table.get_clip_box()
		# add 10 pixel spacing
		# points[0,:] -= 10; points[1,:] += 10


		plt.savefig(config.figure_directory+'si/%s_table_%d.pdf'%(population,index),bbox_inches = 'tight')
