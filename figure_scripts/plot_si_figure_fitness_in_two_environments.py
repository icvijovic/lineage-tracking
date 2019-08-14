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
from matplotlib import rc


from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import lineage.file_parser as file_parser

from lineage.read_clone_data import *
from lineage.plot_utils_clone import * # contains muller and bar chart class
from lineage.tree_utils import *

matplotlib.rcParams['axes.linewidth'] = .25

parser = argparse.ArgumentParser()
parser.add_argument("-dpi", default = 300, type = int, help = 'resolution (dots per inch)')
parser.add_argument("-width", default = 7.2, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 5.76, type = float, help = 'figure height (inches)')
parser.add_argument("-fontsize", default = 6, type = int, help = 'fontsize (pt)')
parser.add_argument("-outfile", default = 'joint_fitness_distribution.pdf', type = str, help = 'figure name')


args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize
OUTFILE = args.outfile

matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize + 1


figure_directory = config.figure_directory

OUTFILE = figure_directory + 'si/' + OUTFILE

fig,ax = plt.subplots(1,2,figsize = (4,2))
outer_grid = gridspec.GridSpec(1,2,hspace = 0.25,wspace = 0.1)

text = {'C1':'a    YPD', 'D1':'b    YPA'}
popindices = {'C1':0,'D1':1}
xlim = {'C1':[-1,5],'D1':[-2,10]}

for population in config.populations:
	popindex = popindices[population]
	

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)
	
	# plot ancestor in black
	clone_dict[''].color = 'k'

	for bcd, clone in clone_dict.iteritems():
		x = clone.evolution_fitness*100

		x_l,x_u = 100*clone.evolution_fitness_CI[0], 100*clone.evolution_fitness_CI[1]

		y = (x + clone.barcoding_fitness)/2.
		sigma_b = (clone.barcoding_fitness_CI[1] - clone.barcoding_fitness_CI[0])/2.
		sigma_e = (x_u - x_l)/2.
		sigma_average = np.sqrt(sigma_b**2 + sigma_e**2)/2.
		y_l, y_u = y - sigma_average, y + sigma_average

		kwargs = dict(color = clone.color,s = 10,marker = 'o',lw = 0.)

		ax[popindex].scatter(x,y,**kwargs)

		ax[popindex].plot([x_l,x_u],[y,y],color = clone.color,lw = 0.5)
		ax[popindex].plot([x,x],[y_l,y_u],color = clone.color,lw = 0.5)


	ax[popindex].plot(xlim[population],xlim[population],lw = 1, color = 'k')
	ax[popindex].plot(xlim[population],[0.5*l for l in xlim[population]],lw = 1,color = 'k',alpha = 0.5)
	ax[popindex].set_xlim(xlim[population])
	ax[popindex].set_ylim(xlim[population])

	ax[popindex].set_xlabel('Evolution fitness (per epoch)',fontsize = 8)
	ax[popindex].text(-0.1,1.1,config.pop_labels[population], fontsize = 8, weight = 'bold',transform=ax[popindex].transAxes)
ax[0].set_ylabel('Average fitness (per epoch)',fontsize = 8)
fig.savefig(OUTFILE,dpi = dpi,bbox_inches = 'tight')			


