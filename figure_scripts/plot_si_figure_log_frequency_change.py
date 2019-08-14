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

fontsize = 8

# set up figures and axes
fig,ax = plt.subplots(2,3,figsize = (8,4),sharey = True)
grid = gridspec.GridSpec(2,3, width_ratios = [72,124,30],hspace = 0.25,wspace = 0.04)

ax_bc = [plt.subplot(grid[0,0]), plt.subplot(grid[0,1])]
ax_ev = [plt.subplot(grid[1,0]), plt.subplot(grid[1,1])]
ax_hist_bc = plt.subplot(grid[0,2])
ax_hist_ev = plt.subplot(grid[1,2])

for axes in [ax_bc[0],ax_bc[1],ax_ev[0],ax_ev[1],ax_hist_bc,ax_hist_ev]:
	axes.patch.set_alpha(0.)
	for axis_position in ['top','right']:
		axes.axes.spines[axis_position].set_visible(False)
	axes.axes.xaxis.set_ticks_position('bottom')
	axes.axes.yaxis.set_ticks_position('left')
	# axes.set_yticklabels([])
	# axes.set_xticklabels([])

for environment, ax, ax_hist in zip(['Evolution', 'Barcoding'],[ax_ev, ax_bc],[ax_hist_ev, ax_hist_bc]):
	ax[1].set_yticklabels([])
	ax_hist.set_yticklabels([])

	popindex = -1
	all_log_changes = []

	ax[0].text(-0.15,1.05,'%s environment' % environment, weight = 'bold', fontsize = 8, transform = ax[0].transAxes)
	for population in config.populations:
		popindex += 1
		print population

		clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

		max_barcode = config.max_barcode[population]


		rank = 0
		num_total = 0
		num_g_1 = 0
		num_g_2 = 0
		num_g_3 = 0
		print len(clone_list)
		for ID in clone_list:
			rank +=1 

			first_epoch = len(ID.split("_")) - 1

			log_changes = []
			relative_fitnesses = []
			for epoch in range(first_epoch,max_barcode-(environment == 'Barcoding')):
				if environment == 'Evolution':
					timepoint_before = epoch * 11
					timepoint_after = epoch * 10
				else:
					timepoint_before=epoch * 11+10
					timepoint_after = epoch * 11+11

				f_b = clone_dict[ID].freqs[timepoint_before]
				f_a = clone_dict[ID].freqs[timepoint_after]

				epsilon = 5*10**-4
				log_freq_change = numpy.log(f_a+epsilon) - numpy.log(f_b+epsilon)
				log_changes.append(log_freq_change)


			log_changes = numpy.asarray(log_changes)
			num_total += sum(1 for i in log_changes)
			num_g_1 += sum(1 for i in log_changes if abs(i) > 1)
			num_g_2 += sum(1 for i in log_changes if abs(i) > 2)
			num_g_3 += sum(1 for i in log_changes if abs(i) > 3)

			ax[popindex].scatter(rank*numpy.ones(len(log_changes)), log_changes,marker = 's',s = 8, color = clone_dict[ID].color,lw = 0,clip_on = False)
			# ax[popindex].scatter(rank, mean_rel_fit, marker  = 'o', color = clone_dict[ID].color,lw = 0)
			all_log_changes.extend(log_changes)
		ax[popindex].set_ylim([-4,4])

		ax[popindex].plot([0,len(clone_list)],[0,0],lw = 1, color = 'k')
		w = 2
		ax[popindex].fill_between([0,len(clone_list)],[-w,-w],[w,w],lw = 0, color = 'k',alpha = 0.1,zorder = -1)
		w = 1
		ax[popindex].fill_between([0,len(clone_list)],[-w,-w],[w,w],lw = 0, color = 'k',alpha = 0.1,zorder = -1)
		ax[popindex].set_xlim([0,len(clone_list)])

		ax_hist.set_ylim([-4,4])
		bin_edges = np.arange(-4,5)
		bins = bin_edges
		
		print num_total
		print num_g_1, " are greater than 1, (", 1.- num_g_1*1./num_total,")"
		print num_g_2, " are greater than 2, (",1. - num_g_2*1./num_total,")"
		print num_g_3, " are greater than 3, (", 1.-num_g_3*1./num_total,")"
	n,x,y = ax_hist.hist(all_log_changes,bins = bins,color = 'k',alpha = 0.4,lw = 0,rwidth = 0.9,orientation='horizontal',normed = True)
	ax_hist.set_xlim([0,1])
	print n, x
	if environment == 'Barcoding':
		ax[0].set_ylabel('Log frequency change\n'+r'during barcoding', fontsize = fontsize)
	else:
		ax[0].set_ylabel('Log frequency change\n'+r'over 100 generations', fontsize = fontsize)
	if environment == 'Evolution':
		ax[0].set_xlabel('Clone rank (YPD)', fontsize = fontsize)
		ax[1].set_xlabel('Clone rank (YPA)', fontsize = fontsize)
		ax_hist.set_xlabel('Percent', fontsize = fontsize)
	# ax.legend(bbox_to_anchor = (1.5,0.2),scatterpoints = 1,labelspacing = .8,handlelength = 1)
	fig.savefig(config.figure_directory+'si/log_freq_change.pdf',bbox_inches = 'tight')



		

			