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

fig, ax = plt.subplots(2,1,figsize = (3,2))
fig_p_vs_f, ax_p_vs_f = plt.subplots(figsize = (2,2))



models = ['frequency','frequency_plus_fitness']
for model in models:

	fig_all, ax_all = plt.subplots(1,2,figsize= (6,1.5))
	grid = gridspec.GridSpec(1,2, width_ratios = [72,124],hspace = 0.06,wspace = 0.1)

	ax_all[0] = plt.subplot(grid[0])
	ax_all[1] = plt.subplot(grid[1])

	for axes in [ax_all[0],ax_all[1]]:
		axes.patch.set_alpha(0.)
		for axis_position in ['top','right']:
			axes.axes.spines[axis_position].set_visible(False)
		axes.axes.xaxis.set_ticks_position('bottom')
		axes.axes.yaxis.set_ticks_position('left')
	# axes.set_yticklabels([])
	# axes.set_xticklabels([])

	popindex = -1
	for population in config.populations:
		popindex += 1

		max_barcode = config.max_barcode[population]


		clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

		clone_list.append("")

		clone_dict[''].color = 'k'

		
		p_values_more_dict = { model:{} for model in models }
		p_values_less_dict = { model:{} for model in models }

		filename = config.clone_data_directory+'%s-p_values_segment_%s.tsv' % (population,model)
		f = open(filename,'r')
		#skip header
		f.readline()
		for row in csv.reader(f,delimiter='\t'):
			bcd = row[0]
			
			p_s = numpy.asarray(row[1:],dtype = float)
			p_values_more_dict[model].update({bcd:p_s[:max_barcode-1]})
			p_values_less_dict[model].update({bcd:p_s[max_barcode-1:]})
		p_values_more = { model:numpy.zeros(len(clone_list)*(max_barcode-1)) for model in models}
		p_values_less = { model:numpy.zeros(len(clone_list)*(max_barcode-1)) for model in models}

		symbols = {'frequency':'o','frequency_plus_fitness':'o'}
		model_labels = {'frequency':'frequency only','frequency_plus_fitness':'frequency and fitness'}

		for i in range(0,len(clone_list)):
			first = i*(max_barcode-1)
			last = (i+1)*(max_barcode-1)
			p_values_more[model][first:last] = p_values_more_dict[model][clone_list[i]]
			p_values_less[model][first:last] = p_values_less_dict[model][clone_list[i]]


		p_values = p_values_more[model]
		FDR = 0.05
		candidate_pvalues = sorted(np.unique(p_values[p_values<1.-10**-4]),reverse=True)
		num = len(p_values[p_values < 1.-10**-4])
		print num
		for p in candidate_pvalues:
			if (num*p/(p_values <= p).sum()) < FDR:
				print FDR, "FDR threshold:", p
				break

		
		for rank in range(0,len(clone_list)):
			ID = clone_list[rank]

			# record which epochs lineage is seen in
			seen = np.zeros((max_barcode-1))
			for epoch in range(3,max_barcode-1):
				begin= inference_params.INTERVALS_PER_EPOCH*epoch
				end = inference_params.INTERVALS_PER_EPOCH*(epoch+1)
				if numpy.sum(clone_dict[ID].freqs[begin:end]) > 10**-5:
					seen[epoch] = 1
			#plot p_values for these epochs
			xs = np.asarray([rank+1]*(max_barcode-1))
			ys = -np.log10(p_values_more_dict[model][ID]+10**-5)

			xs = xs[seen==1]
			ys = ys[seen==1]

			ax_all[popindex].scatter(xs, ys,
							marker = symbols[model],s = 8, color = clone_dict[ID].color,lw = 0,clip_on = False)
		print ""
		ax_all[popindex].scatter([], [], marker = symbols[model],s = 8, color = 'k',lw = 0,clip_on = True,label = model_labels[model])

		ax_all[popindex].set_ylim([-.1,5])
		ax_all[popindex].set_xlim([0,len(clone_list)])
		if (num*p/(p_values <= p).sum()) < FDR:
			line1, = ax_all[popindex].plot([0,len(clone_list)],[-np.log10(p+10**-5),-np.log10(p+10**-5)],lw = 0.5, color= 'k')

			line1, = ax_all[popindex].plot([],[],lw = 0.5, color= 'k',label = r'FDR < %d%%,%s' % (100*FDR,model_labels[model]))

		ax_p_vs_f.set_xscale('log')
		ax_p_vs_f.set_yscale('log')

	ax_all[1].legend(loc = 'best',fontsize = 6,scatterpoints = 1,handlelength = 1.)
	ax_all[0].set_xlabel('Clone rank (YPD)',fontsize = 8)
	ax_all[1].set_xlabel('Clone rank (YPA)',fontsize = 8)
	ax_all[0].set_ylabel(r'$-\log10 P_{m_{i,E,\mathrm{model}} \geq m_{i,E,\mathrm{exp}}}$',fontsize = 8)
	for axis in [ax_all[0],ax_all[1]]:
		axis.yaxis.set_ticks([0,1,2,3,4,5])
		axis.yaxis.set_ticklabels([0,1,2,3,4,'>4'])

	fig_all.savefig(config.figure_directory + "si/pvalues_more_all_%s.pdf" % model, bbox_inches = 'tight')
