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


fig_p_vs_f, ax_p_vs_f = plt.subplots(figsize = (2,2))

popindex = -1
for population in config.populations:
	popindex += 1

	max_barcode = config.max_barcode[population]


	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

	clone_list.append("")

	clone_dict[''].color = 'k'

	models = ['frequency','frequency_plus_fitness']
	p_values_more_dict = { model:{} for model in models }
	p_values_less_dict = { model:{} for model in models }

	for model in models:
		filename = config.clone_data_directory+'%s-p_values_%s.tsv' % (population,model)
		f = open(filename,'r')
		#skip header
		f.readline()
		for row in csv.reader(f,delimiter='\t'):
			bcd = row[0]
			
			p_s = numpy.asarray(row[1:],dtype = float)
			p_values_more_dict[model].update({bcd:p_s[0]})
			p_values_less_dict[model].update({bcd:p_s[1]})

	p_values_more = { model:numpy.zeros(len(clone_list)) for model in models}
	p_values_less = { model:numpy.zeros(len(clone_list)) for model in models}

	symbols = {'frequency':'.','frequency_plus_fitness':'o'}
	linestyle = {'frequency':':','frequency_plus_fitness':'-'}
	model_labels = {'frequency':'frequency only','frequency_plus_fitness':'frequency and fitness'}

	for model in models:
		for i in range(0,len(clone_list)):
			p_values_more[model][i] = p_values_more_dict[model][clone_list[i]]
			p_values_less[model][i] = p_values_less_dict[model][clone_list[i]]


		p_values = p_values_more[model]
		FDR = 0.1
		candidate_pvalues = sorted(np.unique(p_values[p_values<1.-10**-4]),reverse=True)
		num = len(p_values[p_values < 1.-10**-4])
		for p in candidate_pvalues:
			if (num*p/(p_values <= p).sum()) < FDR:
				print FDR, "FDR threshold:", p
				break

		
		for rank in range(0,len(clone_list)):
			ID = clone_list[rank]
			ax_all[popindex].scatter(rank+1, -np.log10(p_values[rank]+10**-5),
							marker = symbols[model],s = 8, color = clone_dict[ID].color,lw = 0,clip_on = False)
		ax_all[popindex].scatter([], [], marker = symbols[model],s = 8, color = 'k',lw = 0,clip_on = True,label = model_labels[model])
		# ax_all[popindex].set_yscale('log')
		ax_all[popindex].set_ylim([-.1,5])
		ax_all[popindex].set_xlim([0,len(clone_list)])
		if (num*p/(p_values <= p).sum()) < FDR:
			ax_all[popindex].plot([0,len(clone_list)],[-np.log10(p+10**-5),-np.log10(p+10**-5)],lw = 0.5, color= 'k',linestyle=linestyle[model])
			ax_all[popindex].plot([],[],lw = 0.5, color= 'k',linestyle=linestyle[model],label = model_labels[model])
		# ax_p_vs_f.plot(average_frequency,p_values,'.',color = config.pop_colors[population])
		ax_p_vs_f.set_xscale('log')
		ax_p_vs_f.set_yscale('log')

ax_all[1].legend(loc = 'best',fontsize = 6,scatterpoints = 1,handlelength = 1.)
ax_all[0].set_xlabel('Clone rank (YPD)',fontsize = 8)
ax_all[1].set_xlabel('Clone rank (YPA)',fontsize = 8)
ax_all[0].set_ylabel(r'$-\log10 P_{m_{i,\mathrm{model}} \geq m_{i,\mathrm{exp}}}$',fontsize = 8)
for axis in [ax_all[0],ax_all[1]]:
	axis.yaxis.set_ticks([0,1,2,3,4,5])
	axis.yaxis.set_ticklabels([0,1,2,3,4,'>4'])

ax_all[1].legend(loc = 'best',fontsize = 7,scatterpoints = 1)
# fig.savefig(config.figure_directory + "pvalues.pdf",bbox_inches = 'tight')
fig_all.savefig(config.figure_directory + "si/pvalues_clone_more_all.pdf", bbox_inches = 'tight')
# fig_p_vs_f.savefig(config.figure_directory + "pvsf.pdf",bbox_inches = 'tight')