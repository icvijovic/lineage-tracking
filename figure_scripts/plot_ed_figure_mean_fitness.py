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
from scipy.ndimage.filters import gaussian_filter


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

fig, ax = plt.subplots(1,1,figsize = (4,2))
popindex = -1


fontsize = 7

popindex = -1
for population in config.populations:
	popindex += 1

	max_barcode = config.max_barcode[population]

	fitness_assay_file = glob.glob('/Users/IvanaCvijovic/Dropbox/InfiniteLineageTracking/Fitness_measurements/fitnesses_only_'+population+'*csv')

	for filename in fitness_assay_file:
		f = open(filename,'r')
		fitness_dict = {}
		f.readline()
		delimiter = '\t'
		for row in csv.reader(f,delimiter= delimiter):
			fitness_dict[int(row[-2])] = float(row[-1])
		f.close()
		fitness_assay_times = list(itertools.chain.from_iterable((110*i, 110*i + 100) for i in xrange(0, max_barcode)))
		# fitness_assay_times = [0,100,110,210,220,320,330,430,440,540,550,650,660,760,770,870,880,980,990,1090]
		mean_fitnesses = numpy.asarray([fitness_dict[key] for key in sorted(fitness_dict.keys())])
		mean_fitnesses = mean_fitnesses[:2*max_barcode]
		print len(mean_fitnesses), len(fitness_assay_times)

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population, assign_colors = False)

	single_inferred_mean_fitness = numpy.zeros(len(times))
	for ID in clone_list:
		single_inferred_mean_fitness += clone_dict[ID].freqs*clone_dict[ID].evolution_fitness

	single_inferred_mean_fitness = single_inferred_mean_fitness[:11*max_barcode]
	times = times[:11*max_barcode]

	plt.plot(fitness_assay_times,mean_fitnesses,'o',color = config.pop_colors[population],clip_on = False)

	offset = - single_inferred_mean_fitness[1*33]+mean_fitnesses[1*5]
	
	inf_xbar_upper = numpy.zeros(110)
	inf_xbar_lower = numpy.zeros(110)

	clone_list.append("")
	for ID in clone_list:
		inf_xbar_upper += clone_dict[ID].freqs * clone_dict[ID].evolution_fitness_CI[1]
		inf_xbar_lower += clone_dict[ID].freqs * clone_dict[ID].evolution_fitness_CI[0]

	for epoch in range(0,max_barcode):
		first = 11*epoch
		last = 11*epoch + 10

		plt.plot(times[first:last],single_inferred_mean_fitness[first:last]+offset,lw = 1,color = config.pop_colors[population])
		plt.fill_between(times[first:last],inf_xbar_lower[first:last] +offset, inf_xbar_upper[first:last] + offset,color =config.pop_colors[population],alpha = 0.5 )

	line1, = plt.plot(times,single_inferred_mean_fitness+offset,lw =1, color = config.pop_colors[population])
	line1.set_dashes((1,1))
	ax.set_xlim([0,110*max_barcode])

	offset = - single_inferred_mean_fitness[popindex*(-1)]+mean_fitnesses[popindex*(-1)]

	if population == 'D1':
		for epoch in range(max_barcode-2,max_barcode):
			first = 11*epoch
			last = 11*epoch + 10

			plt.plot(times[first:last],single_inferred_mean_fitness[first:last]+offset,lw = 1,color = config.pop_colors[population],alpha = 0.5)
			plt.fill_between(times[first:last],inf_xbar_lower[first:last] +offset, inf_xbar_upper[first:last] + offset,color =config.pop_colors[population],alpha = 0.2 )


ax.plot([],[],'o',color ='k',label = 'FACS fitness assay')
ax.plot([],[],lw = 1, color = 'k',label = 'inferred')

ax.plot([],[],lw = 5, color = config.pop_colors['C1'],label = 'YPD')
ax.plot([],[],lw = 5,color = config.pop_colors['D1'],label = 'YPA')

ax.set_xlabel(r'Generations, $t$',fontsize = fontsize)
ax.set_ylabel('Mean evolution\nenvironment fitness,'+r' $\bar{x}(t)$',fontsize = fontsize)
ax.legend(loc = 'best',fontsize = fontsize,ncol = 1,handlelength = 1.)

xticks = []
xticklabels = []
xtick_dt = 110
for i in range(0,max_barcode+1):
	if i == 0:
		xticks.extend([xtick_dt*i,xtick_dt*i+100])
	else:
		xticks.extend([xtick_dt*i,xtick_dt*i+100])
	xticklabels.extend(['%d.0'%(i+1),'%d.100'%(i+1)])
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
ax.set_xlim([0,xtick_dt*max_barcode-10])

i = 0
for tick in ax.get_xticklabels():
		tick.set_rotation(90)
   		if i % 2 == 0:
   			tick.set_horizontalalignment('left')
   		else:
   			tick.set_horizontalalignment('right')
   		i+=1

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
# FIGNAME = 'master_mullers/%s_avg_fitness.pdf' % population
# ax.text(-0.15,1.1,'Extended Data Figure 2',fontsize = 8,weight = 'bold',transform=ax.transAxes)
FIGNAME = config.figure_directory+'ed/EDFigure2_final.pdf'
sys.stderr.write('Saving to %s...\n'%FIGNAME)
fig.savefig(FIGNAME,bbox_inches = 'tight')
