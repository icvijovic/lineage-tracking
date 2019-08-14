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
parser.add_argument("-outfile", default = 'Figure3_final.pdf', type = str, help = 'figure name')
parser.add_argument('--evolution', dest='choice', action='store_const',
                    const='evolution', default='average',
                    help='choice of fitness to use (default: evolution)')
parser.add_argument('--barcoding', dest='choice', action='store_const',
                    const='barcoding', default='average',
                    help='choice of fitness to use (default: evolution)')
parser.add_argument('--average', dest='choice', action='store_const',
                    const='average', default='average',
                    help='choice of fitness to use (default: evolution)')

args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize
OUTFILE = args.outfile

matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize + 1

choice = args.choice

if choice == 'average':
	subfolder = 'main/'
else:
	subfolder = 'si/'

if choice == 'barcoding':
	unit = 'per interval'
else:
	unit = 'percent'

figure_directory = config.figure_directory + subfolder

OUTFILE = figure_directory + OUTFILE.split(".")[0] + "_" + choice + ".pdf"

fig = plt.figure(figsize = (fig_width,fig_height))
outer_grid = gridspec.GridSpec(1,1,hspace = 0.25,wspace = 0.1)

text_panels_AB = {'C1':'a    YPD', 'D1':'b    YPA'}
popindices = {'C1':0,'D1':1}

sys.stderr.write("Plotting traveling waves using %s fitness.\n" % choice)

for population in config.populations:
	popindex = popindices[population]
	
	max_barcode = config.max_barcode[population]


	height = 0.3 # specifies height of longest bar in terms of total figure height6
	left = 0.04 + (1-popindex)*0.08
	YPD_bottom = 0.57
	bottom = YPD_bottom - popindex* (YPD_bottom - 0.05)
	move_bottom = -0.025 
	move_left = 0.05*4.5/7.5
	if population == 'C1':

		min_fitness = -0.01
		if choice == 'average':
			max_fitness = 0.04
		else:
			max_fitness = 0.05
		width = 0.83*0.1*4.5/7.5
	else:
		width = 0.83*0.1*4.5/7.5
		if choice == 'average':
			max_fitness = 0.08
			min_fitness = -0.01
		else:
			max_fitness = 0.1
			min_fitness = -0.02

	fitness_range = max_fitness - min_fitness
	width_multiplier = 2*fitness_range/0.06
	if choice == 'average':
		width_multiplier *= 1.5

	ax_fitness = []

	for i in range(0,max_barcode+1):
		new_axes = fig.add_axes([left,bottom,width*width_multiplier,height])
		ax_fitness.append(new_axes)
		if i > 0:
			ax_fitness[i].set_zorder(ax_fitness[i-1].get_zorder()-1) # put ax[i+1] behind ax[i-1] 
		ax_fitness[i].patch.set_alpha(0.) # make axis transparent
		ax_fitness[i].set_ylim([0,1])
		# ax_fitness[i].set_position([left,bottom,width*width_multiplier,height])
		left += move_left
		if i < max_barcode - 1:
			bottom -= move_bottom
		else:
			bottom -= move_bottom*inference_params.EVOLUTION_INTERVALS_PER_EPOCH/inference_params.INTERVALS_PER_EPOCH
	# if popindex == 0:
		# ax_fitness[0].text(-.35*0.1/fitness_range,1.3,'Figure 3',transform = ax_fitness[0].transAxes,weight = 'bold',fontsize = 8)

	ax_fitness[0].text(-.05*0.1/fitness_range,1.1,text_panels_AB[population],transform = ax_fitness[0].transAxes,weight = 'bold',fontsize = 8)
	ax_fitness[0].set_xlabel('Fitness relative\nto ancestor (%s)'% unit,labelpad = 0)
	ax_fitness[0].set_ylabel('Frequency',labelpad = 1)


	mean_fitness_width = width*width_multiplier+(len(ax_fitness) - 1)*move_left
	mean_fitness_height = -move_bottom * (len(ax_fitness) - 2
										 + 1.*inference_params.EVOLUTION_INTERVALS_PER_EPOCH/inference_params.INTERVALS_PER_EPOCH)

	mean_fitness_left, mean_fitness_bottom, x, y = ax_fitness[0].get_position().bounds

	ax_mean_fitness = fig.add_axes([mean_fitness_left,mean_fitness_bottom,mean_fitness_width,mean_fitness_height])
	ax_mean_fitness.set_zorder(ax_fitness[-1].get_zorder()-1)
	ax_mean_fitness.patch.set_alpha(0.)

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

	single_inferred_mean_fitness_evolution = numpy.zeros(len(times))
	single_inferred_mean_fitness_barcoding = numpy.zeros(len(times))

	for bcd, clone in clone_dict.iteritems():
		single_inferred_mean_fitness_evolution += clone.freqs*clone.evolution_fitness
		single_inferred_mean_fitness_barcoding += clone.freqs*clone.barcoding_fitness

	sequencing_times = times[:11*max_barcode]
	single_inferred_mean_fitness_evolution = single_inferred_mean_fitness_evolution[:11*max_barcode]
	single_inferred_mean_fitness_barcoding = single_inferred_mean_fitness_barcoding[:11*max_barcode]

	if choice == 'evolution':	
		single_inferred_mean_fitness = single_inferred_mean_fitness_evolution
		for bcd,clone in clone_dict.iteritems():
			clone.relative_fitness = numpy.ones(len(clone.relative_fitness))*clone.evolution_fitness
	elif choice == 'barcoding':
		single_inferred_mean_fitness = single_inferred_mean_fitness_barcoding/100.
		for bcd,clone in clone_dict.iteritems():
			clone.relative_fitness = numpy.ones(len(clone.relative_fitness))*clone.barcoding_fitness/100.
	elif choice == 'average':
		single_inferred_mean_fitness = 0.5*(single_inferred_mean_fitness_evolution + single_inferred_mean_fitness_barcoding/100.)
		for bcd,clone in clone_dict.iteritems():
			clone.relative_fitness = numpy.ones(len(clone.relative_fitness))*0.5*(clone.evolution_fitness + clone.barcoding_fitness/100.)

	ax_mean_fitness.set_xlim([min_fitness,max_fitness])
	ax_mean_fitness.set_ylim([min(sequencing_times),max(sequencing_times)])

	x_range = ax_mean_fitness.get_xlim()[1] - ax_mean_fitness.get_xlim()[0]
	y_range = ax_mean_fitness.get_ylim()[1] - ax_mean_fitness.get_ylim()[0]

	width_scale = width *width_multiplier / (width*width_multiplier+(len(ax_fitness)-1)*move_left)

	skewed_fitness = numpy.zeros(len(sequencing_times))
	skewed_fitness += min_fitness 
	skewed_fitness -= min_fitness*width_scale 
	skewed_fitness += (sequencing_times)/y_range*(mean_fitness_width-width*width_multiplier)/mean_fitness_width*x_range
	skewed_fitness += single_inferred_mean_fitness*width_scale
	

	skewed_left_axis = numpy.zeros(len(sequencing_times))
	skewed_left_axis += min_fitness 
	skewed_left_axis -= min_fitness*width_scale 
	skewed_left_axis += (sequencing_times)/y_range*(mean_fitness_width-width*width_multiplier)/mean_fitness_width*x_range
	skewed_right_axis = skewed_left_axis + max_fitness*width_scale
	ax_mean_fitness.plot(skewed_left_axis,sequencing_times,lw = 0.25,color = 'k')
	ax_mean_fitness.plot(skewed_right_axis,sequencing_times,lw = 0.25,color = 'k')


	num = 30
	arrow_x = skewed_right_axis[num] + 0.01
	arrow_dx = skewed_right_axis[-num] - skewed_right_axis[num]
	arrow_y = sequencing_times[num]
	arrow_dy = sequencing_times[-num] - sequencing_times[num]

	# hw = 0.1 * arrow_dx
	# hl = 0.1 * arrow_dy
	# ohg = 0.
	# lw = 0.1
	# line1 = ax_mean_fitness.arrow(arrow_x, arrow_y, arrow_dx,arrow_dy, fc = 'k',  ec='k',  lw = lw, transform = ax_mean_fitness.transData, color = 'k',
	#         head_width=hw, head_length=hl, overhang = ohg, 
 	#         length_includes_head= True, clip_on = True)
	ax_mean_fitness.text(arrow_x, arrow_y, 'Generations', fontsize = fontsize+1)
	for grid_x in numpy.arange(0,round(max_fitness/0.01),1):
		grid = skewed_left_axis + 0.01*(grid_x)*width_scale
		line1, = ax_mean_fitness.plot(grid,sequencing_times,lw = 0.25,color = '0.5',clip_on = False)
		line1.set_dashes((1,1))

	for epoch in numpy.arange(0,max_barcode+1):
		index = epoch*inference_params.INTERVALS_PER_EPOCH

		if epoch == max_barcode:
			index -= 1

		seq_time = times[index]

		grid_left = skewed_left_axis[index]
		grid_right = grid_left + max_fitness *width_scale
		label = '%d.0' % (1+seq_time/110)
		if epoch == max_barcode:
			label = '%d.100' % max_barcode
		ax_mean_fitness.text(grid_right+0.0015,seq_time+30,label, va = 'top', ha = 'left', transform = ax_mean_fitness.transData,fontsize = fontsize)

		if seq_time < 110*max_barcode:
			lw = 0.25
			line1, = ax_mean_fitness.plot([grid_left,grid_right],[seq_time,seq_time],lw = lw,color = '0.5')
			line1.set_dashes((1,1))
		else:
			lw = 0.25
			line1, = ax_mean_fitness.plot([grid_left,grid_right],[seq_time,seq_time],lw = lw,color = 'k')


		
	all_axes = [ax_mean_fitness]
	for axis in all_axes:
		axis.patch.set_alpha(0.)
		for axis_position in ['top','bottom','left','right']:
			axis.axes.spines[axis_position].set_visible(False)
		axis.axes.xaxis.set_ticks_position('none')
		axis.axes.yaxis.set_ticks_position('none')
		axis.set_yticklabels([])
		axis.set_xticklabels([])

	no_divisors = 3*(fitness_range)/0.01 + 1
	if 'C1' in population:		
		bar_chart = [BarChart(fig,ax_fitness[i],numpy.linspace(min_fitness,max_fitness,no_divisors)) for i in range(0,len(ax_fitness))]	
	else:
		bar_chart = [BarChart(fig,ax_fitness[i],numpy.linspace(min_fitness,max_fitness,no_divisors)) for i in range(0,len(ax_fitness))]
	for i in range(0,len(ax_fitness)):
		xticks = ax_fitness[i].get_xticks()
		ax_fitness[i].set_xticks(xticks[1:])
		ax_fitness[i].set_xticklabels(['%d'% (percent*100)  for percent in xticks[1:]])
		ax_fitness[i].tick_params(axis="x",pad=-0.5,direction = "out",which = "both",length = 2)
		# ax_fitness[i].set_ylim([0,1])
		if i > 0:
			ax_fitness[i].set_yticklabels([])

	with_stars = True
	no_error_bar = True

	for epoch in range(0, len(ax_fitness)):
		ax_fitness[epoch].set_ylim([0,1.])
		ax_fitness[epoch].axes.spines['left'].set_visible(False)
		# Only show ticks on the left and bottom spines
		ax_fitness[epoch].axes.yaxis.set_ticks_position('none')
		if epoch > 0:
			ax_fitness[epoch].axes.xaxis.set_ticks_position('none')
			ax_fitness[epoch].axes.spines['bottom'].set_visible(False)
			ax_fitness[epoch].set_xticklabels([])

	for epoch in range(0,len(ax_fitness)):
		if epoch>0:
			ax_mean_fitness.plot(skewed_fitness[(epoch-1)*11:epoch*11+1],sequencing_times[(epoch-1)*11:epoch*11+1],lw = 1,color = 'k')

		if epoch == max_barcode:
			plot_after_start_of_epoch = 10
		else:
			use_epoch = epoch
			plot_after_start_of_epoch = 0

		for ID in clone_list:	
			if use_epoch >= len(ID.split("_"))-1:
				if True:
					if ID in population_tree.keys():
						bar_chart[epoch].add_lineage(clone_dict[ID],use_epoch,no_error_bar = no_error_bar, outline = True, t = plot_after_start_of_epoch)
					else:
						parents = ancestor_list(ID, population_tree)
						for p_ID in parents:
							bar_chart[epoch].add_lineage(clone_dict[ID],
								use_epoch,color = clone_dict[p_ID].color,no_error_bar = no_error_bar, outline = False,update = False, t = plot_after_start_of_epoch)

						bar_chart[epoch].add_lineage(clone_dict[ID],use_epoch,no_error_bar = no_error_bar, outline = True,t = plot_after_start_of_epoch)	
			if with_stars and use_epoch == len(ID.split("_"))-1:
				# place star at establishment time of mutation
				barcode_epoch = len(ID.split("_"))

				# bar_chart[use_epoch].axes.scatter(clone_dict[ID].fitness[use_epoch],0,s = 4,marker = '*',color = clone_dict[ID].color,lw = 0.1,clip_on = False)
				
		if epoch == max_barcode:
			use_epoch -= 1
			plot_after_start_of_epoch = 10
		else:
			use_epoch = epoch
		bar_chart[epoch].add_lineage(clone_dict[''],use_epoch,no_error_bar = no_error_bar, outline = True, t = plot_after_start_of_epoch)	
		
		# fig.savefig(population+'_animated_%d.pdf' % epoch, bbox_inches = 'tight')
		
			# ax_fitness[use_epoch].axes.set_clip_on(False)

	# fig.savefig(OUTFILE,dpi = dpi, bbox_inches = 'tight')		


	# plot panels C and D
	text_panels_CD = {'C1':'c', 'D1':'d'}


	height = 0.15 # specifies height of longest bar in terms of total figure height

	YPD_bottom = 0.57
	bottom = YPD_bottom - popindex* (YPD_bottom - 0.05)
	

	left = 0.69# + 0.4*popindex
	# bottom = 0.05 + 0.4 * popindex
	YPD_bottom = 0.58
	bottom = YPD_bottom - popindex* (YPD_bottom - 0.05)
	height = 0.36*max_barcode/9
	width = 0.83*3/7.5*0.75*fitness_range/0.1
	if choice == 'average':
		width *= 1.2
	ax_mean_fitness = fig.add_axes([left,bottom,width,height])
	ax_mean_fitness.patch.set_alpha(0.)
	ax_mean_fitness.plot(single_inferred_mean_fitness,sequencing_times,lw = 1., color = 'k')

	ax_mean_fitness.text(-0.15*0.1/fitness_range,1.03,text_panels_CD[population],transform = ax_mean_fitness.transAxes,weight = 'bold',fontsize = 8)
 	
	all_axes = [ax_mean_fitness]
	for axis in all_axes:
		axis.axes.xaxis.set_ticks_position('bottom')
		axis.axes.yaxis.set_ticks_position('left')
		for axis_position in []:
			axis.axes.spines[axis_position].set_visible(False)
	

	no_divisors = 3*(fitness_range)/0.01+1

	bar_centers= numpy.linspace(min_fitness,max_fitness,no_divisors)
	ax_mean_fitness.set_ylim([0,110*max_barcode-10])
	
	priority_barcodes = config.highlighted_barcodes[population]

	
	priority_list = []
	other_list = []
	for ID in clone_list:
		if any(ID.startswith(target_ID) for target_ID in priority_barcodes):
			priority_list.append(ID)
		else:
			other_list.append(ID)
	priority_list.extend(other_list)

	grey_intensity = 0.65
	grey_color = (grey_intensity,grey_intensity,grey_intensity,.8)
	clone_dict[''].color = grey_color

	clone_index = 0

	# num_lineages_in_bin
	for ID in priority_list:
		parents = list(ancestor_list(ID,population_tree))
		if len(parents) > 0:
			parent_ID = parents[0]
			last_parent = parents[-1]
		else:
			parent_ID = ''
			last_parent = ''

		if ID in priority_barcodes or parent_ID in priority_barcodes:
			zorder = 2
		else: 
			lineage_dict[ID].color = grey_color #(.2,.2,.2,1.0)#third_color_cycler.get_new_color(alpha = fitness_alpha)
			clone_dict[ID].color = lineage_dict[ID].color
			zorder = 0
		
		#drop shadow
		xval = clone_dict[ID].relative_fitness[0]
		shadow = xval * numpy.ones(len(times))

		shadow_start = 0
		shadow_times = times

		points = np.array([shadow, shadow_times]).T.reshape(-1, 1, 2)
		segments = np.concatenate([points[:-1], points[1:]], axis=1)

		segment_colors = []
		R, G, B, A = clone_dict[ID].color

		min_log10_frequency_to_show = 4.5
		for t in range(0, len(shadow)-1):
			if clone_dict[ID].freqs[shadow_start+t] > 10.**-min_log10_frequency_to_show:
				# alpha = clone_dict[ID].freqs[shadow_start + t]
				alpha = ((min_log10_frequency_to_show+np.log(clone_dict[ID].freqs[t])/np.log(10))**1)/(min_log10_frequency_to_show)**1.

			else:
				alpha = 0.
			# if (t) % 11 > 2 and (t) % 11 <= 9 :
			# 	alpha = 0.
			# else:
			# 	pass
				# print shadow_times[t],
			segment_color = (R,G,B,alpha)
			segment_colors.append(segment_color)
			# ax_mean_fitness.plot(shadow[t:t+2],shadow_times[t:t+2]-5,color = segment_color,lw = 1.)
		lc = LineCollection(segments,colors = segment_colors,zorder = zorder,lw = 1)
		ax_mean_fitness.add_collection(lc)


		if ID in priority_barcodes or parent_ID in priority_barcodes:
			parent_fitness = clone_dict[last_parent].relative_fitness[0]

			dx = xval - parent_fitness
			t_start =  11*(len(ID.split("_"))-2)*10 + 30
			t_end = t_start + 70

			dx = clone_dict[ID].relative_fitness[0] - parent_fitness
			dt = 0
			t_appear = 11*(len(ID.split("_"))-1)*10 + 10*shadow_start
			if dx < 0.002:
				mutation_scale = 2
			else:
				mutation_scale = 4

			lw = 0.5 

			ax_mean_fitness.annotate('',
			            xy=(parent_fitness+dx, t_appear), xycoords='data',
			            xytext=(parent_fitness,t_appear), textcoords='data',
			            arrowprops=dict(arrowstyle="->",
			                            connectionstyle="angle3,angleA=50,angleB=130",lw =lw,color= (R,G,B,.8),shrinkA = 0,shrinkB = 0,mutation_scale = mutation_scale,zorder = 40))
	R, G, B, A = clone_dict[''].color

	xval = clone_dict[''].relative_fitness[0]
	shadow = xval * numpy.ones(len(sequencing_times))
	points = np.array([shadow, sequencing_times]).T.reshape(-1, 1, 2)
	segments = np.concatenate([points[:-1], points[1:]], axis=1)

	segment_colors = []
	for t in range(0, len(shadow)-1):
		if clone_dict[''].freqs[shadow_start+t] > 10.**-min_log10_frequency_to_show:
			# alpha = clone_dict[ID].freqs[shadow_start + t]
			alpha = ((min_log10_frequency_to_show+np.log(clone_dict[''].freqs[t])/np.log(10))**1)/(min_log10_frequency_to_show)**1.
		else:
			alpha = 0.

			# print shadow_times[t],
		segment_color = (R,G,B,alpha)
		segment_colors.append(segment_color)
		# ax_mean_fitness.plot(shadow[t:t+2],shadow_times[t:t+2]-5,color = segment_color,lw = 1.)
	

	lc = LineCollection(segments,colors = segment_colors,zorder = zorder,lw = 1)
	ax_mean_fitness.add_collection(lc)	


	ax_mean_fitness.set_xlim([min_fitness,max_fitness])
	ticks =  numpy.arange(min_fitness,max_fitness+0.01,0.01)
	ax_mean_fitness.set_xticks(ticks)
	ax_mean_fitness.set_xticklabels(['%d' %(round(100*tick)) for tick in ticks])
	ax_mean_fitness.tick_params(axis="x",pad=1,direction = "in",which = "both",length = 2)

	ax_mean_fitness.set_yticks([(110*t) for t in range(0,max_barcode+1)])
	ax_mean_fitness.set_yticklabels(["%d.0"%(i+1) for i in range (0,max_barcode)] + ["%d.100"%(max_barcode)])


	ax_mean_fitness.set_xlabel('Fitness relative\nto ancestor (%s)'% unit,labelpad = 0)
	ax_mean_fitness.set_ylabel('Generations',labelpad = 0)

sys.stderr.write("Saving figure to %s \n"% OUTFILE)
fig.savefig(OUTFILE,format = 'pdf')			


