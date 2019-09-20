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


matplotlib.rcParams['axes.linewidth'] = .25

parser = argparse.ArgumentParser()
parser.add_argument("-dpi", default = 300, type = int, help = 'resolution (dots per inch)')
parser.add_argument("-width", default = 4, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 3.6, type = float, help = 'figure height (inches)')
parser.add_argument("-outfile", default = 'Figure4_final.pdf', type = str, help = 'figure name')
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
OUTFILE = args.outfile

matplotlib.rcParams['xtick.labelsize'] = 6
matplotlib.rcParams['ytick.labelsize'] = 6
matplotlib.rcParams['font.size'] = 7

fontsize = 7

choice = args.choice


if choice in ['evolution', 'average']:
	fitness_unit = 'percent'
else:
	fitness_unit = 'per cycle'

if choice == 'average':
	subdirectory = 'main/'
else:
	subdirectory = 'other/'

pop_labels  = config.pop_labels
pop_colors = config.pop_colors

OUTFILE = config.figure_directory + subdirectory+ OUTFILE.split(".")[0] + "_" + choice + ".pdf"

print OUTFILE
# set up figures and axes
fig = plt.figure(figsize = (fig_width,fig_height))
fig_ed_success,ax_ed_success = plt.subplots(figsize = (2,2))

left, bottom, width, height = 0.0,0.56, 0.3, 0.42
ax_success = fig.add_axes([left,bottom,width,height])

left += 0.35
width = 0.24
ax_no_mutations = fig.add_axes([left, bottom, width, height])

left += 0.38
height = 0.28
ax_hist_YPA = fig.add_axes([left,bottom-0.07+0.18,width,0.256])
ax_hist = fig.add_axes([left,bottom-0.07,width,0.16])
ax_parental_fitness = fig.add_axes([left,0.,width,0.4])


left, bottom, width, height = 0.04, 0., 0.52, 0.4
ax_parental_frequency = fig.add_axes([left,bottom,width,height])


fig_ed, ax_ed = plt.subplots(2,2,figsize = (5,2.5))
ax_number = ax_ed[0][0]
ax_freq_entropy = ax_ed[1][0]
ax_average_fitness_variance = ax_ed[0][1]
ax_fitness_entropy = ax_ed[1][1]

fig_variance, ax_variance = plt.subplots(2,1,figsize=(2,3))
ax_evolution_variance = ax_variance[0]
ax_barcoding_variance = ax_variance[1]


all_axes = [ax_parental_fitness, ax_parental_frequency, ax_hist, ax_hist_YPA, ax_success, ax_no_mutations]
all_axes.extend([ax_number, ax_freq_entropy])
all_axes.extend([ax_evolution_variance, ax_barcoding_variance])
all_axes.extend([ax_average_fitness_variance,ax_fitness_entropy])
all_axes.extend([ax_ed_success])


# ax_success.text(-.36,1.12,'Figure 4',transform = ax_success.transAxes,weight = 'bold',fontsize = 8)
for axes in all_axes:
	axes.patch.set_alpha(0.)
	for axis_position in ['top','right']:
		axes.axes.spines[axis_position].set_visible(False)
	axes.axes.xaxis.set_ticks_position('bottom')
	axes.axes.yaxis.set_ticks_position('left')




popindex = -1
for population in reversed(config.populations):
	popindex += 1
	print population

	max_barcode = config.max_barcode[population]

	clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)
	times = times[:max_barcode*inference_params.INTERVALS_PER_EPOCH]
	mutant_tree = population_tree
	clone_list.append("")

	mean_fitness = {'evolution': numpy.zeros(len(times)),
					'barcoding':numpy.zeros(len(times)),
					'average':numpy.zeros(len(times))}
	fitness_variance = {'evolution':numpy.zeros(len(times)), 
						'barcoding':numpy.zeros(len(times)),
						'average':numpy.zeros(len(times))}

	for bcd, clone in clone_dict.iteritems():
		clone.freqs = clone.freqs[:max_barcode*inference_params.INTERVALS_PER_EPOCH]
		mean_fitness['evolution'] += clone.freqs*clone.evolution_fitness
		mean_fitness['barcoding'] += clone.freqs*clone.barcoding_fitness
		mean_fitness['average'] += clone.freqs*(clone.evolution_fitness*100+clone.barcoding_fitness)/200

		fitness_variance['evolution'] += clone.freqs*(clone.evolution_fitness)**2
		fitness_variance['barcoding'] += clone.freqs*(clone.barcoding_fitness)**2
		fitness_variance['average'] += clone.freqs*((clone.evolution_fitness*100+clone.barcoding_fitness)/200)**2

	for condition in mean_fitness.keys():
		fitness_variance[condition] -= (mean_fitness[condition])**2


	if choice == 'evolution':	
		inferred_mean_fitness = mean_fitness[choice]
		for bcd,clone in clone_dict.iteritems():
			clone.relative_fitness = numpy.ones(len(clone.relative_fitness))*clone.evolution_fitness
	elif choice == 'barcoding':
		inferred_mean_fitness = mean_fitness[choice]/100.
		for bcd,clone in clone_dict.iteritems():
			clone.relative_fitness = numpy.ones(len(clone.relative_fitness))*clone.barcoding_fitness/100.
	elif choice == 'average':
		inferred_mean_fitness = 0.5*(mean_fitness['evolution'] + mean_fitness['barcoding']/100.)
		for bcd,clone in clone_dict.iteritems():
			clone.relative_fitness = numpy.ones(len(clone.relative_fitness))*0.5*(clone.evolution_fitness + clone.barcoding_fitness/100.)

	entropy = numpy.zeros(len(inferred_mean_fitness))
	number_of_segregating_clones = numpy.zeros(len(inferred_mean_fitness))
	number_of_high_frequency_clones = numpy.zeros(len(inferred_mean_fitness))

	max_parents = numpy.zeros(len(inferred_mean_fitness))
	fitness_entropy = numpy.zeros(len(inferred_mean_fitness))

	for ID in clone_dict.keys():
		clone = clone_dict[ID]

		clone_fitness = clone.relative_fitness[0]

		positive_mask = clone.freqs > 10**-4
		number_of_segregating_clones += positive_mask
		
		entropy[positive_mask] += - clone.freqs[positive_mask] * numpy.log(clone.freqs[positive_mask]) 
		if ID == "":
			num_parents = 0
		else:
			num_parents = sum(1 for anc_ID in ancestor_list(ID,population_tree))

		max_parents = numpy.maximum(max_parents,(num_parents+1)*positive_mask) 

		high_frequency_mask = clone.freqs > 10**-2
		number_of_high_frequency_clones += high_frequency_mask

		contribution_to_variance = clone.freqs * (clone_fitness - mean_fitness[choice])**2
		contribution_to_variance /= fitness_variance[choice]

		fitness_entropy[positive_mask] -= contribution_to_variance[positive_mask]*numpy.log(contribution_to_variance[positive_mask])

	
	def find_last_parent(ID,population_tree):
		if ID == "":
			return ''
		else:
			parents = ancestor_list(ID,population_tree)
			i = 0
			for p_ID in parents:
				i+= 1
			if i > 0:
				return p_ID
			else:
				return ''

	fitness_effects = []
	effect_dict = {}

	for ID in clone_dict.keys():
		parent_ID = find_last_parent(ID, population_tree)
		diff = clone_dict[ID].relative_fitness[0] - clone_dict[parent_ID].relative_fitness[0]

		effect_dict.update({ID:diff})


	bins = numpy.linspace(0,0.04,13)

	est_times = [(i-1)*110+55 for i in range(1,11)]

	counts = []

	for i in range(1,11):
		counts.append(sum(1 for ID in effect_dict.keys() if len(ID.split("_")) == i))
		ax_number.axes.bar(est_times[i-1],counts[-1],width = 40,color = pop_colors[population],alpha = 1.0,align = 'center',lw = 0.)

		ax_number.plot(times,number_of_segregating_clones,lw = 1, color = pop_colors[population])
	ax_number.set_ylabel('No. lineages',labelpad = 1)

	if population == config.populations[-1]:
		ax_number.plot([],[],lw = 1,color = 'k',label = r'$f>$'+'0.01' + r'$\%$')
		ax_number.plot([],[],lw = 5,color = '0.4', label = 'new')
	ax_number.legend(handlelength = 1., handletextpad = 0.5,labelspacing = 0.,loc = 2,fontsize = 7)


	ax_freq_entropy.plot(times,entropy,lw = 1, color = pop_colors[population],clip_on = False)
	ax_freq_entropy.set_ylabel('Entropy',labelpad = 4)
	ax_freq_entropy.set_xlabel('Generations',labelpad = 1)

	ax_average_fitness_variance.plot(times,fitness_variance[choice],lw = 1, color = pop_colors[population],clip_on = False)
	ax_average_fitness_variance.set_ylabel('Variance in fitness',labelpad = 1)

	ax_fitness_entropy.plot(times,fitness_entropy,lw = 1, color = pop_colors[population],clip_on = False)
	ax_fitness_entropy.set_ylabel('Fitness entropy',labelpad = 1)
	ax_fitness_entropy.set_xlabel('Generations',labelpad = 1)



	def count_children(ID, population_tree):
		return sum(1 for ch_ID in child_list(ID,population_tree))

	fitness_effects = []
	parental_relative_fitnesses = []
	parental_frequencies = []

	for ID in clone_list:
		parent_ID = find_last_parent(ID, population_tree)
		diff = clone_dict[ID].relative_fitness[0] - clone_dict[parent_ID].relative_fitness[0]

		effect_dict.update({ID:diff})

		no_epochs = len(ID.split("_"))-1

		parental_relative_fitness = clone_dict[parent_ID].relative_fitness[0] - sum(inferred_mean_fitness[no_epochs*11:no_epochs*11+10])/10.
		parental_average_frequency = sum(clone_dict[parent_ID].freqs[no_epochs*11:no_epochs*11+10])/10.
		
		ax_parental_fitness.scatter(diff,parental_relative_fitness,marker = (300,0,0),s = 20,lw = 0, color = pop_colors[population],alpha = 0.8,clip_on = False)

		if no_epochs > 2:
			parental_relative_fitnesses.append(parental_relative_fitness)
			parental_frequencies.append(parental_average_frequency)

	effect_ticks = numpy.arange(-0.01,0.05,0.01)
	ax_parental_fitness.set_xlim([-0.01,0.04])

	effect_labels = ['%.1f' % (100.*percent) for percent in effect_ticks]
	ax_parental_fitness.set_xticks(effect_ticks)
	ax_parental_fitness.set_xticklabels(effect_labels)
	ax_parental_fitness.set_xlabel('Fitness effect (%s)' % fitness_unit)
	if choice == 'evolution':
		ax_parental_fitness.set_ylim([-0.03,0.03])
		parental_fitness_ticks = [-0.03,0,0.03]
	else:
		ax_parental_fitness.set_ylim([-0.02,0.02])
		parental_fitness_ticks = [-0.02,0,0.02]

	effect_labels = ['%.1f' % (100.*percent) for percent in parental_fitness_ticks]
	ax_parental_fitness.set_yticks(parental_fitness_ticks[::2])
	ax_parental_fitness.set_yticklabels(effect_labels[::2])
	# ax_parental_fitness.set_ylabel('Parent relative fitness (%s)' % fitness_unit,labelpad = 1)

	ax_parental_fitness.plot([0.03,-0.01],[-0.03,0.01],color = '0.5',lw = 1,zorder = 0)

	bins = numpy.linspace(-0.01,0.04,16)
	if population == 'C1':
		ax_hist.hist(effect_dict.values(),bins = bins, alpha = 0.8, lw = 0, color = pop_colors[population],zorder = popindex)
	else:
		ax_hist_YPA.hist(effect_dict.values(),bins = bins, alpha = 0.8, lw = 0, color = pop_colors[population],zorder = popindex)
	
	for histogram in [ax_hist,ax_hist_YPA]:
		histogram.set_xlim(-0.01,0.04,13)

		effect_ticks = ax_parental_fitness.get_xticks()
		effect_labels = ['%.1f' % (100.*percent) for percent in effect_ticks]
		histogram.set_xticks(effect_ticks)
		histogram.set_ylabel('Number',labelpad = 1)

	ax_hist.set_xticklabels(effect_labels)
	ax_hist_YPA.set_xticklabels([])
	ax_hist.set_ylim([0,25])
	# ax_hist.set_xlim(ax_parental_fitness.get_ylim())

	# ax_hist.set_yticks([0,10,20,30])
	
	# ax_hist.set_xlabel('Fitness effect (%s)' % fitness_unit,labelpad = 1)


	fitness_dict = {}

	for ID in mutant_tree.keys():
		ancestor_fitness = clone_dict[''].relative_fitness
		diff = clone_dict[ID].relative_fitness[0] - ancestor_fitness[0]
		fitness_dict.update({ID:diff})

	all_fitnesses = []

	lineage_info = []
	first_timepoint = len(times[0:max_barcode*11-1])/2
	num_lineages= sum(1 for thing in fitness_dict.keys() )
	for key in sorted(mutant_tree.keys()):
		all_fitnesses.append(fitness_dict[key])
		no_early_children = len( [x for x in child_list(key,mutant_tree) if len(x.split("_")) <= 5 ] )

		fitness_rank = sum(1 for thing in fitness_dict.keys() if fitness_dict[key] >= fitness_dict[thing] )#and  fitness_dict[thing][1] < fitness_dict[thing][0])
		fitness_rank *= 1./num_lineages

		freq_rank = sum(1 for thing in fitness_dict.keys() if max(lineage_dict[key].freqs[first_timepoint:11*max_barcode]) >= max(lineage_dict[thing].freqs[first_timepoint:11*max_barcode]))
		freq_rank *= 1./num_lineages


		lineage_info.append((fitness_rank,max(lineage_dict[key].freqs[first_timepoint:11*max_barcode]),no_early_children))

		s = min(10*(0.5+no_early_children)**2,300*4.5**2)
		kwargs = dict(marker = (300,0,0),color= config.pop_colors[population],s =s,alpha = 0.8,clip_on=False,zorder = 10-popindex)
		
		ax_ed_success.scatter(fitness_rank,max(lineage_dict[key].freqs[first_timepoint:11*max_barcode])+10**-5,**kwargs)

	ax_ed_success.plot([],[],lw = 8,color = config.pop_colors[population],label = config.pop_labels[population],alpha = 0.8)

	ax_ed_success.set_xlim([0,1])
	ax_ed_success.set_ylim([10**-6,1.1])
	ax_ed_success.set_xlabel('Fitness rank of founding mutation',labelpad = 1)
	ax_ed_success.set_yscale('log')
	ax_ed_success.set_ylabel('Maximum frequency\n in 2nd half of experiment',labelpad = 1)
	if population == 'C1':
		for no_children in range(0,7):
			ax_ed_success.scatter([],[],marker = (300,0,0),color = '#beae8a',alpha =  0.3,  s = min(10*(0.5+no_children)**2,300*4.5**2),label = no_children)

		ax_ed_success.legend(bbox_to_anchor = (1.4,0.),scatterpoints = 1,labelspacing = 1,handlelength = 1,fontsize = 7)

	quantiles = numpy.asarray([0.25, 0.5,0.75,1.],dtype = float)
	original_quantiles = numpy.asarray([0.25, 0.5,0.75,1.],dtype = float)
	means = numpy.zeros(len(quantiles))
	mins = numpy.ones(len(quantiles))
	maxs = numpy.zeros(len(quantiles))
	nums = numpy.zeros(len(quantiles))

	mut_categories = numpy.asarray([0,1,2],dtype = float)
	means_mut = numpy.zeros(len(mut_categories))
	mins_mut = numpy.ones(len(mut_categories))
	maxs_mut = numpy.zeros(len(mut_categories))
	nums_mut = numpy.zeros(len(mut_categories))

	for lin_info in lineage_info:
		fit, freq, no_children = lin_info

		for q in range(0,len(quantiles)):
			if fit < quantiles[q]:
				break

		if freq > maxs[q]:
			maxs[q] = freq

		if freq < mins[q]:
			mins[q] = freq
		means[q] += freq
		nums[q] += 1.

		for muts in range(0,len(mut_categories)):
			if no_children <= mut_categories[muts]:
				break
		if freq > maxs_mut[muts]:
			maxs_mut[muts] = freq

		if freq < mins_mut[muts]:
			mins_mut[muts] = freq
		means_mut[muts] += freq
		nums_mut[muts] += 1.

	means /= nums
	means_mut /= nums_mut

	quantiles -= quantiles[0]/2.
	if population == 'C1':
		quantiles += -0.02
		mut_categories -= 0.1
	else:
		quantiles += 0.02
		mut_categories += 0.1

	ax_success.scatter(quantiles,mins+10**-5,marker = '_', color = pop_colors[population],alpha = 0.8,clip_on = False,zorder = 2,s= 10)
	ax_success.scatter(quantiles,maxs,marker = '_', color = pop_colors[population],alpha = 0.8,clip_on = False, zorder = 2,s = 10)
	ax_success.scatter(quantiles,means,marker = (300,0,0), color = pop_colors[population],alpha = 0.8,clip_on = False,zorder = 2,s = 6)
	for q in range(0,len(quantiles)):
		xpos = quantiles[q]
		ax_success.plot([xpos,xpos],[mins[q]+10**-5,maxs[q]],color = pop_colors[population], lw = 0.5,zorder = 1,clip_on = False)

	ax_no_mutations.scatter(mut_categories,mins_mut+10**-5,marker = '_', color = pop_colors[population],alpha = 0.8,clip_on = False,zorder = 2,s = 10)
	ax_no_mutations.scatter(mut_categories,maxs_mut,marker = '_', color = pop_colors[population],alpha = 0.8,clip_on = False, zorder = 2,s = 10)
	ax_no_mutations.scatter(mut_categories,means_mut,marker = (300,0,0), color = pop_colors[population],alpha = 0.8,clip_on = False,zorder = 2,s = 6)
	for q in range(0,len(mut_categories)):
			
			ax_no_mutations.plot([mut_categories[q],mut_categories[q]],[mins_mut[q]+10**-5,maxs_mut[q]],color = pop_colors[population], lw = 0.5,zorder = 1,clip_on = False)

	ax_success.set_xlim([0,1])
	ax_success.set_xlabel('Fitness rank \nof founding mutation',labelpad = 1)
	ax_success.set_yscale('log')
	ax_no_mutations.set_yscale('log')
	ax_no_mutations.set_yticks([])
	ax_no_mutations.set_xlabel(r'No. subsequent mutations'+'\nin 1st half of experiment',labelpad = 1)
	ax_success.set_ylabel('Max. frequency\n in 2nd half of experiment',labelpad = 1)


	ax_success.set_ylim([10**-6,1.1])
	ax_no_mutations.set_ylim([10**-6,1.1])
	ax_no_mutations.set_xlim([-1,2.5])
	ax_no_mutations.set_xticks(mut_categories)
	ax_no_mutations.set_xticklabels([0,1,r'$\geq$'+'2'])

	xticks = [0]
	xticks.extend(original_quantiles.tolist())
	ax_success.set_xticks(xticks)

	axes_labels = ['a', 'b', 'c', 'd', 'e', 'a','b','a','b','c','d']
	dx = -0.2
	text_left = [1.8*dx,dx,dx,1.5*dx,1.15*dx,dx,dx,1.4*dx,1.4*dx,1.*dx,1.*dx]
	

	axes_index = 0
	labeled_axes = [ax_success,ax_no_mutations,ax_parental_frequency,ax_hist, ax_parental_fitness]
	labeled_axes.extend([ax_number, ax_freq_entropy])
	labeled_axes.extend([ax_evolution_variance,ax_barcoding_variance])
	labeled_axes.extend([ax_average_fitness_variance,ax_fitness_entropy])

	for axes in labeled_axes:
		label_xpos =text_left[axes_index]
		label_ypos = 1.05
		if axes_labels[axes_index] in['a', 'b','c','d']:
			label_ypos = 1.00
		if axes == ax_hist:
			label_ypos = 2.7
		# print label_xpos
		axes.text(label_xpos,label_ypos, axes_labels[axes_index], weight = 'bold', fontsize = 8, transform = axes.transAxes)
		axes_index += 1

	fitnesses, frequencies, no_descendants = [], [], []

	for ID in clone_list:
		clone_epoch = len(ID.split("_")) - 1
		if ID == "":
			clone_children = population_tree.keys()
		else:
			clone_children = child_list(ID, population_tree)

		fit = [clone_dict[ID].relative_fitness[0] - sum(inferred_mean_fitness[epoch*11:epoch*11+10])/10. for epoch in range(0,max_barcode)]
		freq = [ sum(clone_dict[ID].freqs[epoch*11:epoch*11+10])/10. for epoch in range(0,max_barcode)]
		num = [sum([1 for child_ID in clone_children if len(child_ID.split("_"))==epoch + 1 ]) for epoch in range(0,max_barcode)]
		
		fitnesses.append(fit)
		frequencies.append(freq)
		no_descendants.append(num)


	# fitnesses, frequencies, spot_colors, p_values, no_descendants = avg_vs_establishment_per_epoch.calculate_all(population,num_trials = 10**2)
	angle = 45./180*numpy.pi

	for i in range(0,len(clone_list)):
		ID = clone_list[i]

		clone_arising_epoch = len(ID.split("_")) - 1

		contour_data = []
		for epoch in range(clone_arising_epoch,max_barcode):

			shape_kwargs = dict(alpha = 0.5, marker = (300,0,0),zorder = 0,s = 2,lw = 0,clip_on = True)
			max_s = 90
			if no_descendants[i][epoch]>0:
				s = min(5*no_descendants[i][epoch]**2+10,max_s)
				shape_kwargs.update(dict(marker = 's',zorder = 1,s= s, lw = 0., alpha = .8,edgecolor = '#beae8a',clip_on = False))


			# shape_kwargs.update(dict(lw = 0.4, edgecolor = 'r',zorder = 2))

			ax_parental_frequency.scatter(frequencies[i][epoch],fitnesses[i][epoch],color = pop_colors[population], **shape_kwargs)
			if no_descendants[i][epoch]>3.:
				if frequencies[i][epoch]<0.8:
					ax_parental_frequency.text(frequencies[i][epoch]*1.22,fitnesses[i][epoch]+0.0015,'%d'%no_descendants[i][epoch],transform = ax_parental_frequency.transData,color = 'k',fontsize = 5)
				else:
					xy = (frequencies[i][epoch],fitnesses[i][epoch])
					angle -=35./180*numpy.pi + 48./180*numpy.pi * (1-popindex)

					xy_text = (xy[0]*2**(1.2*numpy.cos(angle)),xy[1]+0.008*numpy.sin(angle))
					ax_parental_frequency.annotate('%d'%no_descendants[i][epoch],xy = xy, xytext = xy_text,fontsize = 5,arrowprops = dict(width = 0.002,zorder = 20,shrink = 0,headwidth = 0.,linewidth = 0.7,color = pop_colors[population]))

	ax_parental_frequency.set_xlim([10**-4,1.])
	if choice == 'evolution':
		ax_parental_frequency.set_ylim([-0.03,0.03])
	else:
		ax_parental_frequency.set_ylim([-0.02,0.02])
	ax_parental_frequency.set_yticks(ax_parental_fitness.get_yticks())
	parental_fitness_ticks = ax_parental_fitness.get_yticks()
	effect_labels = ['%.1f' % (100.*percent) for percent in parental_fitness_ticks]
	ax_parental_frequency.set_yticklabels(effect_labels)
	ax_parental_frequency.set_xscale('log')

	ax_parental_frequency.set_xlabel('Mean frequency of background',labelpad = 1)
	ax_parental_frequency.set_ylabel('Parent fitness\nrelative to mean (%s)' % fitness_unit,labelpad = 1)

	### plot_variance in fitness
	dt={'evolution':10.,'barcoding':1.}

	for environment, axis in zip(['evolution','barcoding'],[ax_evolution_variance,ax_barcoding_variance]):
		mean_fitness_increase = numpy.diff(mean_fitness[environment])
		if environment =='evolution':
			mean_fitness_increase = gaussian_filter(mean_fitness_increase,1)

		axis.plot(times,fitness_variance[environment],color = pop_colors[population],lw = 1.5,alpha = 0.8)
		ylim = axis.get_ylim()

		axis.plot(times[1:],mean_fitness_increase/dt[environment],color = pop_colors[population],lw = 1,alpha = 0.5)
		
		axis.set_ylim(ylim)


for axes in [ax_number,ax_freq_entropy,ax_evolution_variance,
			 ax_barcoding_variance,ax_average_fitness_variance,ax_fitness_entropy]:
	axes.set_xlim([0,1100-10])
	axes.set_xticks([110*x for x in range(0,10)])
	axes.set_xticklabels(['%d.0' % (x+1) for x in range(0,10)])


ax_number.set_xticklabels([])
ax_average_fitness_variance.set_xticklabels([])


ax_parental_frequency.scatter([],[], marker=(100, 0, 0),s = 5,lw = 0,alpha = 1.,clip_on = True,edgecolor = 'w',color = '#beae8a',label = r'0')
for no_descendants in range(1,5):
	s = min(5*no_descendants**2+10,max_s)
	if s < max_s:
		label = '%d' % no_descendants
	else:
		label = r'$\geq$'+'%d'%(no_descendants)
	ax_parental_frequency.scatter([],[], marker = 's',s = s,lw = 0,alpha = 1.,clip_on = True,edgecolor = '#beae8a',color = '#beae8a',label = label)
	if s == max_s:
		break
leg = ax_parental_frequency.legend(scatterpoints = 1,loc = 4,bbox_to_anchor = (1.23,0.6),handletextpad = 0.5,frameon = True,borderpad = 0.5,labelspacing = 0.55,ncol = 1,fontsize = 6)
leg.get_frame().set_linewidth(0.)

for population in config.populations:
	ax_hist.plot([],[],lw = 8,color = pop_colors[population],label = pop_labels[population],alpha = 0.8)
leg = ax_hist.legend(bbox_to_anchor = (1.2,2.4),handlelength = 1.2,title= 'Population',frameon=True,borderpad = .6,ncol = 1,fontsize = 7)
leg.get_frame().set_linewidth(0.25)




ax_number.set_yticks([30*x for x in range(0,4)])

ax_success.legend(bbox_to_anchor = (1.28,0.2),scatterpoints = 1,labelspacing = .8,handlelength = 1,fontsize = 6)

# set up variance in fitness axes etc.
for population in config.populations:
	ax_evolution_variance.plot([],[],lw = 5, color =config.pop_colors[population], label = config.pop_labels[population])
ax_evolution_variance.plot([],[],lw = 1.5, color ='k',alpha = 0.8, label = 'direct')
ax_evolution_variance.plot([],[],lw = 1, color = 'k',alpha = 0.5,label = 'rate of adaptation')

ax_evolution_variance.legend(handlelength = 1.2, handletextpad = 0.5,labelspacing = 0.,loc =2,fontsize = 7)

ax_evolution_variance.set_yticks([0,1*10**-4,2*10**-4,3*10**-4])
ax_barcoding_variance.set_yticks([0,1,2,3])
# ax_evolution_variance.set_ylim([10**-6,10**-4])
# ax_evolution_variance.set_yscale('log')

for axis in [ax_evolution_variance, ax_barcoding_variance,ax_average_fitness_variance]:
	axis.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
	mf.set_powerlimits((0,0))
	axis.yaxis.set_major_formatter(mf)
	# axis.set_xticklabels([])
ax_average_fitness_variance.set_ylim([0,0.00025])
ax_barcoding_variance.set_xlabel('Generation',labelpad=0)



# ax_number.text(0.,1.1,'Extended Data Figure 4',fontsize = 8,weight = 'bold',transform = ax_number.transAxes)
# ax_ed_success.text(0.,1.1,'Extended Data Figure 3',fontsize = 8,weight = 'bold',transform = ax_ed_success.transAxes)
ax_evolution_variance.set_ylabel('Evolution fitness\nvariance (percent'+ r'$^{\mathrm{2}}$)')
ax_barcoding_variance.set_ylabel('Barcoding fitness\nvariance (per procedure'+ r'$^{\mathrm{2}}$)')
fig_variance.savefig(config.figure_directory+'si/fitness_variance.pdf',bbox_inches = 'tight')
fig.savefig(OUTFILE,bbox_inches = 'tight')
fig_ed.savefig(config.figure_directory+'ed/EDFigure4_final.pdf',bbox_inches = 'tight')
fig_ed_success.savefig(config.figure_directory+'ed/EDFigure3_final.pdf',bbox_inches = 'tight')



	

# 		