import sys, os, glob
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col

# add custom modules to path
sys.path.insert(0,'../modules/') 
sys.path.insert(0,'../wgs_scripts/')

# import custom modules
import config 
import local_matplotlibrc
import produce_annotation_map

input_directory = config.wgs_data_directory+'filtered/'

color_dict = {}

fig,axes = plt.subplots(1,len(config.populations),figsize = (2.625*len(config.populations),2.2))
fontsize = 7


color_set = plt.get_cmap('Set1') 


def get_new_color(color_index):
	return scalarMap.to_rgba(color_index)

popindex = -1
for population in config.populations:
	popindex += 1

	max_barcode = config.max_barcode[population]

	sys.stderr.write('Plotting timecourses from population %s...\n' % (population))
	filenames =  glob.glob(input_directory+population+'*filtered_timecourse.txt')

	min_cross_time = 0
	trajectories = []

	for filename in filenames:
		f = open(filename,'r')
		print filename
		f.readline()

		for line in f:
			items = line.strip().split(', ')
			label = ";".join(items[0:3])
			times = numpy.asarray([(float(item.split('.')[0])-1)*110 + float(item.split('.')[1]) for item in items[3].split(' ')],dtype = float) 

			counts = numpy.asarray([float(item.strip()) for item in items[4].split(' ')],dtype = float) 
			depths = numpy.asarray([int(item.strip()) for item in items[5].split(' ')],dtype = long)

			pvalues = numpy.asarray([float(item.strip()) for item in items[6].split(' ')],dtype = float)

			chromosome = items[0]
			position = int(items[1])

			freqs = counts/depths
			min_cross_time = next(f for f in range(0,len(freqs)) if sum(freqs[:f])>0.3)
			min_cross_time += sum(freqs*times)/sum(freqs)/1.1
			label, mut_type = produce_annotation_map.annotate_mutation(', '.join(items[0:3]))
			print label, mut_type
			
			trajectories.append([min_cross_time,freqs,label,mut_type])



		f.close()
	trajectories = sorted(trajectories, key=lambda x: x[0])

	max_cross_time = trajectories[-1][0]
	color_wheel_size = max_cross_time
	cNorm  = col.Normalize(vmin=0, vmax=color_wheel_size*(1))
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=color_set)

	i = 0

	for traj in trajectories:
		i += 1

		min_cross_time, freqs, label, mut_type = traj
		color = get_new_color((min_cross_time-400)%color_wheel_size)
		line1, = axes[popindex].plot(times,freqs,label = label,lw = 1,color= color)
		if mut_type == 'syn' or mut_type == 'intergenic':
			line1.set_dashes((3,2))

	axes[popindex].set_ylim([0,1])
	axes[popindex].set_xticks([0]+[110*x-10 for x in range(1,max_barcode+1)])
	axes[popindex].set_xticklabels(["1.0"]+["%d.%d" % (epoch+1,100) for epoch in range(0,max_barcode)])
	axes[popindex].set_xlim([0,110*max_barcode])
	for tick in axes[popindex].get_xticklabels():
   		tick.set_rotation(90)
   		tick.set_horizontalalignment('center')
	# axes[popindex].set_xlim([0,1100])

	
	axes[popindex].set_xlabel('Generations',fontsize = fontsize)
	axes[popindex].spines['right'].set_visible(False)
	axes[popindex].spines['top'].set_visible(False)
	# Only show ticks on the left and bottom spines
	axes[popindex].yaxis.set_ticks_position('left')
	axes[popindex].xaxis.set_ticks_position('bottom')
	
	if population == 'C1':
		# axes[popindex].text(-0.2,1.15,'Extended Data Figure 1',fontsize = 7,weight = 'bold',transform = axes[popindex].transAxes)
		panel_label = 'a  '
	else:
		panel_label = 'b  '
	axes[popindex].text(-100,1.07, panel_label + config.pop_labels[population],fontsize = 8,weight = 'bold')
	
	line1,= axes[popindex].plot([0,1100],[0.5,0.5],color = 'k',lw = 0.5)
	line1.set_dashes((1,1))
	if population == 'C1':
		legend_kw = dict(loc = 'best')
	else:
		legend_kw = dict(bbox_to_anchor= (2.2,0.25),ncol = 2)
	axes[popindex].legend(fontsize = fontsize - 1,**legend_kw)
axes[0].set_ylabel('Allele frequency',fontsize = fontsize)
plt.savefig(config.figure_directory+'/ed/EDFigure1_final.pdf',bbox_inches = 'tight')

sys.stderr.write('Done!\n')