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
import subprocess

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
# import local_matplotlibrc
import lineage.file_parser as file_parser
from lineage.read_clone_data import *
from lineage.plot_utils_clone import * # contains muller and bar chart class
from lineage.tree_utils import *


parser = argparse.ArgumentParser()
parser.add_argument("-dpi", default = 600, type = int, help = 'resolution (dots per inch)')
parser.add_argument("-width", default = 6, type = float, help = 'figure width (inches)')
parser.add_argument("-height", default = 2.6, type = float, help = 'figure height (inches)')
parser.add_argument("-fontsize", default = 8, type = int, help = 'fontsize (pt)')
parser.add_argument("-outfile", default = 'Figure2.png', type = str, help = 'figure name')


args = parser.parse_args()

dpi = args.dpi
fig_width = args.width
fig_height = args.height
fontsize = args.fontsize
OUTFILE = args.outfile

matplotlib.rcParams['xtick.labelsize'] = fontsize
matplotlib.rcParams['ytick.labelsize'] = fontsize
matplotlib.rcParams['font.size'] = fontsize + 1


OUTFILE = config.figure_directory + OUTFILE

io_directory = config.clone_data_directory+'simulated_wgs_data/'

cohort_size = 5
gene_conversion_rate = 10**-2

num_sims = 4

fig,ax = plt.subplots(num_sims,2,figsize = (fig_width,num_sims*fig_height))

numpy.random.seed(seed = 1)

for it_sim in range(0,num_sims):

	popindex = -1

	for population in config.populations:
		wgs_sequencing_directory = config.wgs_data_directory + 'cluster_output/'
		depth_timecourse_filename = '%s%s_depth_timecourse.bz2' % (wgs_sequencing_directory, population)
		pre_annotated_filename = '%s%s_simulated_pre_annotated.txt' % (io_directory,population)
		annotated_timecourse_filename = '%s%s_simulated_annotated.txt' % (io_directory, population)

		depth_line = subprocess.check_output("bzcat %s" % depth_timecourse_filename, shell=True)
		median_sequencing_depth = [float(item) for item in depth_line.strip().split(', ')[-1].split(' ')]
		
		text_depth_filename = '%s%s_depth_timecourse.txt' % (io_directory, population)
			
		with open(text_depth_filename, 'w') as f: 
			f.write(depth_line)

		old_stdout = sys.stdout
		sys.stdout = open(pre_annotated_filename, 'w')
		popindex += 1

		max_barcode = config.max_barcode[population]

		clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)

		population_tree = {"":population_tree}
		
		ax[it_sim][popindex].set_xticks([110*i for i in range(0,max_barcode+1)])
		ax[it_sim][popindex].set_xlim([0,110*max_barcode-10])
		ax[it_sim][popindex].set_ylim([0,1.])


		sys.stderr.write("subsampling...\n")

		subsampled_trajectories = {}
		sequencing_times = [110*i-10 for i in range(1,max_barcode+1)]
		num_printed_lines = 0
		num_rejected = 0

		for clone_ID in clone_list:

			total_mutations = 1 + numpy.random.poisson(cohort_size)

			for num in range(0,total_mutations):

				mut_ID = clone_ID + ";"+str(num)

				depth_timecourse = numpy.zeros(len(sequencing_times))
				alternate_timecourse = numpy.zeros(len(sequencing_times))

				freqs  = numpy.asarray([ lineage_dict[clone_ID].freqs[t/10] for t in sequencing_times])
				converted = numpy.random.binomial(1,gene_conversion_rate)

				for i in range(0,len(sequencing_times)):

					sequencing_depth = numpy.random.poisson(median_sequencing_depth[i])
					if converted:
						mutated_chromosome_reads = sequencing_depth
					else:
						mutated_chromosome_reads = numpy.random.binomial(sequencing_depth,0.5)

					depth_timecourse[i] = sequencing_depth

					if freqs[i]> 1e-5:
						number_of_draws = numpy.random.binomial(mutated_chromosome_reads, freqs[i])
					else:
						number_of_draws = 0

					alternate_timecourse[i] = number_of_draws
				
				# check if passes basic quality filters	
				if ((alternate_timecourse>=2).sum() > 1) and ((alternate_timecourse >= 2)*(depth_timecourse>=10)*(alternate_timecourse >= 0.05*depth_timecourse)).sum() > 0:
					#f it does, print to file
					print ", ".join([clone_ID, str(num), '', " ".join(str(t) for t in sequencing_times), " ".join(str(a) for a in alternate_timecourse), " ".join(str(d) for d in depth_timecourse)])
					num_printed_lines+=1

				else:
					num_rejected += 1

		sys.stderr.write("%d passed, %d rejected\n" % (num_printed_lines, num_rejected))

		sys.stdout = old_stdout

		# Annotate pvalues
		sys.stdout.write("Calculating pvalues...\n")
		return_val = os.system('cat %s %s | ./annotate_pvalues > %s' % (text_depth_filename,  pre_annotated_filename, annotated_timecourse_filename))
		if return_val==0:
			sys.stderr.write('Done!\n')
		else:
			sys.stderr.write("Error!\n")

		#now plot

		for line in open(annotated_timecourse_filename,'r'):
			items = line.strip().split(', ')

			lineage_ID = items[0]
			mut_ID = items[1]

			times = numpy.asarray([float(item.strip()) for item in items[3].split(' ')],dtype = float) 

			counts = numpy.asarray([float(item.strip()) for item in items[4].split(' ')],dtype = float) 
			depths = numpy.asarray([float(item.strip()) for item in items[5].split(' ')],dtype = long)

			pvalues = numpy.asarray([float(item.strip()) for item in items[6].split(' ')],dtype = float)

			if pvalues[1] < 0.05:

				R,G,B,A= clone_dict[lineage_ID].color
				line_color = R,G,B,1.

				freq_traj = counts/depths
				ax[it_sim][popindex].plot(sequencing_times, freq_traj, color = line_color, lw = 1.5)
		
		if it_sim == 0:
			ax[it_sim][popindex].text(0.,1.05,config.pop_labels[population],fontsize = 10,weight = 'bold')

		line1,= ax[it_sim][popindex].plot([0,1100],[0.5,0.5],color = 'k',lw = 0.5)
		line1.set_dashes((1,1))

		ax[it_sim][popindex].set_xticks([0]+[110*x-10 for x in range(1,max_barcode+1)])
		ax[it_sim][popindex].set_xticklabels(["1.0"]+["%d.%d" % (epoch+1,100) for epoch in range(0,max_barcode)])
		ax[it_sim][popindex].set_xlim([0,110*max_barcode])
		for tick in ax[it_sim][popindex].get_xticklabels():
	   		tick.set_rotation(90)
	   		tick.set_horizontalalignment('center')
	   	ax[it_sim][popindex].tick_params(axis='x', which='major', pad=3)
		ax[it_sim][popindex].set_ylabel('Allele frequency',fontsize = fontsize)
		ax[it_sim][popindex].set_xlabel('Generations',fontsize = fontsize)
		ax[it_sim][popindex].spines['right'].set_visible(False)
		ax[it_sim][popindex].spines['top'].set_visible(False)
		# Only show ticks on the left and bottom spines
		ax[it_sim][popindex].yaxis.set_ticks_position('left')
		ax[it_sim][popindex].xaxis.set_ticks_position('bottom')


		fig.savefig(config.figure_directory+'si/simulated_trafic.pdf', bbox_inches = 'tight')

