import sys, os, glob
import numpy
import argparse 

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config

parser = argparse.ArgumentParser()
parser.add_argument("-pvalue_threshold", default = 0.05, type = float, help = 'max combined pvalue')


args = parser.parse_args()

input_directory = config.wgs_data_directory + 'annotated/'

pvalue_threshold = args.pvalue_threshold

old_stdout = sys.stdout

for population in config.populations:
	filenames =  glob.glob(input_directory+population+'*.txt')
	filtered_filename = config.wgs_data_directory + 'filtered/%s_filtered_timecourse.txt' % population
	sys.stdout = open(filtered_filename, 'w')

	for filename in filenames:
		f = open(filename,'r')
		
		#print header
		print ", ".join(['Chromosome','Position','Allele','Timepoints','Alt reads','Depths','Autocorrelation Combined_Pvalue'])

		passed_lines = 0
		all_lines = 0

		for line in f:
			all_lines += 1

			items = line.split(', ')
			
			counts = numpy.asarray([int(item.strip()) for item in items[4].split(' ')],dtype = long) 
			depths = numpy.asarray([float(item.strip()) for item in items[5].split(' ')],dtype = long)
			pvalues = numpy.asarray([float(item.strip()) for item in items[6].split(' ')],dtype = float)


			if (pvalues[1] < pvalue_threshold):
				passed_lines+=1
				print line,


		f.close()
		sys.stderr.write("%d trajectories passed filter\n" % (passed_lines))
sys.stdout = old_stdout

