import sys, os, glob
import numpy

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config

for population in config.populations:

	max_barcode = config.max_barcode[population]
	input_directory = config.wgs_data_directory + 'cluster_output/'
	output_directory = config.wgs_data_directory + 'annotated/'

	merged_timecourse_filename = '%s%s_merged_timecourse.bz2' % (input_directory, population)
	snp_timecourse_filename = '%s%s_snp_timecourse.bz2' % (input_directory, population)
	indel_timecourse_filename = '%s%s_indel_timecourse.bz2' % (input_directory, population)
	depth_timecourse_filename = '%s%s_depth_timecourse.bz2' % (input_directory, population)

	likelihood_timecourse_filename = '%s%s_timecourse.txt' % (output_directory, population)

	sys.stdout.write("\nProcessing %s...\n" % population)

	# Filter SNPs and calculate avg depth per sample
	sys.stdout.write('Filtering SNPS and calculating depth...\n')
	return_val = os.system('python filter_snps_and_calculate_depth.py %s %s %s %s %d' % (merged_timecourse_filename, depth_timecourse_filename, snp_timecourse_filename, indel_timecourse_filename,max_barcode))
	if return_val==0:
		sys.stdout.write('Done!\n')
	else:
		sys.stdout.write("Error!\n")

	os.system('bzcat %s \n' % (depth_timecourse_filename))

	# Annotate pvalues
	sys.stdout.write("Calculating pvalues...\n")
	return_val = os.system('bzcat %s %s %s | ./annotate_pvalues > %s' % (depth_timecourse_filename, snp_timecourse_filename, indel_timecourse_filename, likelihood_timecourse_filename))
	if return_val==0:
		sys.stderr.write('Done!\n')
	else:
		sys.stderr.write("Error!\n")
