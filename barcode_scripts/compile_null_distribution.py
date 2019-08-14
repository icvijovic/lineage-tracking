import sys, os, glob, csv, re
import string, math, numpy
from copy import deepcopy
import matplotlib.gridspec as gridspec
from scipy import stats

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import lineage.inference_params as inference_params
import lineage.file_parser as file_parser
from lineage.fitness_estimator import *

for population in config.populations:

	timepoints, data, counts = file_parser.get_data(population, config.barcode_data_root_directory)
	kappas = file_parser.read_kappas_from_file(config.error_model_directory + population + '-kappas.tsv')

	fitness_estimator = FitnessEstimator(counts,kappas)

	max_barcode = config.max_barcode[population]	

	for environment in ['barcoding', 'evolution']:
		barcoding = environment == 'barcoding'
		q_values, empirical_null, t_statistic_95_percent_cutoff = determine_q_values(data,fitness_estimator,max_barcode - barcoding,barcoding = barcoding)

		with open(config.error_model_directory+'%s-q_values_%s.tsv' % (population,environment), 'w') as csvfile:
			  	out_writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			   	out_writer.writerow(["%f" %x for x in q_values])
			   	out_writer.writerow([inference_params.FDR])

		with open(config.error_model_directory+'%s-empirical_null_%s.tsv' % (population,environment), 'w') as csvfile:
			  	out_writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			   	out_writer.writerow(["%f" %x for x in sorted(empirical_null)])


