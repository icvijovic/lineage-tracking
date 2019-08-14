import sys, os, glob, csv, re
import string, math, numpy

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import lineage.file_parser as file_parser

from lineage.read_clone_data import *
from lineage.plot_utils_clone import * # contains muller and bar chart class
from lineage.tree_utils import *


with open(config.clone_data_directory+'all_clone_data.csv', 'w') as csvfile:

	out_writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
	for population in config.populations:
		print config.pop_labels[population]

		clone_list, times, clone_dict, lineage_dict, population_tree = read_clone_data(population)
		clone_list = ['']+ clone_list

		if population == 'C1':
			header = ['Population', 'Barcodes']
			header.extend(['EvolutionFitness(percent)','EvolutionFitnessConfidenceInterval(percent)'])
			header.extend(['BarcodingFitness(per barcoding procedure)','BarcodingFitnessConfidenceInterval(per barcoding procedure)'])
			header.extend(["Frequency(t=%d.%d)" % (t/110+1, t%110) for t in times ])

			out_writer.writerow(header)

		for ID in clone_list:
			row = [config.pop_labels[population]]
			lineage = clone_dict[ID]
			if lineage.ID == '':
				row.extend(['ancestor'])
			else:
				row.extend([ID])

			evolution_fitness_info = [str(lineage.evolution_fitness*100),"%.10f,%.10f" % (lineage.evolution_fitness_CI[0]*100,lineage.evolution_fitness_CI[1]*100)]
			barcoding_fitness_info = [str(lineage.barcoding_fitness),"%.10f,%.10f" % (lineage.barcoding_fitness_CI[0],lineage.barcoding_fitness_CI[1])]
			
			row.extend(evolution_fitness_info)
			row.extend(barcoding_fitness_info)

			row.extend((["%.10f" % f for f in lineage.freqs[:11*config.max_barcode[population]]]))
			row.extend((["NA" for f in lineage.freqs[11*config.max_barcode[population]:]]))
			out_writer.writerow(row)


			