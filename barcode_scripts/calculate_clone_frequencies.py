import sys, os, glob, csv, re
import string, math, numpy
import argparse
from copy import deepcopy

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import lineage.file_parser as file_parser

def insert_clone(barcode,parental_barcode, clone_dict, full_dataset):

	if parental_barcode in clone_dict.keys():
		# there is a parental lineage with this barcode: insert barcode into dict
		epoch = len(barcode.split("_")) - 1
		p_ID = "_".join(barcode.split("_")[:-1])
		last_ID = barcode.split("_")[-1]
		clone_dict[barcode] = deepcopy(full_dataset[epoch][p_ID][last_ID])

		# then, reduce frequency of parent clone by frequency of this child lineage
		parent_epoch = len(parental_barcode.split("_")) - 1
		difference = epoch - parent_epoch
		clone_dict[parental_barcode].freqs -= clone_dict[barcode].freqs
	else:
		if len(parental_barcode.split("_"))>1:
			# we may not have tried the right parent barcode. look for an earlier parent before inserting lineage
			insert_clone(barcode,"_".join(parental_barcode.split("_")[:-1]),clone_dict,full_dataset)
		else:
			# parent not found. simply insert clone into dict
			epoch = len(barcode.split("_")) - 1
			p_ID = "_".join(barcode.split("_")[:-1])
			last_ID = barcode.split("_")[-1]
			clone_dict[barcode] = deepcopy(full_dataset[epoch][p_ID][last_ID])


for population in config.populations:

	timepoints, data, counts = file_parser.get_data(population,config.barcode_data_root_directory,as_matrix = False)

	clone_dict = {}
	clone_list = []

	clone_list_files = [config.clone_data_directory+population+'_clone_list.tsv']
	for filename in clone_list_files:
		f = open(filename,'r')
		for row in csv.reader(f,delimiter = '\t'):
			barcodes = row[0]
			clone_list.append(barcodes)
			parental_ID = "_".join(barcodes.split("_")[:-1])
			insert_clone(barcodes,parental_ID,clone_dict,data)
		f.close()

	with open(config.clone_data_directory+'%s-clone_frequencies.tsv' % (population), 'w') as csvfile:
		out_writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
   		out_writer.writerow(['POP']+['BC'] + ['%d' % (10*d) for d in xrange(0,110)])

 		for clone_ID in clone_list:
 			row = [population]+[clone_ID]
 			row.extend(clone_dict[clone_ID].freqs.tolist())
 			out_writer.writerow(row)



