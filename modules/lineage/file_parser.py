import sys, os, glob, csv, re
import string, math, numpy
from itertools import takewhile

from lineage import *

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def add_nested_dict_item(dict_item,dictionary):
	for item_key, item_value in dict_item.iteritems():
		break
	if item_key in dictionary:
		if (isinstance(item_value,dict)):
			add_nested_dict_item(item_value,dictionary[item_key])
		else:			# found duplicate, display warning
			sys.stderr.write("Warning: found duplicate ---> ignoring item (last barcode is %s) \n" % item_key)
	else:
		dictionary[item_key] = item_value

def parse_frequency_file(filename,prepend_zeros = True):
	######### parses file and returns dictionary of the form {parentID:{lastBarcode:lineage_data}} #######
	######		 NOTE: skips first (population) population barcode. #######
	f = open(filename,'r')
	num_bcd_read = False
	dictionary = {}
	for row in csv.reader(f,delimiter='\t'):
		if not num_bcd_read:
			no_barcodes = 0
			for bcd in takewhile(lambda x: not x.isdigit(), row):
				no_barcodes += 1
			timepoints = numpy.asarray(row[no_barcodes:],dtype = float)
			num_bcd_read = True
		else:
			
			if prepend_zeros:
				frequencies = numpy.zeros(11*(no_barcodes-2))
				frequencies = numpy.concatenate((frequencies,numpy.asarray(row[no_barcodes:],dtype = float)))
			else:
				frequencies = numpy.asarray(row[no_barcodes:],dtype = float)
			dict_item = {"_".join(row[1:no_barcodes-1]): \
							{row[no_barcodes-1]: \
								Lineage('_'.join(row[1:no_barcodes]),frequencies)}}
			add_nested_dict_item(dict_item,dictionary)

	f.close()
		
	return timepoints, dictionary

def parse_frequency_file_as_matrix(filename):
	# find number of barcodes
	f = open(filename,'r')
	num_bcd_read = False
	for row in csv.reader(f,delimiter='\t'):
		if not num_bcd_read:
			no_barcodes = 0
			for bcd in takewhile(lambda x: not x.isdigit(), row):
				no_barcodes += 1
			timepoints = numpy.asarray(row[no_barcodes:],dtype = float)
			num_bcd_read = True
		else:
			pass
	f.close()
	usecols = tuple(numpy.arange(len(timepoints)) + no_barcodes)
	freqs = np.loadtxt(filename, dtype=float, delimiter='\t', skiprows=1, usecols=usecols)
	return timepoints, freqs
		

def parse_count_file(filename):
	f = open(filename,'r')
	count_dict = {}
	for row in csv.reader(f):
		row =  row[0].split()
		count_dict[int(row[0])] = int(row[1])
	return count_dict
	f.close()


def get_data(population,root_path,fitness_files = None,as_matrix = False):
	# input location 
	directory = root_path + population +'/'
	
	# check that this directory exists
	if not os.path.exists(directory):
		sys.stderr.write("Error: directory not found \n \t attempted: "+directory+'\n')
		return [], {}, []
	else:
		sys.stderr.write("processing population %s...\n" % population)
		sys.stderr.write("\treading files...")

		freq_files = natural_sort(glob.glob(directory+population+'*BC*_frequencies.txt'))
		count_files = glob.glob(directory+population+'*read_coverage.txt')		

		# read total read depth from file
		for filename in count_files:
			count_dict = parse_count_file(filename)
		counts = numpy.asarray([count_dict[key] for key in sorted(count_dict.keys())])

		if not as_matrix:
			# read barcodes and their frequencies and return a dictionary of lineages
			data, timepoints = [], []
			for filename in freq_files:
				times, dataset = parse_frequency_file(filename)
				data.append(dataset)
				timepoints.append(times)
			sys.stderr.write("\tcompleted!\n")

			if fitness_files is not None:
				sys.stderr.write(" reading fitness inferences...\n")
				for filename in fitness_files:
					f = open(filename,'r')
					f.readline()
					#scale all fitnesses by number of generations in environment
					if "evolution" in filename:
						key = "evolution"
						scale = 100.
					else:
						key = "barcoding"
						scale = 1.

					for row in csv.reader(f,delimiter='\t'):
						bcd = row[0]
						fitnesses = numpy.asarray(row[1:],dtype = float)
						bcd_order = len(bcd.split('_'))
						parental_barcode = "_".join(bcd.split("_")[:-1])

						if not isinstance(data[bcd_order - 1][parental_barcode][bcd.split("_")[-1]].relative_fitness,dict):
							zero_vector = numpy.zeros(len(data[bcd_order - 1][parental_barcode][bcd.split("_")[-1]].relative_fitness))
							data[bcd_order - 1][parental_barcode][bcd.split("_")[-1]].relative_fitness = {"evolution":zero_vector,
																											"barcoding":zero_vector}
							data[bcd_order - 1][parental_barcode][bcd.split("_")[-1]].relative_fitness_CI = {"evolution":zip(zero_vector,zero_vector),
																											"barcoding":zip(zero_vector,zero_vector)}																					
						data[bcd_order - 1][parental_barcode][bcd.split("_")[-1]].relative_fitness.update({key:fitnesses[::3]*scale})
						data[bcd_order - 1][parental_barcode][bcd.split("_")[-1]].relative_fitness_CI.update({key:zip(fitnesses[1::3]*scale,fitnesses[2::3]*scale)})
					f.close()
			return timepoints, data, counts

		else:

			freq_matrix, timepoints = [], []
			for filename in freq_files:
				times, freqs = parse_frequency_file_as_matrix(filename)
				timepoints.append(times)
				freq_matrix.append(freqs)
			sys.stderr.write("\tcompleted!\n")
			return timepoints, freq_matrix, counts

def read_kappas_from_file(kappa_filename):
	kappa_file = open(kappa_filename,'r')
	kappa_file.readline()	
	kappas = np.asarray([float(item) for item in kappa_file.readline().strip().split('\t')], dtype = float)
	return kappas
	

def read_empirical_null_from_file(empirical_distribution_filename):
	empirical_distribution_file = open(empirical_distribution_filename,'r')
	
	empirical_null = np.asarray([float(item) for item in empirical_distribution_file.readline().strip().split('\t')], dtype = float)
	num_all = len(empirical_null)

	for t in np.sort(empirical_null):
		if ((empirical_null >= t).sum()*1.+1.)/(num_all+1.) <=0.025:
			break

	t_statistic_95_percent_cutoff = t

	return empirical_null, t_statistic_95_percent_cutoff

def read_q_values_from_file(q_value_filename):
	q_value_file = open(q_value_filename,'r')
	q_values = np.asarray([float(item) for item in q_value_file.readline().strip().split('\t')], dtype = float)

	return q_values



