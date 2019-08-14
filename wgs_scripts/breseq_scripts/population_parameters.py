import numpy
import sys
from math import log10,exp,ceil

output_directory = 'data/breseq_output'

populations = ['C1','D1']

chrom_dict = {'chr01':230218, 'chr02':808695, 'chr03':306042, 'chr04':1531935, 'chr05':575777, 'chr06':270161, 'chr07': 1090941, 'chr08': 562641,'chr09': 439891, 'chr10': 745752,'chr11': 666818, 'chr12':  1078176,'chr13':  924758,'chr14': 784328, 'chr15': 1091101, 'chr16': 948063 }

#initialize sample names, times, and fastq directories
sample_names = {}
sample_times = {}

#initialize filenames and clone accessions
population_filenames = set()

# read entire file
population_file = open("ilt_wgs_population_samples.csv","r")
population_file.readline() # skip header
for line in population_file:
    items = line.split(",")
    sample_name = items[3].strip()
    fastq_directory = items[4].strip()
    population_filenames.add(sample_name) # create set of filenames for this population
    populations = [subitem.strip() for subitem in items[0].split(";")] # if file relevant to multiple pops (e.g. 0-timepoint) split these
    timepoint = float(items[1])*100 # timestamp
    for population in populations: # compile fastq directories
        if population not in sample_names: 
            sample_names[population] = []
            sample_times[population] = []
        sample_names[population].append((sample_name,fastq_directory))
        sample_times[population].append(timepoint)
population_file.close()

population_filenames = list(population_filenames)
population_filenames.sort()

# sort by timepoint
for population in sample_names.keys():
    sample_times[population], sample_names[population] = (list(x) for x in zip(*sorted(zip(sample_times[population],sample_names[population]))))


def print_parameters():
   for sample in sample_list:
      print sample

def parameters_to_string(parameters):
   return "_".join([str(param) for param in parameters])

def print_split_parameters(param_string):
    for item in param_string.split("_"):
        print item,


if __name__ == '__main__':
    if True:
        if sys.argv[1] == "samples_directories":
            populations = [item.strip() for item in sys.argv[2].split(",")]

            for population in populations:
                desired_samples = sample_names[population]

                for (sample_name, directory) in desired_samples:
                    print sample_name, directory

        if sys.argv[1] == "samples":
            populations = [item.strip() for item in sys.argv[2].split(",")]

            for population in populations:
                desired_samples = sample_names[population]

                for (sample_name, directory) in desired_samples:
                    print sample_name


        elif sys.argv[1] == 'times':
            population = sys.argv[2]
            for t in sample_times[population]:
                print t,
            print ""
        elif sys.argv[1] == 'directory':
            print output_directory
        elif sys.argv[1] == "split":
            print_split_parameters(sys.argv[2])
        elif sys.argv[1] == "population_filenames":
            for filename in population_filenames:
                print filename
        elif sys.argv[1] == "chromosomes":
            round_at = long(sys.argv[2])
            for chromosome in sorted(chrom_dict.keys()):
                print chromosome, int(ceil(float(chrom_dict[chromosome])/round_at)*round_at)
        else:
	    print "Usage error!"
   