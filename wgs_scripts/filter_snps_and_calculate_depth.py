import numpy
import sys
import bz2

import produce_annotation_map

input_filename = sys.argv[1]
depth_filename = sys.argv[2]
snp_filename = sys.argv[3]
indel_filename = sys.argv[4]
max_barcode = int(sys.argv[5])

input_file = bz2.BZ2File(input_filename,"r")
snp_file = bz2.BZ2File(snp_filename,"w")
depth_file = bz2.BZ2File(depth_filename,"w")
indel_file = bz2.BZ2File(indel_filename,"w")

avg_depths = None
times = None
alts = None

depth_records = []

num_rejected = 0
num_all = 0
num_indels = 0

for line in input_file:
    num_all += 1

    items = line.split(", ")
    chromosome = items[0]
    position = int(items[1])
    allele = items[2]

    # calculate depths and add them 
    times = numpy.array([float(subitem) for subitem in items[3].split()])   
    alts = numpy.array([float(subitem) for subitem in items[4].split()])
    depths = [float(subitem) for subitem in items[5].split()]

    # discard all data following epoch max_barcode 
    times = times[:max_barcode]
    alts = alts[:max_barcode]
    depths = depths[:max_barcode]

    #print times in new format
    times_array = " ".join([str(int(t)/100)+".100"for t in times])
    alts_array = " ".join([str(int(a)) for a in alts])
    depths_array = " ".join([str(int(d)) for d in depths])

    new_line = ", ".join(items[0:3]+[times_array]+[alts_array]+[depths_array]) + '\n'

    if produce_annotation_map.in_repeat(chromosome,position):
        num_rejected += 1
        continue # in repeat region

    if allele[:5] == 'indel':
        num_indels += 1
        indel_file.write(new_line)
        continue

    if allele[1:3]!='->':
        continue # not a snp!
    
    snp_file.write(new_line)
    
    depth_records.append(depths)
    
depths = numpy.array(depth_records)

# Could do median or mean
#avg_depths = depths.mean(axis=0)
avg_depths = numpy.median(depths, axis=0)

alts = numpy.array([0 for t in times])

depth_line = ", ".join(["BY4742", "0", "Depth", " ".join([str(t) for t in times]), " ".join([str(alt) for alt in alts]), " ".join([str(avg_depth) for avg_depth in avg_depths])])
depth_file.write(depth_line)
depth_file.write("\n")

sys.stderr.write("%d lines processed, %d of them indels, %d of all calls in repeat regions \n" %(num_all,num_indels, num_rejected))

input_file.close()
snp_file.close()
indel_file.close()
depth_file.close()