
import numpy, itertools
from copy import deepcopy

def insert_into_tree(barcode,tree):
	ancestral_barcode_options = ["_".join(barcode.split("_")[0:l]) for l in range (1,len(barcode.split("_")))]
	ancestors_in_root = [item for item in ancestral_barcode_options if item in tree.keys()]
	if len(ancestors_in_root) == 1:
		ancestral_barcode = ancestors_in_root[0]
		insert_into_tree(barcode,tree[ancestral_barcode])
	elif len(ancestors_in_root) == 0:
		tree.update({barcode:{}})
	else:
		print "Found more than one ancestor in root!"

def child_list(parent_ID, tree):
	candidate = next(item for item in tree.keys() if parent_ID.startswith(item))
	if candidate == parent_ID:
		return tree[candidate].keys()
	else:
		return child_list(parent_ID, tree[candidate])

def ancestor_list(child_ID, tree):
	ancestor = next(item for item in tree.keys() if child_ID.startswith(item))
	if ancestor == child_ID:
		return []
	else:
		return itertools.chain([ancestor],ancestor_list(child_ID,tree[ancestor]))

def compile_max_freqs(key_list,dictionary):
	max_freqs = numpy.zeros(len(key_list))
	it = 0
	for key in key_list:
		max_freqs[it] = max(dictionary[key].freqs)
		it+=1
	return max_freqs

def assign_muller_lower(sorted_list_of_indices,lineage_dictionary,population_tree,baseline,range_available):
	
	middle_item = sorted_list_of_indices[0]
	list_below = deepcopy(sorted_list_of_indices[1::2])
	list_above = deepcopy(sorted_list_of_indices[2::2])

	freq_below,freq_above = numpy.zeros(len(baseline)), numpy.zeros(len(baseline))
	for ID in list_below:
		freq_below += lineage_dictionary[ID].freqs
	for ID in list_above:
		freq_above += lineage_dictionary[ID].freqs

	middle_item_midline = 0.5*range_available

	space_below = middle_item_midline - lineage_dictionary[middle_item].freqs/2. - freq_below
	if any(space_below < 0.):
		middle_item_midline[space_below < 0.] -= space_below[space_below < 0.]

	space_above = (range_available - middle_item_midline) - lineage_dictionary[middle_item].freqs/2. - freq_above
	if any(space_above < 0.):
		middle_item_midline[space_above < 0.] += (space_above[space_above < 0.])

	lineage_dictionary[middle_item].muller_lower = baseline + middle_item_midline - lineage_dictionary[middle_item].freqs/2.
	
	if len(child_list(middle_item,population_tree)) >0:

		max_freqs = compile_max_freqs(child_list(middle_item,population_tree),lineage_dictionary)
		sorted_child_ids = [x for y,x in reversed(sorted( zip(max_freqs,child_list(middle_item,population_tree)) )) ]

		assign_muller_lower(sorted_child_ids,lineage_dictionary,population_tree,
			baseline = lineage_dictionary[middle_item].muller_lower, 
				range_available = lineage_dictionary[middle_item].freqs)

	if len(list_below) > 0:
		assign_muller_lower(list_below,lineage_dictionary,population_tree,
			baseline = baseline,
				range_available=middle_item_midline - lineage_dictionary[middle_item].freqs/2.)

	if len(list_above) > 0:
		assign_muller_lower(list_above,lineage_dictionary,population_tree,
			baseline = baseline + middle_item_midline + lineage_dictionary[middle_item].freqs/2.,
				range_available=range_available - middle_item_midline - lineage_dictionary[middle_item].freqs/2.)
