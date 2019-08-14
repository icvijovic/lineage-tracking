import sys, os
import numpy
import bz2
import subprocess

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config

MIN_DISTANCE_TO_REPEAT = 100

input_filename = config.wgs_data_directory+'genome_feature_files/reference.gff3.bz2'

genome_annotation_file = bz2.BZ2File(input_filename,"r")

repeat_types = ['telomere', 'centromere_DNA_Element_I', 'centromere_DNA_Element_II', 
	'centromere_DNA_Element_III', 'transposable_element_gene','LTR_retrotransposon',
	'centromere', 'long_terminal_repeat','ARS']
genes = ['gene']
coding_sequence = ['CDS']
RNA_types = ['rRNA','ncRNA','snRNA','snoRNA','tRNA']

codon_table = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
	'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
	'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
	'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
	'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
	'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
	'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
	'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
	'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
	'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 
	'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 
	'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
	'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
	'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
	'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 
	'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
}

base_table = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])


gene_map = {}
CDS_map = {}
temp_repeat_map = {}

gene_names = []
start_positions = []
end_positions = []

for line in genome_annotation_file:

	items = line.split('\t')
	
	if len(items) < 9:
		pass
	else:
		seqid = items[0]

		seq_type = items[2]
		start = long(items[3])
		end = long(items[4])
		strand = items[6]
		phase = items[7]

		if seq_type in repeat_types:
			if seqid in temp_repeat_map:
				temp_repeat_map[seqid].append((long(start),long(end)))
			else:
				temp_repeat_map[seqid] = []

		if seq_type in 'gene':
			has_alias = False
			for info_item in items[8].split(';'):
				if 'gene=' in info_item:
					has_alias = True
					gene_name = info_item.split('=')[1]

			if not has_alias:
				for info_item in items[8].split(';'):
					if 'Name=' in info_item:
						gene_name = info_item.split('=')[1]
			
			if gene_name in gene_map.keys():
				pass
			else:
				gene_map[gene_name] = (seqid,start,end)

		if seq_type == 'CDS':
			phase = long(phase)
			for info_item in items[8].split(';'):
					if 'Name=' in info_item:
						CDS_name = info_item.split('=')[1]
			if CDS_name in CDS_map.keys():
				pass
			else:
				CDS_map[CDS_name] = (seqid,start,end,phase,strand)


final_map = {}
for chrom in temp_repeat_map.keys():
	final_map[chrom] = {}
	for pairs in temp_repeat_map[chrom]:
		for position in xrange(pairs[0]-MIN_DISTANCE_TO_REPEAT, pairs[0]+ MIN_DISTANCE_TO_REPEAT):
			final_map[chrom][position] = final_map[chrom].get(position,'repeat')

for gene_name in gene_map.keys():
	chrom, start, end = gene_map[gene_name]
	if chrom == 'chr17':
		pass
	else:
		for position in xrange(start,end):
			final_map[chrom][position] = final_map[chrom].get(position,gene_name)

pos_to_CDS_map = {}
for CDS_name in CDS_map.keys():

	chrom, start, end, phase, strand = CDS_map[CDS_name]

	if chrom in pos_to_CDS_map.keys():
		pass
	else:
		pos_to_CDS_map[chrom] = {}

	for position in xrange(start,end):
		if strand == '+':
			this_position_phase = (phase+(position-start))%3
		else:
			this_position_phase = (phase + (end - position))%3

		pos_to_CDS_map[chrom][position] = (CDS_name,this_position_phase,strand)

def in_repeat(chrom, pos):
	if final_map[chrom].get(pos,'not-annotated') == 'repeat':
		return True
	else:
		return False

def annotate_mutation(mut_string):
	print mut_string
	chrom, pos, allele = mut_string.split(', ')
	pos= long(pos)

	# first check if mutation is in gene
	if 'intergenic' in final_map[chrom].get(pos,'intergenic (%s)'%chrom):
		return final_map[chrom].get(long(pos),'intergenic (%s)'%chrom), 'intergenic'
	else:
		if '->' in allele:
			# this is a SNP, classify as synonymous or nonsynymous
			# pull out base change in this SNP
			original_base, mutant_base = allele.split('->')

			# lookup relevant CDS, phase and strand
			CDS_name, position_phase, strand = pos_to_CDS_map[chrom][long(pos)]
			# find codon start and end. if strand is reverse also translate the base change
			if strand == '+':
				codon_start = pos - position_phase
				codon_end = codon_start + 3
			else:
				codon_start = pos + position_phase
				codon_end = codon_start - 3
				original_base = calculate_reverse_complement_sequence(original_base)
				mutant_base = calculate_reverse_complement_sequence(mutant_base)

			# print "position_phase = ", position_phase
			# now read the codon
			# identify line that this codon will be on among the FASTA lines for this contig
			start_line_no, end_line_no = codon_start/60 + 1, codon_end/60 + 1
			conting_start_string = ">"+chrom
			read_command = 'bzcat %s | grep -A %d \'%s\' | tail -n %d' % (input_filename, max(start_line_no, end_line_no), chrom, 1 + abs(end_line_no - start_line_no))
			p = subprocess.Popen(read_command, stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			line = "".join(out.split('\n'))


			within_line_start_position = min(codon_start,codon_end)%60 - 1*(strand == '+')
			within_line_end_position = within_line_start_position + 3
			
			reference_codon = line[within_line_start_position:within_line_end_position]
			if strand == "-":
				reference_codon = calculate_reverse_complement_sequence(reference_codon)

				# print start_line_no, end_line_no
				# print within_line_start_position
				# print within_line_end_position
				# print line
				# print reference_codon
				# print position_phase, reference_codon, len(reference_codon)

			if reference_codon[position_phase] == original_base:
				mutant_codon = list(reference_codon)
				mutant_codon[position_phase] = mutant_base
				mutant_codon = "".join(mutant_codon)

				print "codon change: ",  reference_codon, "(", codon_table[reference_codon], ") -> ", mutant_codon, "(", codon_table[mutant_codon], ")"
				
				if codon_table[reference_codon] == codon_table[mutant_codon]:
					mutation_type = 'syn'
				else:
					if codon_table[mutant_codon] == '*':
						mutation_type = 'nonsense'
					else:
						mutation_type = 'nonsyn'
				
			else:
				print "annotation error"
				mutation_type = 'null'
			
		else:
			mutation_type = 'indel'

	return final_map[chrom].get(long(pos),'intergenic (%s)'%chrom), mutation_type

if __name__ == '__main__':
	chrom = sys.argv[1]
	position = long(sys.argv[2])
	print annotate_mutation(', '.join([chrom,str(position),'G->T']))
