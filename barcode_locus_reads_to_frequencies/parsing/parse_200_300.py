#! /usr/bin/python3 -u

# ##############################################################################
# Performs parsing of the cleaned fastq files to identify the barcodes.
# Written and designed by Jose Rojas Echenique
# Modified by Alex N Nguyen Ba
# rojasechenique@fas.harvard.edu and nnguyenba@fas.harvard.edu
# 2016
#
# Version 1
#
# To run, use clean script on the Fastq file.
# Then run this script on the cleaned fastq file.
#
# LICENCE 
#
# ##############################################################################

import sys
from collections import OrderedDict
sys.path.append('/n/home00/nnguyenba/lib/python')

import regex

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(seq):
	return "".join(complements.get(base, base) for base in reversed(seq))

molecular_barcode = ".{8}"
lineage_barcode = "(.{3})(?e)(?:TA){e<=1}(.{5})(?:AT){e<=1}(.{5})(?:TA){e<=1}(.{3})"

KanP1 = "(?e)(?:GCAATTAACCCTCACTAAAGG)"
HygP3 = "(?:GGATTATTCATACCGTCCCA)"


S5index = [
	"GCGATCTA",
	"ATAGAGAG",
	"AGAGGATA",
	"TCTACTCT",
	"CTCCTTAC",
	"TATGCAGT",
	"TACTCCTT",
	"AGGCTTAG",
	"GAGTAGCC",
	"GTCTGAGG",
	"CGTAAGGA",
	"CCACGCGT",
	"GGAGTTCC",
	"CATGGCCA",
	"AATCTCTC",
	"TAACCGCG",
	"TGGCGGTC",
	"CCATCTTA",
	"ATGTCAAT",
	"AGTTGGCT",
	"ACCTAGTA",
	"AACCGTGA",
	"TCATTACA",
	"CTGACGTG",
	"GAATTCAG",
	"CCGGTACG",
	"CCGTCATC",
	"CGTCTATA",
	"TCAATGAC",
	"AACGATGC",
	"GTCAACCT",
	"CAGTTTCA",
	"TGTGATTG",
	"TTGCATGT",
	"GGCGCGAT",
	"TTAACCGA"]

S5indexDict = {idx:(i+1) for i,idx in enumerate(S5index)}

ilt_col_index = [
	"AGTGATT",
	"GAATCCCC",
	"CCTGGGAAA",
	"TACCTCCCAG",
	"ATTTGTGGTAT",
	"CCCGAGAGATCG",
	"GAGTATA",
	"TCTTAATT",
	"CAGCGCTGG",
	"GTTGACGAAT",
	"TGCATCAGCGC",
	"AGCCATATGAGC"]

def make_corrector(options):
	checkers = [regex.compile("("+o+"){e<=2}") for o in options]
	def corrector(match):
		for (i,c) in enumerate(checkers):
			if c.match(match):
				return i+1
	return corrector

correct_ilt_col = make_corrector(ilt_col_index[0:12])

ilt_row_index = [
	"CCGTAAA",
	"GTAGTCGA",
	"TGTAAGAAT",
	"AGGCCCGCCG",
	"CTCCACGCAAG",
	"GTCTTTGATTAG",
	"TAGAGCC",
	"AGAAATGG"]

correct_ilt_row = make_corrector(ilt_row_index[0:8])


N7index = [
	"TAAGGCGA",
	"CGTACTAG",
	"AGGCAGAA",
	"TCCTGAGC",
	"GGACTCCT",
	"TAGGCATG",
	"CTCTCTAC",
	"CAGAGAGG",
	"GCTACGCT",
	"CGAGGCTG",
	"AAGAGGCA",
	"GTAGAGGA",
	"ATTGTAAT",
	"GATCATTC",
	"ACCGATCG",
	"CCGTTATT",
	"TTCTTCTA",
	"TACCTGAC",
	"AGGACCGC",
	"GTCCGATT",
	"CACGAGTT",
	"CCACGGCC",
	"ACATGTAA",
	"TGTTAACT"]

N7indexDict = {idx:(i+1) for i,idx in enumerate(N7index)}

Lox5171L = "(?:GAAAGGGGCTATAATGTGTACTATACGAAGTTAT)"
Lox5171RL = "(?e)(?:GAAAGGGGCTATAATGTGTACTATAGCCCCTTTC)"
Lox5171RLrc = "(?e)(?:GAAAGGGGCTATAGTACACATTATAGCCCCTTTC)"

pTEFRrc = "(?:GGGGACGAGGCAAGCT)"
LoxP2Lrc = "(?:ATAACTTCGTATAATGTATGCTATAGCCCCTTTC)"
LoxPLRrc = "(?e)(?:GAAAGGGGCTATAATGTATGCTATAGCCCCTTTC)"
Lox5171Lrc = "(?e)(?:ATAACTTCGTATAGTACACATTATAGCCCCTTTC)"

def OR(xs):
	return "(" + "|".join(["(?:"+x+")" for x in xs]) + ")"

BC200_300_R1 = (
	"(" + molecular_barcode + ")" +
	OR(ilt_col_index[0:12]) + "{e<=2}" +
	KanP1 + "{e<=4}" + "(?e)(?:TACT){e<=1}" +
	lineage_barcode +
	"(?e)(?:CGCT){e<=1}" + Lox5171RL + "{e<=4}" + "(?e)(?:TACT){e<=1}" +
	lineage_barcode +
	"(?e)(?:CGCT){e<=1}" )

BC200_300_R2 = (
	"(" + molecular_barcode + ")" +
	OR(ilt_row_index[0:8]) + "{e<=2}" + 
	HygP3 + "{e<=4}" + "(?e)(?:AGCG){e<=1}" +
	lineage_barcode +
	"(?e)(?:AGTA){e<=1}" + Lox5171RLrc + "{e<=4}" + "(?e)(?:AGCG){e<=1}" +
	lineage_barcode +
	"(?e)(?:AGTA){e<=1}" )

BC200_300_R1re = regex.compile(BC200_300_R1)
BC200_300_R2re = regex.compile(BC200_300_R2)


ldiscarded = 0
rdiscarded = 0
nreads = 0
for line in sys.stdin:
	nreads += 1
	row=OrderedDict([
		(x,"") for x in ["illumina_index_1","illumina_index_2","row","col","lineage_barcode_1","lineage_barcode_2","lineage_barcode_3","lineage_barcode_4","mb1","mb2"]
	])
	
	splitline = line.split("	")
	if len(splitline) < 3:
		continue

	row["illumina_index_1"] = N7indexDict[splitline[0]]
	row["illumina_index_2"] = S5indexDict[splitline[1]]
	R1 = splitline[2]
	R2 = splitline[3]
	
	m1 = BC200_300_R1re.match(R1)
	if m1:
		row["mb1"] = m1.groups()[0]
		row["col"] = correct_ilt_col(m1.groups()[1])
		row["lineage_barcode_1"] = m1.groups()[2] + m1.groups()[3] + m1.groups()[4] + m1.groups()[5]
		row["lineage_barcode_2"] = m1.groups()[6] + m1.groups()[7] + m1.groups()[8] + m1.groups()[9]

		m2 = BC200_300_R2re.match(R2)
		if m2:
			row["mb2"] = m2.groups()[0]
			row["row"] = correct_ilt_row(m2.groups()[1])
			row["lineage_barcode_4"] = rc(m2.groups()[2] + m2.groups()[3] + m2.groups()[4] + m2.groups()[5])
			row["lineage_barcode_3"] = rc(m2.groups()[6] + m2.groups()[7] + m2.groups()[8] + m2.groups()[9])
		else:
			rdiscarded+=1
			continue
	else:
		#print(R1)
		ldiscarded+=1
		continue
		
	print(*row.values(), sep="	", file=sys.stdout)
print("Discarded " + str(ldiscarded) + "+" + str(rdiscarded) + " out of " + str(nreads) + " reads.", file=sys.stderr)
