import string, math
import numpy as np


class Lineage:

	def __init__(self,lineage_ID,freq_array):
		self.ID = lineage_ID
		self.freqs = freq_array

		self.evolution_fitness = 0.
		self.barcoding_fitness = 0.
		self.evolution_fitness_CI = np.zeros(2)
		self.barcoding_fitness_CI = np.zeros(2)

		self.relative_fitness = np.zeros(len(freq_array)/11)
		self.relative_fitness_CI = [(0,0) for i in xrange(0,len(self.relative_fitness))]

		self.color = None
		self.muller_lower = np.zeros(len(freq_array))
		self.flags = []




	