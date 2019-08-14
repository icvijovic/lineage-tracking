import string, math
import numpy as np
from scipy.optimize import minimize
import random
import sys
from scipy.stats import chi2

from lineage import *
import inference_params

class FitnessEstimator:	

	def __init__(self, counts, kappas, 
					qvals = None, 
					t_statistic_95_percent_cutoff = None, 
					empirical_null = None):
		# initialize population-level observables
		self.total_reads = counts*1.0
		self.ratio_of_read_depths = self.total_reads[1:]/self.total_reads[:-1]
		self.population_fitness = np.zeros(len(self.total_reads))

		# initialize error model parameters
		self.k = kappas
		self.qvals = qvals
		self.t_statistic_95_percent_cutoff = t_statistic_95_percent_cutoff
		self.empirical_null = empirical_null

	def get_interval_endpoints(self, in_epoch, barcoding = False):
		if barcoding:
			begin = in_epoch*inference_params.INTERVALS_PER_EPOCH + inference_params.EVOLUTION_INTERVALS_PER_EPOCH
			end = begin + inference_params.BARCODING_INTERVALS_PER_EPOCH + 1
		else:
			begin = in_epoch * inference_params.INTERVALS_PER_EPOCH
			end = begin + inference_params.EVOLUTION_INTERVALS_PER_EPOCH + 1

		return begin, end

	def get_counts(self, lineage, in_epoch, barcoding = False):
		""" Picks out observed counts in appropriate epoch  """

		begin, end = self.get_interval_endpoints(in_epoch,barcoding)
		observed = lineage.freqs[begin:end] * self.total_reads[begin:end]

		return observed, begin, end

	def get_branching_process_log_likelihoods(self, observed, begin, end, barcoding = False):
		
		if barcoding:
			parameters = inference_params.fitness_grid['barcoding']
			multiplier = inference_params.scale_fitness_per_interval['barcoding']
		else:
			parameters = inference_params.fitness_grid['evolution']
			multiplier = inference_params.scale_fitness_per_interval['evolution']

		expectations = np.ones((len(parameters),len(observed)-1))
		expectations *= self.ratio_of_read_depths[begin:end-1]
		expectations = (expectations.T * np.exp(parameters*multiplier)).T
		expectations *= observed[:-1]
		expectations[expectations < 1] = 1

		#evaluate negative log likelihood
		llhs = np.sum((np.sqrt(expectations) - np.sqrt(observed[1:]))**2/self.k[begin:end-1] 
						+ 0.75* np.log(observed[1:]) - 0.25* np.log(expectations)
						+ 0.5 * np.log(4*np.pi*self.k[begin:end-1])
					, axis = 1)	
		return llhs


	def get_log_likelihoods(self, observed, begin, end, barcoding = False):

		if barcoding:
			parameters = inference_params.fitness_grid['barcoding']
			multiplier = inference_params.scale_fitness_per_interval['barcoding']
		else:
			parameters = inference_params.fitness_grid['evolution']
			multiplier = inference_params.scale_fitness_per_interval['evolution']

		freqs = observed*1.0/self.total_reads[begin:end]

		if max(freqs) < 10**-1:
			expectations = np.ones((len(parameters),len(observed)-1)) 
			expectations *= self.ratio_of_read_depths[begin:end-1] \
							* np.exp(-multiplier*self.population_fitness[begin:end-1])
			expectations = (expectations.T * np.exp(parameters*multiplier)).T
			expectations *= observed[:-1]
			expectations[expectations < 1]  = 1
		else:
			expectations = np.exp(np.outer(multiplier*parameters,np.arange(1,len(observed))) - np.outer(0.*parameters,np.cumsum(freqs[:-1]))) * self.total_reads[begin+1:end]/self.total_reads[begin]
			expectations *= observed[0] * np.exp(-multiplier*np.cumsum(self.population_fitness[begin:end-1]))
			expectations[expectations < 1] = 1

		llhs = np.sum((np.sqrt(expectations) - np.sqrt(observed[1:]))**2/self.k[begin:end-1] + 0.75* np.log(observed[1:]) - 0.25* np.log(expectations) + 0.5 * np.log(4*np.pi*self.k[begin:end-1]), axis = 1)
		
		return llhs


	def estimate_relative_fitness(self, observed, begin, end, barcoding = False):

		if barcoding:
			parameters = inference_params.fitness_grid['barcoding']
			multiplier = inference_params.scale_fitness_per_interval['barcoding']
		else:
			parameters = inference_params.fitness_grid['evolution']
			multiplier = inference_params.scale_fitness_per_interval['evolution']

		freqs = observed*1./self.total_reads[begin:end]

		if max(freqs) < 10**-1:

			a = (freqs[:-1] * np.exp(-multiplier*self.population_fitness[begin:end-1]) * self.total_reads[begin+1:end]/self.k[begin:end-1] ).sum()		
			b = -np.sum(np.sqrt(freqs[:-1]*freqs[1:] * np.exp(-multiplier*self.population_fitness[begin:end-1])) * self.total_reads[begin+1:end]/self.k[begin:end-1] )
			c = -len(freqs)/4.
			y = (-b + np.sqrt(b**2 - 4 * a * c))/2/a

			best_fitness = 2* np.log(y)/multiplier
			two_sigma = 1./np.sqrt((2* a - b)*y)
			
			return [best_fitness], best_fitness - two_sigma/2., best_fitness + two_sigma/2.
			
		else:
			llhs = self.get_log_likelihoods(observed, begin, end, barcoding)
		
			best_fitness = parameters[llhs == min(llhs)]
			fitness_CI_range = llhs - min(llhs) < self.t_statistic_95_percent_cutoff/2.

			return best_fitness, min(parameters[fitness_CI_range]), max(parameters[fitness_CI_range])


	def empirical_pvalue(self, lineage, in_epoch, barcoding = False):
		observed, begin, end = self.get_counts(lineage, in_epoch, barcoding)

		if any(observed < inference_params.threshold_lineage_size):
			# do not attempt to infer fitness of miniscule lineages
			return 2., 0
		else:
			# otherwise calculate llh ratio statistic
			llhs = self.get_branching_process_log_likelihoods(observed, begin, end, barcoding)
			min_llh, neutral_llh = min(llhs), llhs[inference_params.num_fitnesses]

			test_statistic = 2*(neutral_llh - min_llh)
			num_larger = (self.empirical_null >= test_statistic).sum()*1. 

			empirical_pvalue = (num_larger+1.)/(len(self.empirical_null)+1.)

			return empirical_pvalue,test_statistic


	def update_lineage_fitness(self, lineage, in_epoch, barcoding = False, masked_epochs = []):
		pval, t = self.empirical_pvalue(lineage, in_epoch, barcoding)
		if pval <= self.qvals[in_epoch]:
			# infer fitness if there's evidence that the lineage is selected or if it has been selected in a previous epoch
			observed, begin, end = self.get_counts(lineage, in_epoch, barcoding)

			fitness, fitness_CI_lower, fitness_CI_upper= self.estimate_relative_fitness(observed, begin, end, barcoding)

			lineage.relative_fitness[in_epoch] = fitness[len(fitness)/2]
			lineage.relative_fitness_CI[in_epoch] = fitness_CI_lower, fitness_CI_upper
			updated = True
		else:
			lineage.relative_fitness[in_epoch] = 0.0
			updated = False

		return updated, lineage.relative_fitness[in_epoch]
	
	def update_mean_fitness(self,dataset,in_epoch,barcoding = False):
		""""Updates the mean fitness in the experiment in epoch 'in_epoch' based on fitnesses of lineages in 'dataset'."""

		begin, end = self.get_interval_endpoints(in_epoch, barcoding)
		self.population_fitness[begin:end] = np.zeros(end-begin)
		for parent_id in dataset.keys():
			for bcd, lineage in dataset[parent_id].iteritems():
				if lineage.relative_fitness is not None:
					self.population_fitness[begin:end] = self.population_fitness[begin:end]+lineage.freqs[begin:end] * lineage.relative_fitness[in_epoch]
					self.population_fitness[begin:end] -= (self.population_fitness[begin]) * np.ones(end-begin)


def determine_q_values(data, fitness_estimator, max_barcodes, barcoding = False):
	threshold_values = np.zeros(max_barcodes)

	empirical_null = []
	
	for target_epoch in range(0,max_barcodes):
			
		test_statistic = []
		
		for dataset in (range(0,target_epoch+1)):
			for parental_barcode in data[dataset].keys():
				for bcd, lineage in data[dataset][parental_barcode].iteritems():

					observed, begin, end = fitness_estimator.get_counts(lineage, target_epoch, barcoding)

					if all(observed > inference_params.threshold_lineage_size):
						llhs = fitness_estimator.get_branching_process_log_likelihoods(observed, begin, end, barcoding)
						min_llh = min(llhs)
						neutral_llh = llhs[inference_params.num_fitnesses]

						test_statistic.append(2*(neutral_llh - min_llh))
						if dataset == 0 and target_epoch == 0:
							empirical_null.append(2*(neutral_llh - min_llh))

		if target_epoch == 0:
			empirical_null = np.asarray(empirical_null)
			num_all = len(empirical_null) 

		for t in np.sort(empirical_null):
			if ((empirical_null >= t).sum()*1.+1.)/(num_all+1.) <=0.025:
				break
		t_statistic_95_percent_cutoff = t

		test_statistic = np.asarray(test_statistic)

		num_larger = np.asarray([(empirical_null >= t).sum()*1. for t in test_statistic])

		empirical_p_values = (num_larger+1.)/(num_all+1.)

		p_values = empirical_p_values
		candidate_pvalues = sorted(np.unique(p_values),reverse=True)
		num = len(p_values)
		for p in candidate_pvalues:
			if (num*p/(p_values <= p).sum()) < inference_params.FDR:
				break

		threshold_values[target_epoch] = p

	return threshold_values, empirical_null, t_statistic_95_percent_cutoff

