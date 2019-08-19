import os.path
import logging
import random
import math
import time

from itertools import repeat

import multiprocessing as mp

import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

from mapping import MatrixMapping
from method.factory import create_method

def sampling(params):
	method, data, count = params
	return method.estimator(random.sample(data, count))

class OncodriveFmAnalysis(object):

	def __init__(self, log_name, num_samplings, mut_threshold, num_cores=None):
		self.num_samplings = num_samplings
		self.mut_threshold = mut_threshold
		self.num_cores = num_cores

		self.log = logging.getLogger(log_name)

	def _count_mutations(self, valid_data, slice, mapping):
		"""
		Counts the number of mutations for each group in the given slice of valid_data
		:param valid_data: [slice, gene, sample] masked for nan's
		:param slice: the slice to count
		:param mapping: a MatrixMapping with groups which mutations has to be counted
		:return: [(mut_count, [group_indices])]
		"""

		# Create a dictionary for each slice where each key represents a mutations count over the mutations threshold
		# and each value a list of group indices with that mutations count [{mut_count : [group_indices]}]

		group_counts = {}

		for group_index, group_rows in mapping.group_row_indices.items(): #TODO parallelize

			mut_count = valid_data[slice, group_rows, :].count()

			if mut_count >= self.mut_threshold:
				if mut_count in group_counts:
					group_counts[mut_count] += [group_index]
				else:
					group_counts[mut_count] = [group_index]

		# Convert the group_counts dictionary {mut_count : [group_indices]}
		# into an ordered list of tuples sorted by mut_count [(mut_count, [group_indices])]

		total_group_indices = set()

		mut_counts = []
		for mut_count, group_indices in sorted(group_counts.items(), key=lambda p:p[0]):
			mut_counts += [(mut_count, group_indices)]
			total_group_indices.update(group_indices)

		return mut_counts, sorted(total_group_indices)

	def compute(self, matrix, mapping, method, slices):

		num_slices = len(slices)

		valid_data = np.ma.masked_invalid(matrix.data)

		method = create_method(method)
		if method is None:
			raise Exception("Unknown test method: {0}".format(method))

		results = method.create_results(num_slices, mapping.num_groups)

		pool = mp.Pool(self.num_cores)

		for slice_results_index, slice in enumerate(slices):

			self.log.info("[{0}]".format(matrix.slice_names[slice]))

			self.log.info("  Counting mutations ...")

			mut_counts, group_indices = self._count_mutations(valid_data, slice, mapping)

			self.log.info("  Calculating observed estimator ...")

			observed = np.empty(mapping.num_groups)
			observed[:] = np.nan

			for group_index in group_indices: #TODO parallelize
				group_rows = mapping.group_row_indices[group_index]
				observed[group_index] = method.observed(valid_data[slice, group_rows, :])

			self.log.info("  Bootstrapping with {0} repetitions ...".format(self.num_samplings))

			background = valid_data[slice].compressed().tolist()

			estimations = np.empty(self.num_samplings)

			for mut_count, group_indices in mut_counts:
				self.log.debug("    With {0} mutations -> {1} groups ...".format(mut_count, len(group_indices)))


				#for i in xrange(self.num_samplings):
				#	sample = random.sample(background, mut_count)
				#	sampling_estimations[i] = method.estimator(sample)

				iter = repeat((method, background, mut_count), self.num_samplings)
				res = pool.map_async(sampling, iter)

				while not res.ready():
					try:
						time.sleep(0.1)
					except KeyboardInterrupt:
						pool.terminate()
						raise

				estimations[:] = np.array(res.get())

				results[slice_results_index, group_indices] = method.compare(observed[group_indices], estimations)

		return results

	def combine(self, results, method):
		method = create_method(method)

		return method.combine(results)
