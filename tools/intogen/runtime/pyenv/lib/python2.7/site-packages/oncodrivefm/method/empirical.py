import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

from base import PVALUE_EPSILON

class EmpiricalTest(object):
	NAME = "empirical"

	def __init__(self):
		self.name = self.NAME
		self.results_columns = ["PVALUE"] #, "N"]
		self.combination_columns = ["PVALUE", "QVALUE"]

	def observed(self, data):
		raise Exception("Abstract method")

	def estimator(self, sample):
		raise Exception("Abstract method")

	def create_results(self, num_slices, num_rows):
		results = np.empty((num_slices, num_rows))
		results[:] = np.nan
		return results

	def compare(self, observed, samples):
		def _count_ge_samples(obs):
			return np.ma.masked_less(samples, obs).count() if not np.isnan(obs) else 0

		count_ge_samples = np.vectorize(_count_ge_samples)

		pvalues = count_ge_samples(observed) / float(samples.size)

		pvalues[pvalues < PVALUE_EPSILON] = PVALUE_EPSILON

		return pvalues

	def combine(self, results):
		"""
		Fisher's combination of pvalues
		:param results:
		:return:
		"""

		results = np.copy(results)
		results[results < PVALUE_EPSILON] = PVALUE_EPSILON

		log = np.ma.log(results)
		s = log.sum(axis=1)
		count = log.count(axis=1)

		size = s.size
		is_nan = s.mask
		valid_indices = np.where(~is_nan)
		invalid_indices = np.where(is_nan)

		pv = 1.0 - stats.chi2.cdf(-2.0 * s[valid_indices], 2 * count[valid_indices])
		pvalues = np.empty(size)
		pvalues[valid_indices] = pv
		pvalues[invalid_indices] = np.nan

		if pv.size != 0:
			qv = multipletests(pv, method='fdr_bh')[1]
		else:
			qv = np.array([])
		qvalues = np.empty(size)
		qvalues[valid_indices] = qv
		qvalues[invalid_indices] = np.nan

		return np.array([pvalues, qvalues])

class MeanEmpiricalTest(EmpiricalTest):
	NAME = "mean-empirical"

	def observed(self, data):
		return np.ma.mean(data)

	def estimator(self, sample):
		return np.mean(sample)

class MedianEmpiricalTest(EmpiricalTest):
	NAME = "median-empirical"

	def observed(self, data):
		return np.ma.median(data)

	def estimator(self, sample):
		return np.median(sample)