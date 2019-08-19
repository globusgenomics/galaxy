import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

class ZscoreTest(object):
	NAME = "zscore"

	def __init__(self):
		self.name = self.NAME
		self.results_columns = ["ZSCORE"] #, "N"]
		self.combination_columns = ["ZSCORE", "PVALUE", "QVALUE"]

	def observed(self, data):
		raise Exception("Abstract method")

	def estimator(self, sample):
		raise Exception("Abstract method")

	def create_results(self, num_slices, num_rows):
		results = np.empty((num_slices, num_rows))
		results[:] = np.nan
		return results

	def compare(self, observed, samples):
		estimated_std = np.std(samples)
		if estimated_std > 0.0:
			estimated_mean = np.mean(samples)
			return (observed - estimated_mean ) / estimated_std
		else:
			return np.nan

	def combine(self, results):
		"""
		Stouffer combination of zscores
		:param results:
		:return:
		"""

		zscores = results.sum(axis=1) / np.sqrt(results.count(axis=1))

		size = zscores.size
		is_nan = zscores.mask
		valid_indices = np.where(~is_nan)
		invalid_indices = np.where(is_nan)

		pv = stats.norm.sf(zscores[valid_indices])

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

		return np.array([zscores, pvalues, qvalues])

class MeanZscoreTest(ZscoreTest):
	NAME = "mean-zscore"

	def observed(self, data):
		return np.ma.mean(data)

	def estimator(self, sample):
		return np.mean(sample)

class MedianZscoreTest(ZscoreTest):
	NAME = "median-zscore"

	def observed(self, data):
		return np.ma.median(data)

	def estimator(self, sample):
		return np.median(sample)
