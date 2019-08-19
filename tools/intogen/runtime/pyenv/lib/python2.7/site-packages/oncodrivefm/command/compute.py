import os
import os.path

import numpy as np
from datetime import datetime

from bgcore import tsv
from bgcore.labelfilter import LabelFilter

from oncodrivefm import tdm
from oncodrivefm.analysis import OncodriveFmAnalysis
from oncodrivefm.method.factory import create_method
from oncodrivefm.method.empirical import EmpiricalTest
from oncodrivefm.method.zscore import ZscoreTest
from oncodrivefm.mapping import MatrixMapping

from base import Command

class ComputeCommand(Command):
	def __init__(self):
		Command.__init__(self, prog="oncodrivefm-compute", desc="Compute the FM bias")

	def _add_arguments(self, parser):
		Command._add_arguments(self, parser)

		parser.add_argument("data_path", metavar="DATA",
							help="File containing the data matrix in TDM format")

		parser.add_argument("-N", "--samplings", dest="num_samplings", type=int, default=10000, metavar="NUMBER",
							help="Number of samplings to compute the FM bias pvalue")

		parser.add_argument("-t", "--threshold", dest="mut_threshold", type=int, default=2, metavar="THRESHOLD",
							help="Minimum number of mutations to compute the FM bias")

		parser.add_argument("-s", "--slices", dest="slices", metavar="SLICES",
							help="Slices to process separated by commas")

		parser.add_argument("-e", "--estimator", dest="estimator", metavar="ESTIMATOR",
							choices=["mean", "median"], default="mean",
							help="Test estimator for computation.")

		parser.add_argument("-m", "--mapping", dest="mapping", metavar="PATH",
							help="File with mappings between rows and features to be analysed")

		parser.add_argument("-f", "--filter", dest="filter", metavar="PATH",
							help="File containing the features to be filtered. By default labels are includes,"
								 " labels preceded with - are excludes.")

		parser.add_argument("--save-data", dest="save_data", default=False, action="store_true",
							help="The input data matrix will be saved")

	def _check_args(self):
		Command._check_args(self)

		if self.args.analysis_name is None:
			self.args.analysis_name, ext = os.path.splitext(os.path.basename(self.args.data_path))

		if self.args.num_samplings < 1:
			self._error("Number of samplings out of range [2, ..)")

		if self.args.mut_threshold < 1:
			self._error("Minimum number of mutations out of range [1, ..)")

		if self.args.filter is not None:
			if not os.path.exists(self.args.filter):
				self._error("Filter file not found: {0}".format(self.args.filter))

	def run(self):
		Command.run(self)

		# Load data

		self.log.info("Loading data ...")
		self.log.debug("  > {0}".format(self.args.data_path))

		#TODO: Support loading plain matrices: /file.tsv#name=SIFT

		self.matrix = tdm.load_matrix(self.args.data_path)

		self.log.debug("  {0} rows, {1} columns and {2} slices".format(
		self.matrix.num_rows, self.matrix.num_cols, self.matrix.num_slices))

		# Load filter

		self.filter = LabelFilter()
		if self.args.filter is not None:
			self.log.info("Loading filter ...")
			self.log.debug("  > {0}".format(self.args.filter))

			self.filter.load(self.args.filter)

			self.log.debug("  {0} includes, {1} excludes".format(
				self.filter.include_count, self.filter.exclude_count))

		# Load mapping

		if self.args.mapping is not None:
			self.log.info("Loading mapping ...")
			self.log.debug("  > {0}".format(self.args.mapping))

			self.mapping = self.load_mapping(self.matrix, self.args.mapping, self.filter)

			self.log.debug("  {0} features".format(self.mapping.num_groups))

			method_name = "{0}-{1}".format(self.args.estimator, ZscoreTest.NAME)
		else: # One to one mapping
			map = {}
			for row_name in self.matrix.row_names:
				if self.filter.valid(row_name):
					map[row_name] = (row_name,)
			self.mapping = MatrixMapping(self.matrix, map)
			method_name = "{0}-{1}".format(self.args.estimator, EmpiricalTest.NAME)

		# Get selected slice indices

		if self.args.slices is not None:
			slices = []
			for name in self.args.slices.split(","):
				name = name.strip()
				if name not in self.matrix.slice_name_index:
					self.log.warn("Skipping slice not found: {0}".format(name))
					continue
				slices += [self.matrix.slice_name_index[name]]
		else:
			slices = range(self.matrix.num_slices)

		col_names = [self.matrix.slice_names[i] for i in slices]

		if self.args.save_data:
			for i in slices:
				slice_name = self.matrix.slice_names[i]
				self.log.info("Saving {0} data matrix ...".format(slice_name))
				self.save_matrix(self.args.output_path, self.args.analysis_name, self.args.output_format,
								 self.matrix.row_names, self.matrix.col_names, self.matrix.data[i],
								 suffix="data-{0}".format(slice_name))

		# Run the analysis

		self.log.info("Running the analysing using '{0}' ...".format(method_name))

		analysis = OncodriveFmAnalysis(
			"oncodrivefm.compute",
			num_samplings = self.args.num_samplings,
			mut_threshold = self.args.mut_threshold,
			num_cores=self.args.num_cores)

		results = analysis.compute(self.matrix, self.mapping, method_name, slices)

		method = create_method(method_name)

		self.log.info("Saving results ...")

		#TODO: Have an option to save in TDM instead of splited
		self.save_splited_results(
			self.args.output_path, self.args.analysis_name, self.args.output_format,
			self.matrix, self.mapping, method, results, slices)

def main():
	ComputeCommand().run()

if __name__ == "__main__":
	main()
