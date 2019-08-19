import os
import os.path

import numpy as np
from datetime import datetime

from bgcore import tsv
from oncodrivefm import tdm

from oncodrivefm.analysis import OncodriveFmAnalysis
from oncodrivefm.method import create_method, method_names, PVALUE_EPSILON
from oncodrivefm.mapping import MatrixMapping

from base import Command

class CombineCommand(Command):
	def __init__(self):
		Command.__init__(self, prog="oncodrivefm-combine", desc="Combine FM bias results")

	def _add_arguments(self, parser):
		Command._add_arguments(self, parser)

		parser.add_argument("data_paths", metavar="DATA", nargs="+",
							help="Files with the results to be combined")

		parser.add_argument("-m", dest="method", metavar="NAME",
							choices=method_names(),
							help="The NAME of the method to combine values")

		parser.add_argument("--save-data", dest="save_data", default=False, action="store_true",
							help="The input data matrix will be saved")

	def _check_args(self):
		Command._check_args(self)

		if self.args.analysis_name is None:
			names = [os.path.splitext(os.path.basename(path))[0] for path in self.args.data_paths]
			prefix = os.path.commonprefix(names)
			if prefix.endswith("-"):
				prefix = prefix[:-1]
			if len(prefix) == 0:
				prefix = "oncodrivefm"
			self.args.analysis_name = prefix

	def load_data(self, data_paths, method=None):
		columns = []
		col_names = []
		row_name_index = {}
		for col_index, data_file in enumerate(data_paths):
			self.log.debug("  > {0}".format(data_file))
			names = []
			values = []
			with tsv.open(data_file, "r") as f:
				col_name, ext = os.path.splitext(os.path.basename(data_file))
				params = tsv.params(f)
				if "slice" in params:
					col_name = params["slice"]
				if "method" in params:
					if method is None:
						method = params["method"]
					elif method != params["method"]:
						self.log.warn("Different method of computation used for file {0}".format(data_file))

				for name, value in tsv.lines(f, (str, float), header=True, null_value="-"):
					if len(name) == 0:
						self.log.warn("Empty identifier detected")
						continue

					if name not in row_name_index:
						row_name_index[name] = len(row_name_index)

					names += [name]
					values += [value]
			col_names += [col_name]
			columns += [(names, values)]

		num_cols = len(columns)
		num_rows = len(row_name_index)
		row_names = [None] * num_rows
		for name, index in row_name_index.items():
			row_names[index] = name

		data = np.empty((num_rows, num_cols))
		data[:] = np.nan

		for col_index, (names, values) in enumerate(columns):
			for i, name in enumerate(names):
				data[row_name_index[name], col_index] = values[i]

		return row_names, col_names, data, method

	def run(self):
		Command.run(self)

		# Load data

		self.log.info("Loading data ...")

		#TODO: Allow to specify the name of the column to load from data files: --data-column=PVALUE && /file.tsv,column=PVALUE
		#TODO: Allow TDM format

		row_names, col_names, data, method = self.load_data(self.args.data_paths, self.args.method)

		self.log.debug("  {0} rows, {1} columns to combine with method '{2}'".format(
							len(row_names), len(col_names), method or "unknown"))

		if method is None:
			self.log.error("Method of combination not defined. Use -m to define it.")
			exit(-1)

		method = create_method(method)

		if self.args.save_data:
			self.log.info("Saving data matrix ...")
			self.save_matrix(self.args.output_path, self.args.analysis_name, self.args.output_format,
						   row_names, col_names, data, suffix="data")

		self.log.info("Combining data using method '{0}' ...".format(method.name))

		combined_results = method.combine(np.ma.masked_invalid(data))

		self.log.info("Saving combined results ...")
		self.save_matrix(self.args.output_path, self.args.analysis_name, self.args.output_format,
						 row_names, method.combination_columns, combined_results.T,
						 params=[("slices", ",".join(col_names)), ("method", method.name)],
						 valid_row=lambda row: sum([1 if np.isnan(v) else 0 for v in row]) == 0)

def main():
	CombineCommand().run()

if __name__ == "__main__":
	main()