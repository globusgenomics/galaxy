import os
import os.path
import argparse
import logging

import numpy as np
from datetime import datetime

from oncodrivefm import VERSION

from bgcore import tsv

from oncodrivefm.mapping import MatrixMapping

_LOG_LEVELS = {
	"debug" : logging.DEBUG,
	"info" : logging.INFO,
	"warn" : logging.WARN,
	"error" : logging.ERROR,
	"critical" : logging.CRITICAL,
	"notset" : logging.NOTSET }

class Command(object):
	def __init__(self, prog=None, desc=""):
		parser = argparse.ArgumentParser(prog=prog, description=desc)

		self._add_arguments(parser)

		parser.add_argument("-j", "--cores", dest="num_cores", type=int, metavar="CORES",
							help="Number of cores to use for calculations. Default is 0 that means all the available cores")

		parser.add_argument("-D", dest="defines", metavar="KEY=VALUE", action="append",
							help="Define external parameters to be saved in the results")

		parser.add_argument("-L", "--log-level", dest="log_level", metavar="LEVEL", default=None,
							choices=["debug", "info", "warn", "error", "critical", "notset"],
							help="Define log level: debug, info, warn, error, critical, notset")

		self.args = parser.parse_args()

		logging.basicConfig(
			format = "%(asctime)s %(name)s %(levelname)-5s : %(message)s",
			datefmt = "%Y-%m-%d %H:%M:%S")

		if self.args.log_level is None:
			self.args.log_level = "info"
		else:
			self.args.log_level = self.args.log_level.lower()

		self.log = logging.getLogger("oncodrivefm")

		self.log.setLevel(_LOG_LEVELS[self.args.log_level])

		self._arg_errors = []

	def _add_arguments(self, parser):
		parser.add_argument("-o", "--output-path", dest="output_path", metavar="PATH",
							help="Directory where output files will be written")

		parser.add_argument("-n", dest="analysis_name", metavar="NAME",
							help="Analysis name")

		parser.add_argument("--output-format", dest="output_format", metavar="FORMAT",
							choices=["tsv", "tsv.gz", "tsv.bz2"], default="tsv",
							help="The FORMAT for the output file")

	def _check_args(self):
		if self.args.output_path is None:
			self.args.output_path = os.getcwd()

	def _error(self, msg):
		self._arg_errors += [msg]

	def load_mapping(self, matrix, path, filt=None):
		map = {}

		if path is None:
			return map

		with open(path) as f:
			for line in f:
				line = line.rstrip()
				if len(line) == 0 or line.startswith("#"):
					continue

				fields = line.split('\t')
				if len(fields) < 2:
					continue

				gene, pathway = fields[0:2]

				if filt is None or filt.valid(gene):
					if pathway in map:
						map[pathway].add(gene)
					else:
						map[pathway] = set([gene])

		for pathway, genes in map.items():
			map[pathway] = sorted(genes)

		return MatrixMapping(matrix, map)

	def save_splited_results(self, output_path, analysis_name, output_format,
								matrix, mapping, method, results, slices, suffix=""):

		if len(suffix) > 0:
			suffix = "-{0}".format(suffix)

		for slice_results_index, slice in enumerate(slices):
			slice_name = matrix.slice_names[slice]
			path = os.path.join(output_path, "{0}{1}-{2}.{3}".format(
				analysis_name, suffix, slice_name, output_format))
			self.log.debug("  > {0}".format(path))

			with tsv.open(path, 'w') as f:
				tsv.write_line(f, "## version={0}".format(VERSION))
				tsv.write_line(f, "## slice={0}".format(slice_name))
				tsv.write_line(f, "## method={0}".format(method.name))
				tsv.write_line(f, "## date={0}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
				for key, value in self.parameters:
					tsv.write_line(f, "## {0}={1}".format(key, value))
				tsv.write_line(f, "ID", *method.results_columns)
				for row_index, row_name in enumerate(mapping.group_names):
					value = results[slice_results_index, row_index]
					if not np.isnan(value):
						tsv.write_line(f, row_name, value, null_value="-")

	def save_matrix(self, output_path, analysis_name, output_format,
						 row_names, col_names, data,
						 suffix="", params=None, valid_row=lambda row: True):

		if len(suffix) > 0:
			suffix = "-{0}".format(suffix)

		if params is None:
			params = []

		path = os.path.join(output_path, "{0}{1}.{2}".format(analysis_name, suffix, output_format))
		self.log.debug("  > {0}".format(path))

		with tsv.open(path, 'w') as f:
			tsv.write_line(f, "## version={0}".format(VERSION))
			tsv.write_line(f, "## date={0}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
			for key, value in params + self.parameters:
				tsv.write_line(f, "## {0}={1}".format(key, value))
			tsv.write_line(f, "ID", *col_names)
			for row_index, row_name in enumerate(row_names):
				if len(row_name) == 0:
					self.log.warn("Empty identifier detected")
					continue

				row = data[row_index, :]
				if valid_row(row):
					values = [v if not np.isnan(v) else None for v in row]
					tsv.write_line(f, row_name, *values, null_value="-")

	#DEPRECATED
	def save_combined_results(self, output_path, analysis_name, output_format,
								method, row_names, col_names, data, suffix="combination"):

		self.log.info("Saving combination results ...")

		path = os.path.join(self.args.output_path,
							"{0}-{1}.{2}".format(self.args.analysis_name, suffix, self.args.output_format))

		self.log.debug("  > {0}".format(path))
		with tsv.open(path, 'w') as f:
			tsv.write_line(f, "## slices={0}".format(",".join(col_names)))
			tsv.write_line(f, "## method={0}".format(method.name))
			tsv.write_line(f, "## date={0}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
			for key, value in self.parameters:
				tsv.write_line(f, "## {0}={1}".format(key, value))
			tsv.write_line(f, "ID", *method.combination_columns)
			for row_index, row_name in enumerate(row_names):
				if not np.isnan(data[row_index, 0]):
					values = [v if not np.isnan(v) else None for v in data[row_index, :]]
					tsv.write_line(f, row_name, *values, null_value="-")

	def run(self):

		# Check arguments

		self._check_args()

		if len(self._arg_errors) > 0:
			for msg in self._arg_errors:
				self.log.error(msg)
			exit(-1)

		# Parse defined parameters

		self.parameters = []

		if self.args.defines is not None:
			for define in self.args.defines:
				try:
					pos = define.index("=")
					key = define[:pos]
					value = define[pos + 1:]

					self.parameters += [(key, value)]
				except:
					self.log.warn("Wrong parameter definition: {0}".format(define))
					continue
