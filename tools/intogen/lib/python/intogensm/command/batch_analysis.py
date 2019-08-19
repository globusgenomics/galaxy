import os
import os.path
import re

from intogensm.command import RunCommand
from intogensm import common

from intogensm.command.analysis import PROJECT_EPILOG
		
class BatchAnalysisCommand(RunCommand):

	cmd_name = "batch-analysis"

	def __init__(self):
		RunCommand.__init__(self,
			args_usage = "<projects-path> [<projects-path> ...]",
			epilog = PROJECT_EPILOG)
		
		self.scan_paths = []
		
		# Gather scan paths from arguments

		for scan_path in self.args.paths:
			if not os.path.isabs(scan_path):
				scan_path = os.path.join(os.getcwd(), scan_path)

			if not os.path.exists(scan_path):
				self.log.error("Path not found: {}".format(scan_path))
				exit(-1)

			self.scan_paths += [scan_path]

	def add_arguments(self, parser):
		RunCommand.add_arguments(self, parser)

		parser.add_argument("paths", nargs="+",
							help="Projects paths")

		parser.add_argument("-i", "--include", dest = "include", metavar = "ID", action="append",
			help = "Include project which id=ID")

		parser.add_argument("-I", "--include-regex", dest = "include_regex", metavar = "REGEX", action="append",
			help = "Include project which id matches REGEX")

		parser.add_argument("--include-from", dest = "include_from", metavar = "FILE", action="append",
			help = "Include project which id found in file FILE")

		parser.add_argument("-e", "--exclude", dest = "exclude", metavar = "ID", action="append",
			help = "Exclude project which id=ID")

		parser.add_argument("-E", "--exclude-regex", dest = "exclude_regex", metavar = "REGEX", action="append",
			help = "Exclude project which id matches REGEX")

		parser.add_argument("--exclude-from", dest = "exclude_from", metavar = "FILE", action="append",
			help = "Exclude project which id found in file FILE")
		
	def build_conf(self):
		RunCommand.build_conf(self)

	def execute(self):

		# Gather includes and excludes from options

		includes = []

		if self.args.include is not None:
			for inc in self.args.include:
				includes += ["^{0}$".format(re.escape(inc))]
		if self.args.include_regex is not None:
			for inc_regex in self.args.include_regex:
				includes += [inc_regex]
		if self.args.include_from is not None:
			for file in self.args.include_from:
				with open(file, "r") as f:
					for line in f:
						line = line.strip()
						if line.startswith("#") or len(line) == 0:
							continue
						includes += ["^{0}$".format(re.escape(line))]

		excludes = []

		if self.args.exclude is not None:
			for exc in self.args.exclude:
				excludes += ["^{0}$".format(re.escape(exc))]
		if self.args.exclude_regex is not None:
			for exc_regex in self.args.exclude_regex:
				excludes += [exc_regex]
		if self.args.exclude_from is not None:
			for file in self.args.exclude_from:
				with open(file, "r") as f:
					for line in f:
						line = line.strip()
						if line.startswith("#") or len(line) == 0:
							continue
						excludes += ["^{0}$".format(re.escape(line))]

		# Create the wok engine and the workflow instance
		
		self.conf_builder, self.flow = common.projects_analysis(
					self.scan_paths,
					includes=includes,
					excludes=excludes,
					conf_builder=self.conf_builder)
		
		self._wok_run(self.instance_name, self.conf_builder, self.flow)

