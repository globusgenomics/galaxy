import os
import os.path
import shutil

from intogensm.command import RunCommand
from intogensm.utils import normalize_id
from intogensm import common

PROJECT_EPILOG="""
Configurable parameters:

  You will find all the available configuration parameters from

  http://www.intogen.org/sm/help/configuration.html
"""

class AnalysisCommand(RunCommand):

	cmd_name = "analysis"

	def __init__(self):
		RunCommand.__init__(self,
			args_usage = "<variants-file> [<variants-file> ...]",
			epilog = PROJECT_EPILOG)
		
		self.variants_files = []
		
		# Gather variants files from arguments

		for var_file in self.args.files:
			if not os.path.isabs(var_file):
				var_file = os.path.join(os.getcwd(), var_file)

			if not os.path.exists(var_file):
				self.log.error("Variants file not found: {}".format(var_file))
				exit(-1)

			if not os.path.isfile(var_file):
				self.log.error("A file is expected: {}".format(var_file))
				exit(-1)

			self.variants_files += [var_file]

		# Get project id

		if self.args.project_id is None:
			self.log.error("Project identifier not specified.")
			exit(-1)

		self.project_id = normalize_id(self.args.project_id)
		
	def add_arguments(self, parser):
		RunCommand.add_arguments(self, parser)

		parser.add_argument("files", nargs="+",
			help="Variants files")

		parser.add_argument("-p", "--id", dest="project_id", metavar="NAME",
			help="Define the project identifier. Required.")
		
		parser.add_argument("-a", "--assembly", dest="assembly", metavar="ASSEMBLY",
			choices=["hg18", "hg19"], default="hg19",
			help="Define the assembly [hg18, hg19]. Default is hg19.")

		parser.add_argument("--single-tumor", dest="single_tumor", action="store_true", default=False,
			help="Run a single tumor analysis instead of the regular cohort analysis.")

	def build_conf(self):
		if self.instance_name is None:
			self.instance_name = self.project_id

		RunCommand.build_conf(self)
		
		self.conf_builder.add_value("project.id", self.project_id)

		if self.args.single_tumor:
			self.conf_builder.add_value("variants_only", True)
			self.conf_builder.add_value("skip_oncodrivefm", True)
			self.conf_builder.add_value("skip_oncodriveclust", True)

	def execute(self):
		# Copy variants files into temp_path
		
		dest_path = os.path.join(self.results_path, self.workspace, "projects", self.project_id)
		if not os.path.exists(dest_path):
			os.makedirs(dest_path)

		var_files = []
		for i, var_file in enumerate(self.variants_files):
			file_name = "{0:02}-{1}".format(i, os.path.basename(var_file))
			dest_mut_file = os.path.join(dest_path, file_name)
			shutil.copy2(var_file, dest_mut_file)
			var_files += [dest_mut_file]

		self.variants_files = var_files

		self.conf_builder, self.flow = common.project_analysis(
					self.variants_files,
					assembly = self.args.assembly,
					conf_builder = self.conf_builder)

		self._wok_run(self.instance_name, self.conf_builder, self.flow)
