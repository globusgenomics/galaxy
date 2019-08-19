import sys
import os
import os.path
import subprocess
import json
import getpass
import signal

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from datetime import datetime

from wok import logger
from wok.config.data import Data
from wok.config.builder import ConfigBuilder
from wok.engine import WokEngine

from intogensm import VERSION
from intogensm.utils import normalize_id

def keyboardinterrupt_handler(*args):
	raise KeyboardInterrupt

class Command(object):
	def __init__(self, args_usage="", epilog="", logger_name=None):

		self.args_usage = args_usage
		self.epilog = epilog

		if logger_name is None:
			if hasattr(self, "cmd_name"):
				logger_name = self.cmd_name
			else:
				logger_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]

		# Determine binary path

		script_path = os.path.dirname(sys.argv[0])
		if script_path:
			self.launchers_path = os.path.abspath(script_path)
		else:
			self.launchers_path = os.getcwd()

		# Determine default paths

		self.root_path = os.path.normpath(os.path.join(self.launchers_path, ".."))

		self.bin_path = os.path.join(self.root_path, "bin")

		self.data_path = os.path.join(self.root_path, "data")

		self.conf_path = os.path.join(self.root_path, "conf")

		self.workflows_path = os.path.join(self.root_path, "workflows")

		self.runtime_path = os.path.join(self.root_path, "runtime")

		# Parse arguments

		parser = ArgumentParser(
			prog = "run " + self.cmd_name,
			epilog = self.epilog,
			formatter_class=RawDescriptionHelpFormatter)

		parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

		self.add_arguments(parser)

		parser.add_argument("-L", "--log-level", dest="log_level",
						  default=None, choices=["debug", "info", "warn", "error", "critical", "notset"],
						  help="Define the logging level")

		self.args = parser.parse_args()

		# Create preliminary logger

		level = self.args.log_level or "info"
		logger.initialize(level=level)


		self.log = logger.get_logger(logger_name, level=level)

	def add_arguments(self, parser):
		"Called from __init__() to let the command to add more options to the OptionsParser"
		pass

	def execute(self):
		"Called from run() when the configuration is ready"
		raise Exception("Abstract method called")

	def run(self):
		"Main entry point"

		self.execute()

class RunCommand(Command):

	DEFAULT_CONF_FILES = [
		"system.conf"
	]

	DEFAULT_REQUIRED_CONF = [
#		"work_path", "temp_path"]
	]

	def __init__(self, args_usage="", epilog="", flow_file=None, conf_files=None, conf_keys=None, logger_name=None):

		Command.__init__(self, args_usage, epilog, logger_name)

		signal.signal(signal.SIGINT, keyboardinterrupt_handler)
		signal.signal(signal.SIGTERM, keyboardinterrupt_handler)

		if conf_files is None:
			conf_files = []
		if conf_keys is None:
			conf_keys = []

		self.flow_file = flow_file
		self.conf_files = self.DEFAULT_CONF_FILES + conf_files
		self.conf_keys = self.DEFAULT_REQUIRED_CONF + conf_keys

		# Workspace

		self.workspace = self.args.workspace

		# Instance name

		self.instance_name = self.args.instance_name

		if self.instance_name is not None:
			self.instance_name = normalize_id(self.instance_name)

		'''
		# Override configuration path if required

		if self.args.conf_path is not None:
			self.conf_path = os.path.abspath(self.args.conf_path)

		# Get required configuration files and override system.conf if required

		if self.args.system_conf is not None:
			req_conf_files = []
			for cf in self.conf_files:
				if cf == "system.conf":
					req_conf_files += [self.args.system_conf]
				else:
					req_conf_files += [cf]
		else:
			req_conf_files = self.conf_files
		'''

		req_conf_files = self.conf_files
		
		# Determine required and user configuration files and data

		self.required_conf_files = [os.path.join(self.conf_path, cf)
										for cf in req_conf_files]

		if self.args.conf_files is not None:
			self.user_conf_files = []
			for cf in self.args.conf_files:
				if not os.path.isabs(cf):
					cf = os.path.join(os.getcwd(), cf)
				self.user_conf_files += [cf]
		else:
			self.user_conf_files = []

		if self.args.conf_data is not None:
			self.user_conf_data = self.args.conf_data
		else:
			self.user_conf_data = []

		# max cores

		if self.args.max_cores is None:
			self.max_cores = 0
		else:
			self.max_cores = self.args.max_cores

		# Prepare extra configuration data

		self.extra_conf_data = []

	def add_arguments(self, parser):
		"Called from __init__() to let the command to add more options to the OptionsParser"

		parser.add_argument("-w", "--workspace", dest = "workspace", metavar = "NAME", default="default",
						  help = "Define the workspace name.")

		parser.add_argument("-u", "--user", dest="user", metavar="NAME",
						  help="Define the user name.")

		g = parser.add_argument_group("Wok Options")

		g.add_argument("-n", "--instance", dest = "instance_name", metavar = "NAME",
						  help = "Define the instance name")

		g.add_argument("-C", "--conf", action="append", dest="conf_files", metavar="FILE",
					 help="Load configuration from a file. Multiple files can be specified")

		g.add_argument("-D", action="append", dest="conf_data", metavar="PARAM=VALUE",
					 help="Define configuration parameter. Multiple definitions can be specified. Example: -D option1=value")

		g.add_argument("-j", dest="max_cores", metavar="NUM", type=int,
					 help="Define the maximum number of cores to use. Default all the cores available.")

	def build_conf(self):
		"Called from run() to prepare configuration before expansion of values"

		conf_files = self.required_conf_files

		user_conf_file = os.path.join(self.conf_path, "user.conf")
		if os.path.exists(user_conf_file):
			conf_files += [user_conf_file]

		conf_files += self.user_conf_files
		
		# Check that configuration files exist

		missing_conf_files = [cf for cf in conf_files if not os.path.exists(cf)]
		if len(missing_conf_files) > 0:
			self.log.error("Configuration files not found:\n{}".format(
							"\n".join("  {}".format(cf) for cf in missing_conf_files)))
			exit(-1)
		
		# Build the configuration
			
		self.conf_builder = ConfigBuilder()

		if self.instance_name is None:
			if self.workspace is not None:
				self.instance_name = self.workspace
			else:
				self.instance_name = datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f")

		self.conf_builder.add_value("wok.instance.name", self.instance_name)

		# Files
		
		for cf in conf_files:
			self.conf_builder.add_file(cf)

		# General conf

		self.conf_builder.add_value("workspace", self.workspace)

		self.conf_builder.add_value("workflows_path", self.workflows_path)

		if self.max_cores > 0:
			self.conf_builder.add_value("wok.platform.jobs.max_cores", self.max_cores)

		if self.args.log_level is not None:
			self.conf_builder.add_value("log.level", self.args.log_level)
			self.conf_builder.add_value("wok.engine", self.args.log_level)

		# Incorporate user and extra configuration data

		for data in self.extra_conf_data + self.user_conf_data:
			try:
				pos = data.index("=")
				key = data[0:pos]
				value = data[pos+1:]
				try:
					v = json.loads(value)
				except:
					v = value
				self.conf_builder.add_value(key, Data.create(v))
			except:
				raise Exception("Wrong configuration data: KEY=VALUE expected but found '{}'".format(data))

	def process_conf(self):
		"Called from run() to process the configuration after expansion of values"

		# Validate that required keys exist

		mk = self.conf.missing_fields(self.conf_keys)
		if len(mk) > 0:
			sb = ["The following configuration parameters were not found:\n",
					"\n".join("* " + k for k in mk),
					"\nin any of the following configuration files:\n",
					"\n".join("* " + k for k in conf_files)]
			self.log.error("".join(sb))
			exit(-1)

		# Configurable paths

		if "runtime_path" in self.conf:
			self.runtime_path = self.conf["runtime_path"]
		else:
			self.conf_builder.add_value("runtime_path", self.runtime_path)

		# Configure wok work path

		if "wok.work_path" not in self.conf:
			self.conf_builder.add_value("wok.work_path", os.path.join(self.runtime_path, "wok"))

		# Prepare paths by user_id

		user_id = self.args.user or getpass.getuser()

		rpath = self.conf.get("results_path", os.path.join(self.runtime_path, "results"))
		self.results_path = os.path.join(rpath, user_id)
		self.conf_builder.add_value("results_path", self.results_path)
		if not os.path.exists(self.results_path):
			os.makedirs(self.results_path)

		tpath = self.conf.get("temp_path", os.path.join(self.runtime_path, "temp"))
		self.temp_path = os.path.join(tpath, user_id)
		self.conf_builder.add_value("temp_path", self.temp_path)
		if not os.path.exists(self.temp_path):
			os.makedirs(self.temp_path)
		
	def run(self):
		"Run the command and execute the command"

		self.build_conf()

		# Expand configuration variables

		self.conf = self.conf_builder.get_conf()
		
		self.conf.expand_vars()
		#self.log.debug(repr(self.conf))

		# Validate and process configuration

		self.process_conf()
		
		# Regenerate configuration
		
		self.conf = self.conf_builder.get_conf()

		# Final logging configuration

		log = logger.get_logger("")
		log.removeHandler(log.handlers[0])
		logging_conf = self.conf.get("wok.logging")
		logger.initialize(logging_conf)

		# Show some debugging information

		self.log.debug("Root path = {}".format(self.root_path))
		self.log.debug("Conf path = {}".format(self.conf_path))
		self.log.debug("Data path = {}".format(self.data_path))
		self.log.debug("Workflows path = {}".format(self.workflows_path))
		self.log.debug("Runtime path = {}".format(self.runtime_path))
		self.log.debug("Results path = {}".format(self.results_path))
		self.log.debug("Temp path = {}".format(self.temp_path))

		if len(self.user_conf_files) > 0:
			self.log.debug("User defined configuration files:\n{}".format(
							"\n".join("  {}".format(cf) for cf in self.user_conf_files)))

		self.log.debug("Effective configuration: " + str(self.conf))
		
		# Execute

		try:
			self.execute()
		except Exception as ex:
			self.log.exception(ex)
			return -1

		return 0

	def _wok_run(self, case_name, conf_builder, flow_uri):

		wok = WokEngine(self.conf)

		wok.start(wait=False, single_run=True)

		try:
			case = wok.create_case(case_name, conf_builder, flow_uri)
			wok.wait()
		except KeyboardInterrupt:
			pass
		finally:
			try:
				wok.stop()
			except KeyboardInterrupt:
				self.log.warn("Ctrl-C pressed while stopping Wok engine.")