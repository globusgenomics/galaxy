#! /usr/bin/env python

import sys

from intogensm.command.analysis import AnalysisCommand
from intogensm.command.batch_analysis import BatchAnalysisCommand
from intogensm.command.analysis_genes import AnalysisGenesCommand
from intogensm.command.results import AnalysisResultsCommand
from intogensm.command.combination import CombinationCommand
from intogensm.command.summary import SummaryCommand
from intogensm.command.qc import QcCommand
from intogensm.command.qc_report import QcReportCommand

_COMMANDS = [
	AnalysisCommand,
	BatchAnalysisCommand,
	AnalysisGenesCommand,
	AnalysisResultsCommand,
	CombinationCommand,
	SummaryCommand,
	QcCommand,
	QcReportCommand
]

_COMMANDS_MAP = {}
for cmd_class in _COMMANDS:
	_COMMANDS_MAP[cmd_class.cmd_name] = cmd_class

if __name__ == "__main__":
	supported_commands = [cmd.cmd_name for cmd in _COMMANDS]

	if len(sys.argv) < 2:
		print "usage: run <command> [<args>]"
		print
		print "Supported commands: {0}".format(", ".join(supported_commands))
		print
		exit(-1)

	cmd = sys.argv[1]
	
	del sys.argv[1]
	
	if cmd not in _COMMANDS_MAP:
		print "Unknown command: {0}".format(cmd)
		print "Supported commands: {0}".format(", ".join(supported_commands))
		exit(-1)

	command = _COMMANDS_MAP[cmd]()

	exit(command.run())
