import os
import os.path

from intogensm.command import RunCommand
		
class AnalysisResultsCommand(RunCommand):

	cmd_name = "results"

	def execute(self):

		self._wok_run(self.instance_name, self.conf_builder, "intogen-mutations:results")
