import os
import os.path

from intogensm.command import RunCommand
		
class AnalysisGenesCommand(RunCommand):

	cmd_name = "analysis-genes"

	def execute(self):

		self._wok_run(self.instance_name, self.conf_builder, "intogen-mutations:genes")
