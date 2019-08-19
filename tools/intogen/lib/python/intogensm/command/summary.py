from intogensm.command import RunCommand

class SummaryCommand(RunCommand):

	cmd_name = "summary"

	def execute(self):

		self._wok_run(self.instance_name, self.conf_builder, "intogen-mutations:summary")