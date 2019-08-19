from intogensm.command import RunCommand

class QcCommand(RunCommand):

	cmd_name = "qc"

	def execute(self):

		self._wok_run(self.instance_name, self.conf_builder, "intogen-mutations:qc")