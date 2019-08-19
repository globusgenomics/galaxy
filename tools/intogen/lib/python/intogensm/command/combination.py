from intogensm.command import RunCommand
		
class CombinationCommand(RunCommand):

	cmd_name = "combination"

	def execute(self):

		self._wok_run(self.instance_name, self.conf_builder, "intogen-mutations:combination")