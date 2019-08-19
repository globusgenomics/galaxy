import os
import shutil
import subprocess

from wok.task import task

from intogensm.projres import ProjectResults

@task.foreach()
def finish(project):
	log = task.logger
	conf = task.conf

	project_id = project["id"]

	log.info("Cleaning project {0} ...".format(project_id))

	projres = ProjectResults(project)

	projres.clean(conf)

	# save project.conf

	projres.save_def()

task.run()