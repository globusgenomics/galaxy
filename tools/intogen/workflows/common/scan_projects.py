import os
import subprocess

from wok.task import task

from intogensm.utils import list_projects
from intogensm.paths import get_projects_path, get_temp_path

@task.source()
def scan_projects(project_out):
	log = task.logger
	conf = task.conf

	projects_path = get_projects_path(conf)

	log.info("Scanning projects ...")

	count = 0

	sent_projects = []

	for path, project in list_projects(log, projects_path):
		project_path = os.path.dirname(path)

		if "id" not in project:
			log.warn("Discarding project that doesn't have 'id': {0}".format(path))
			continue

		project_id = project["id"]

		if "name" in project:
			log.info("--- [{0}: {1}] ---------------------------------".format(project_id, project["name"]))
		else:
			log.info("--- [{0}] ---------------------------------".format(project_id))

		if "db" not in project:
			project["db"] = os.path.join(project_path, "project.db.gz")
		elif not os.path.isabs(project["db"]):
			project["db"] = os.path.join(project_path, project["db"])

		if not os.path.exists(project["db"]):
			log.error("Project database not found at {0}".format(os.path.relpath(project["db"], project_path)))
			continue

		if project["db"].endswith(".gz"):
			log.info("Uncompressing project database ...")
			retcode = subprocess.call("gunzip -fc {0} >{1}".format(project["db"], project["db"][:-3]), shell=True)
			if retcode != 0:
				log.error("Unexpected error while uncompressing the project database at {0}".format(os.path.relpath(project["db"], projects_path)))
				continue

			project["db"] = project["db"][:-3]

		temp_path = get_temp_path(conf, project_id)
		if not os.path.exists(temp_path):
			os.makedirs(temp_path)

		project["path"] = project_path
		project["temp_path"] = temp_path

		sent_projects += [project_id]
		project_out.send(project)

		count += 1

	log.info("Found {0} projects:\n  - {1}".format(count, "\n  - ".join(sent_projects)))

task.run()