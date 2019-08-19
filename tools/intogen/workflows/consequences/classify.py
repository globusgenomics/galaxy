import os.path

from intogensm.projdb import ProjectDb

from wok.task import task

@task.begin()
def begin():
	task.context = {}

@task.foreach()
def join_partitions(partition):

	project = partition["project"]

	project_id = project["id"]

	project_parts = task.context

	if project_id not in project_parts:
		project_parts[project_id] = []

	project_parts[project_id] += [partition]

@task.end()
def end():
	log = task.logger

	projects_port = task.ports("projects")

	for project_id, partitions in task.context.items():
		log.info("Project {0}: {1} partitions".format(project_id, len(partitions)))

		project = partitions[0]["project"]
		project["partitions"] = part_list = list()
		for part in partitions:
			del part["project"]
			part_list += [part]

		projects_port.send(project)

task.run()