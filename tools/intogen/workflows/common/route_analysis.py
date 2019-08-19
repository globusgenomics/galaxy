import os

from wok.task import task

@task.begin()
def begin():
	log = task.logger
	conf = task.conf

	task.context["skip_recurrences"] = conf.get("skip_recurrences", False, dtype=bool)
	task.context["skip_oncodrivefm"] = conf.get("skip_oncodrivefm", False, dtype=bool)
	task.context["skip_oncodriveclust"] = conf.get("skip_oncodriveclust", False, dtype=bool)

	log.info("Routing projects into analysis ...")
	log.info("  skip_recurrences = {0}".format(task.context["skip_recurrences"]))
	log.info("  skip_oncodrivefm = {0}".format(task.context["skip_oncodrivefm"]))
	log.info("  skip_oncodriveclust = {0}".format(task.context["skip_oncodriveclust"]))
	log.info("")

@task.foreach()
def route(project):
	log = task.logger
	conf = task.conf

	projects_port = task.ports("projects_out")
	recurrences_projects_port = task.ports("recurrences_projects")
	oncodrivefm_projects_port = task.ports("oncodrivefm_projects")
	oncodriveclust_projects_port = task.ports("oncodriveclust_projects")

	log.info("  {0}".format(project["id"]))

	projects_port.send(project)

	if not task.context["skip_recurrences"]:
		recurrences_projects_port.send(project)

	if not task.context["skip_oncodrivefm"]:
		oncodrivefm_projects_port.send(project)

	if not task.context["skip_oncodriveclust"]:
		oncodriveclust_projects_port.send(project)

task.run()