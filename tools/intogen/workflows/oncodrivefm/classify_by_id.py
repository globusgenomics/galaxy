import os
import subprocess

from wok.task import task

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.model import Gene, Pathway

@task.begin()
def begin():
	log = task.logger
	conf = task.conf

	log.info("Classifying results by project id ...")

	task.context["combinations"] = {}

@task.foreach()
def classify(project):
	log = task.logger
	conf = task.conf

	combinations = task.context["combinations"]

	project_id = project["id"]

	ofm = project["oncodrivefm"]
	feature = ofm["feature"]

	log.info("Project {0}, feature {1}".format(project_id, feature))

	key = project_id

	if key not in combinations:
		project["oncodrivefm"] = [(feature, ofm["results"])]
		combinations[key] = project
	else:
		project = combinations[key]
		project["oncodrivefm"] += [(feature, ofm["results"])]

@task.end()
def end():
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	log.info("Sending projects ...")

	combinations = task.context["combinations"]

	for key, project in sorted(combinations.items(), key=lambda v: v[0]):
		log.debug(key)
		projects_out_port.send(project)

task.run()