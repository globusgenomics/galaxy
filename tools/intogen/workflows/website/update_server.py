import os
import re
import fcntl

from wok.task import task

from intogensm.onexus import Onexus

@task.foreach()
def update_server(project):
	log = task.logger
	conf = task.conf

	project_id = project["id"]

	log.info("--- [{0}] --------------------------------------------".format(project_id))

	projects_list = conf.get("website.projects_list")
	if projects_list is None:
		log.warn("The project can not be registered on the website server as the configuration 'website.projects_list' is not defined")
		return

	user_id = conf.get("website.user_id")
	if user_id is None:
		log.warn("The project can not be registered on the website server as the 'user_id' was not found in the configuration")
		return

	website_path = project["website"]

	onexus = Onexus(projects_list, logger=log)
	onexus.add_project(user_id, project_id, website_path)

task.run()