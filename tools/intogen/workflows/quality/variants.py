import os

from wok.task import task

from intogensm.projdb import ProjectDb
from intogensm.projres import ProjectResults

@task.foreach()
def variants(project):
	log = task.logger
	conf = task.conf

	log.info("--- [{0}] --------------------------------------------".format(project["id"]))

	log.info("Calculating number of variants processed in each step ...")

	proj_res = ProjectResults(project)

	projdb = ProjectDb(project["db"])

	counts = projdb.count_variants()

	proj_res.save_quality_control("variants", counts)

	projdb.close()

task.run()