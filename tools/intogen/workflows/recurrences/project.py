import os
import json

from wok.task import task

from intogensm.projdb import ProjectDb

@task.foreach()
def projects(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	log.info("--- [{0}] --------------------------------------------".format(project["id"]))

	projdb = ProjectDb(project["db"])

	total_samples = projdb.get_total_affected_samples()

	if total_samples == 0:
		log.warn("There are no samples, recurrences cannot be calculated.")
		projdb.close()
		return

	log.info("Calculating project recurrences for variant genes ...")

	projdb.compute_affected_genes_recurrences(total_samples)

	if not conf.get("variants_only", False):

		log.info("Calculating project recurrences for genes ...")

		projdb.compute_gene_recurrences(total_samples)

		log.info("Calculating project recurrences for pathways ...")

		projdb.compute_pathway_recurrences(total_samples)

	projdb.commit()
	projdb.close()

	projects_out_port.send(project)

task.run()