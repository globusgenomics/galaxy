import os
import subprocess

from wok.task import task

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.model import Gene, Pathway

@task.foreach()
def update_db(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	project_id = project["id"]
	log.info("--- [{0}] --------------------------------------------".format(project_id))

	oclust = project["oncodriveclust"]
	del project["oncodriveclust"]

	if not os.path.exists(oclust["results"]):
		log.warn("No results have been found. Skipping it.")
		return

	log.info("Updating the project database ...")

	projdb = ProjectDb(project["db"])

	exc_path = os.path.join(project["temp_path"], "oncodriveclust-excluded-cause.tsv")

	log.info("  Excluded gene causes ...")
	log.debug("    > {0}".format(exc_path))

	count = 0
	with tsv.open(exc_path, "r") as exf:
		for gene, cause in tsv.lines(exf, (str, str), header=True):
			projdb.update_gene(Gene(id=gene, clust_exc_cause=cause))
			count += 1

	log.debug("    {0} genes excluded".format(count))

	log.info("  OncodriveCLUST results ...")

	with tsv.open(oclust["results"], "r") as f:
		types = (str, str, float, float, float)
		columns = ("GENE", "CLUST_COORDS", "ZSCORE", "PVALUE", "QVALUE")
		for gene, coords, zscore, pvalue, qvalue in tsv.lines(f, types, columns=columns, header=True, null_value="NA"):
			projdb.update_gene(Gene(id=gene, clust_coords=coords, clust_zscore=zscore, clust_pvalue=pvalue,
									clust_qvalue=qvalue, clust_exc_cause=ProjectDb.NO_GENE_EXC))

	projdb.commit()

	projdb.close()

	projects_out_port.send(project)

task.run()