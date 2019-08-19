import os
import subprocess

from wok.task import task

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.model import Gene, Pathway

@task.begin()
def begin():
	task.logger.info("Classifying slices by project ...")

@task.foreach()
def update_db(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	project_id = project["id"]
	ofm = project["oncodrivefm"]

	if project_id in task.context:
		task.context[project_id] += [project]
	else:
		task.context[project_id] = [project]

@task.end()
def end():
	log = task.logger

	projects_out_port = task.ports("projects_out")

	log.info("Updating the projects database ...")

	for project_id, projects in task.context.items():

		log.info("[{0}]".format(project_id))

		for index, project in enumerate(projects):
			projdb = ProjectDb(project["db"])

			if index == 0:
				log.info("  Functional impact ...")

				projdb.delete_sample_gene_fimpact()

				with tsv.open(project["sample_gene_fi_data"], "r") as f:
					types = (int, str, float, float, int, float, float, int, float, float, int)
					for fields in tsv.lines(f, types, header=True, null_value="-"):
						projdb.add_sample_gene_fimpact(*fields)

			ofm = project["oncodrivefm"]
			del project["oncodrivefm"]

			exc_path = os.path.join(project["temp_path"], "oncodrivefm-excluded-cause.tsv")

			log.info("  Excluded gene causes ...")
			log.debug("    > {0}".format(exc_path))

			count = 0
			with tsv.open(exc_path, "r") as exf:
				for gene, cause in tsv.lines(exf, (str, str), header=True):
					projdb.update_gene(Gene(id=gene, fm_exc_cause=cause))
					count += 1

			log.debug("    {0} genes excluded".format(count))

			for feature, results_path in ofm:

				log.info("  {0} ...".format(feature))
				log.debug("    > {0}".format(results_path))

				if feature == "genes":
					with tsv.open(results_path, "r") as f:
						count = 0
						for gene, pvalue, qvalue in tsv.lines(f, (str, float, float), header=True):
							projdb.update_gene(Gene(id=gene, fm_pvalue=pvalue,
													fm_qvalue=qvalue, fm_exc_cause=ProjectDb.NO_GENE_EXC))
							count += 1
						log.info("    {0} genes".format(count))
				elif feature == "pathways":
					with tsv.open(results_path, "r") as f:
						count = 0
						for pathway, zscore, pvalue, qvalue in tsv.lines(f, (str, float, float, float), header=True):
							projdb.update_pathway(Pathway(id=pathway, fm_zscore=zscore,
														  fm_pvalue=pvalue, fm_qvalue=qvalue))
							count += 1
						log.info("    {0} pathways".format(count))

			projdb.commit()

			projdb.close()

		projects_out_port.send(projects[0])

task.run()