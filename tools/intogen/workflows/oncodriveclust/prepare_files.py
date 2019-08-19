import os
import json

from wok.task import task

from bgcore import tsv
from bgcore.labelfilter import LabelFilter

from intogensm.projdb import ProjectDb
from intogensm.projres import ProjectResults
from intogensm.utils import get_project_conf
from intogensm.constants.oncodriveclust import *
from intogensm.paths import get_data_gene_filter_path
from intogensm.oncodriveclust import get_oncodriveclust_configuration, load_cds_len, retrieve_data, NON_SYN, SYN

@task.foreach()
def prepare_files(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	project_id = project["id"]
	log.info("--- [{0}] --------------------------------------------".format(project_id))

	project_results = ProjectResults(project)

	mutations_threshold, genes_filter_enabled, genes_filter, filt = get_oncodriveclust_configuration(log, conf, project)

	log.info("Loading transcripts CDS length ...")

	cds_len = load_cds_len(conf)

	log.info("Retrieving gene alterations ...")

	projdb = ProjectDb(project["db"])

	data = retrieve_data(projdb, cds_len)

	projdb.close()

	data_paths = [
		os.path.join(project["temp_path"], "oncodriveclust-non-syn-data.tsv"),
		os.path.join(project["temp_path"], "oncodriveclust-syn-data.tsv")]

	log.info("Saving data ...")
	log.debug("> {0}".format(data_paths[NON_SYN]))
	log.debug("> {0}".format(data_paths[SYN]))

	df = [tsv.open(path, "w") for path in data_paths]

	gene_sample_count = {}

	for key, value in data.items():
		findex, gene, sample = key
		transcript, transcript_len, protein_pos = value

		if findex == NON_SYN:
			if gene not in gene_sample_count:
				gene_sample_count[gene] = 1
			else:
				gene_sample_count[gene] += 1

			if genes_filter_enabled and not filt.valid(gene):
				continue

		tsv.write_line(df[findex], gene, sample, protein_pos)

	for f in df:
		f.close()

	exc_path = os.path.join(project["temp_path"], "oncodriveclust-excluded-cause.tsv")

	log.info("Saving excluded gene causes ...")
	log.debug("> {0}".format(exc_path))

	with tsv.open(exc_path, "w") as exf:
		tsv.write_line(exf, "GENE", "EXCLUDED_CAUSE")
		for gene, sample_count in gene_sample_count.items():
			causes = []
			if genes_filter_enabled and not filt.valid(gene):
				causes += [ProjectDb.GENE_EXC_FILTER]
			if sample_count < mutations_threshold:
				causes += [ProjectDb.GENE_EXC_THRESHOLD]
			if len(causes) > 0:
				tsv.write_line(exf, gene, "".join(causes))

	log.info("Sending project ...")

	projects_out_port.send(dict(project,
								oncodriveclust=dict(
									data_paths=data_paths,
									mutations_threshold=mutations_threshold,
									genes_filter_enabled=genes_filter_enabled, # not used
									genes_filter=genes_filter))) # not used

task.run()