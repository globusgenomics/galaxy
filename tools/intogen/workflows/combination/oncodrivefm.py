import os
import subprocess

from wok.task import task

from intogensm.woktmp import make_temp_dir, remove_temp
from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.utils import normalize_id
from intogensm.paths import get_combination_path

@task.foreach()
def combination_oncodrivefm(projects_set):
	log = task.logger
	conf = task.conf

	classifier, projects = projects_set

	classifier_id = classifier["id"]

	group_values = classifier["group_values"]
	short_values = classifier["group_short_values"]
	long_values = classifier["group_long_values"]

	group_name = classifier["group_name"]
	group_short_name = classifier["group_short_name"]
	group_long_name = classifier["group_long_name"]

	if len(group_values) == 0:
		group_file_prefix = classifier_id
	else:
		group_file_prefix = "{0}-{1}".format(classifier_id, group_short_name)

	group_file_prefix = normalize_id(group_file_prefix)

	log.info("--- [{0} ({1}) ({2}) ({3})] {4}".format(classifier["name"], group_long_name, group_short_name, group_name, "-" * 30))

	log.info("Exporting project data ...")

	base_path = make_temp_dir(task, suffix=".{0}".format(group_file_prefix))

	log.debug("> {0}".format(base_path))

	project_ids = []
	gene_files = []
	pathway_files = []
	for project in projects:
		project_id = project["id"]
		project_ids += [project_id]

		log.info("  Project {0}:".format(project["id"]))

		projdb = ProjectDb(project["db"])

		log.info("    Genes ...")

		count = 0
		file_path = os.path.join(base_path, "{0}-genes.tsv".format(project_id))
		gene_files += [file_path]
		with open(file_path, "w") as f:
			tsv.write_param(f, "classifier", classifier_id)
			tsv.write_param(f, "group_id", group_name)
			tsv.write_param(f, "slice", project_id)
			tsv.write_line(f, "GENE_ID", "PVALUE")
			for gene in projdb.genes():
				if gene.fm_pvalue is not None:
					tsv.write_line(f, gene.id, gene.fm_pvalue, null_value="-")
					count += 1

		log.info("      {0} genes".format(count))

		log.info("    Pathways ...")

		count = 0
		file_path = os.path.join(base_path, "{0}-pathways.tsv".format(project_id))
		pathway_files += [file_path]
		with open(file_path, "w") as f:
			tsv.write_param(f, "classifier", classifier_id)
			tsv.write_param(f, "group_id", group_name)
			tsv.write_param(f, "slice", project_id)
			tsv.write_line(f, "PATHWAY_ID", "ZSCORE")
			for pathway in projdb.pathways():
				if pathway.fm_zscore is not None:
					tsv.write_line(f, pathway.id, pathway.fm_zscore, null_value="-")
					count += 1

		log.info("      {0} pathways".format(count))

		projdb.close()

	log.info("Combining ...")

	combination_path = get_combination_path(conf, "oncodrivefm")

	log.info("  Genes ...")

	cmd = " ".join([
			"oncodrivefm-combine",
			"-m median-empirical",
			"-o '{0}'".format(combination_path),
			"-n 'gene-{0}'".format(group_file_prefix),
			"-D 'classifier={0}'".format(classifier_id),
			"-D 'group_id={0}'".format(group_name),
			"-D 'group_short_name={0}'".format(group_short_name),
			"-D 'group_long_name={0}'".format(group_long_name),
			"--output-format tsv.gz"
	] + ["'{0}'".format(name) for name in gene_files])

	log.debug(cmd)

	ret_code = subprocess.call(cmd, shell=True)
	if ret_code != 0:
		#log.error("OncodriveFM error while combining gene pvalues:\n{0}".format(cmd))
		#return -1
		raise Exception("OncodriveFM error while combining gene pvalues:\n{0}".format(cmd))

	log.info("  Pathways ...")

	cmd = " ".join([
			"oncodrivefm-combine",
			"-m median-zscore",
			"-o '{0}'".format(combination_path),
			"-n 'pathway-{0}'".format(group_file_prefix),
			"-D 'classifier={0}'".format(classifier_id),
			"-D 'group_id={0}'".format(group_name),
			"-D 'group_short_name={0}'".format(group_short_name),
			"-D 'group_long_name={0}'".format(group_long_name),
			"--output-format tsv.gz"
	] + ["'{0}'".format(name) for name in pathway_files])

	log.debug(cmd)

	ret_code = subprocess.call(cmd, shell=True)
	if ret_code != 0:
		#log.error("OncodriveFM error while combining pathway zscores:\n{0}".format(cmd))
		#return -1
		raise Exception("OncodriveFM error while combining pathway zscores:\n{0}".format(cmd))

	remove_temp(task, base_path)

task.run()