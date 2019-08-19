import os
import subprocess

from wok.task import task

from bgcore import tsv

from intogensm.projdb import ProjectDb

@task.foreach()
def oncoclust(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	project_id = project["id"]
	log.info("--- [{0}] --------------------------------------------".format(project_id))

	gene_transcripts_path = os.path.join(conf["data_path"], "Ensembl", "gene_transcripts.tsv")

	oclust = project["oncodriveclust"]

	data_paths = oclust["data_paths"]

	mutations_threshold = oclust["mutations_threshold"]

	oclust["results"] = os.path.join(project["temp_path"], "oncodriveclust-results.tsv")

	cmd = " ".join([
		"oncodriveclust",
		"-c",
		"-m", str(mutations_threshold),
		"-o", oclust["results"],
		data_paths[0],
		data_paths[1],
		gene_transcripts_path
	])

	log.debug(cmd)

	ret_code = subprocess.call(cmd, shell=True)
	if ret_code == 1:
		log.warn("No results were generated")
	elif ret_code != 0:
		log.error("Error while executing OncodriveCLUST:\n{0}".format(cmd))

	projects_out_port.send(project)

task.run()