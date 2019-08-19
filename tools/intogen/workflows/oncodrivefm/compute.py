import os
import subprocess

from wok.task import task
from wok.config.data import Data

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.model import Gene, Pathway
from intogensm.paths import get_data_kegg_path

@task.foreach()
def compute(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	project_id = project["id"]
	log.info("--- [{0}] --------------------------------------------".format(project_id))

	ofm = Data.element(project["oncodrivefm"])

	feature = ofm["feature"]
	slice_name = ofm["slice"]
	estimator = ofm.get("estimator")
	num_cores = ofm.get("num_cores", dtype=str)
	num_samplings = ofm.get("num_samplings", dtype=str)
	threshold = ofm.get("threshold", dtype=str)
	filter_enabled = ofm.get("filter_enabled", dtype=bool)
	filt = ofm.get("filter", dtype=str)

	log.info("feature = {0}".format(feature))
	log.info("slice = {0}".format(slice_name))
	log.info("estimator = {0}".format(estimator))
	log.info("num_cores = {0}".format(num_cores))
	log.info("num_samplings = {0}".format(num_samplings))
	log.info("threshold = {0}".format(threshold))
	log.info("filter_enabled = {0}".format(filter_enabled))
	log.info("filter = {0}".format(os.path.basename(filt)))

	cmd = [
		"oncodrivefm-compute",
		"-o", project["temp_path"],
		"-n oncodrivefm-{0}".format(feature),
		"-N", num_samplings,
		"--threshold", threshold,
		"-e {0}".format(estimator),
		"-j", num_cores,
		"--slices '{0}'".format(slice_name)]

	if filter_enabled:
		cmd += ["--filter", filt]

	if feature == "pathways":
		cmd += ["-m", get_data_kegg_path(conf, "ensg_kegg.tsv")]

	cmd += [ofm["data"]]

	project["oncodrivefm"] = dict(
		feature=feature,
		slice=slice_name,
		results=os.path.join(project["temp_path"], "oncodrivefm-{0}-{1}.tsv".format(feature, slice_name)))

	cmd = " ".join(cmd)

	log.debug(cmd)

	ret_code = subprocess.call(cmd, shell=True)
	if ret_code != 0:
		raise Exception("OncodriveFM error while computing {0}:\n{1}".format(feature, cmd))

	projects_out_port.send(project)

task.run()