import os
import subprocess

from wok.task import task

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.model import Gene, Pathway

@task.foreach()
def compute(project):
	log = task.logger
	conf = task.conf

	projects_out_port = task.ports("projects_out")

	project_id = project["id"]

	ofm = project["oncodrivefm"]
	feature = ofm["feature"]

	log.info("--- [{0} @ {1}] --------------------------------------------".format(project_id, feature))

	cmd = [
		"oncodrivefm-combine",
		"-o", project["temp_path"],
		"-n oncodrivefm-{0}".format(feature)]

	cmd += ofm["data"]

	ofm["results"] = os.path.join(project["temp_path"], "oncodrivefm-{0}.tsv".format(feature))

	cmd = " ".join(cmd)

	log.debug(cmd)

	ret_code = subprocess.call(cmd, shell=True)
	if ret_code != 0:
		raise Exception("OncodriveFM error while combining {0}:\n{1}".format(feature, cmd))

	projects_out_port.send(project)

task.run()