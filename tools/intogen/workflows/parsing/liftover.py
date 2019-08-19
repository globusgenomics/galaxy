import os
import subprocess

from wok.task import task

from bgcore import tsv
from intogensm.woktmp import make_temp_file, remove_temp

from intogensm.projdb import ProjectDb

@task.foreach()
def liftover(project):
	log = task.logger
	conf = task.conf

	lifted_project_port = task.ports("lifted_projects")

	log.info("--- [{0}] --------------------------------------------".format(project["id"]))

	log.info("Preparing liftOver files ...")

	in_path = make_temp_file(task, suffix=".bed")
	in_file = open(in_path, "w")
	out_path = make_temp_file(task, suffix=".bed")
	unmapped_path = os.path.join(project["temp_path"], "liftover_unmapped.bed")

	projdb = ProjectDb(project["db"])

	for var in projdb.variants(order_by="position"):
		in_file.write(tsv.line_text("chr" + var.chr, var.start, var.start + len(var.ref), var.id))

	in_file.close()

	log.info("Running liftOver ...")

	project["from_assembly"] = project["assembly"]
	project["assembly"] = "hg19"

	cmd = " ".join([
		conf["liftover_bin"],
		in_path,
		os.path.join(conf["liftover_chain_path"], "hg18ToHg19.over.chain"),
		out_path,
		unmapped_path
	])

	log.debug(cmd)

	subprocess.call(cmd, shell=True)

	log.info("Annotating unmapped variants ...")

	count = 0
	with open(unmapped_path, "r") as f:
		for line in f:
			if line.lstrip().startswith("#"):
				continue
			fields = line.rstrip().split("\t")
			var_id = int(fields[3])
			projdb.update_variant_start(var_id, start=None)
			count += 1

	log.info("  {0} unmapped variants annotated".format(count))

	log.info("Updating variants ...")

	count = 0
	with open(out_path, "r") as f:
		for line in f:
			fields = line.rstrip().split("\t")
			chr, start, end, var_id = fields
			projdb.update_variant_start(var_id, start=start)
			count += 1

	log.info("  {0} variants".format(count))

	remove_temp(task, in_path, out_path)

	projdb.commit()
	projdb.close()

	lifted_project_port.send(project)

task.run()