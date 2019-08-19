import os.path

from wok.task import task

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.variants.utils import var_to_tab
from intogensm.pathutils import ensure_path_exists

@task.foreach()
def variants(project):
	log = task.logger
	conf = task.conf

	partition_port = task.ports("partitions")

	log.info("--- [{0}] --------------------------------------------".format(project["id"]))

	projdb = ProjectDb(project["db"])

	log.info("Preparing variants for VEP ...")

	base_path = os.path.join(project["temp_path"], "vep_ma")

	ensure_path_exists(base_path)

	partition_size = conf.get("vep_partition_size", 100, dtype=int)
	partition = -1
	f = None

	count = 0
	for var in projdb.variants(order_by="position"):
		start, end, ref, alt = var_to_tab(var)

		if count % partition_size == 0:
			if f is not None:
				f.close()

			partition += 1
			partition_path = os.path.join(base_path, "{0:08d}.vep.bed".format(partition))
			f = open(partition_path, "w")
			partition_port.send({
				"project" : project,
				"index" : partition,
				"bed_path" : partition_path,
				"base_path" : base_path})

		tsv.write_line(f, var.chr, start, end, ref + "/" + alt, var.strand, var.id)

		count += 1

	if f is not None:
		f.close()

	log.info("{0} variants prepared".format(count))

	projdb.close()

task.run()