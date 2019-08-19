import os.path

from bgcore import tsv

from intogensm.ma.local import MaLocal
from intogensm.ma.service import MaService
from intogensm.projdb import ProjectDb
from intogensm.variants.utils import var_to_tab
from intogensm import so

from wok.task import task

@task.foreach()
def ma_run(partition):
	log = task.logger
	conf = task.conf

	results_port = task.ports("results")

	project = partition["project"]

	log.info("--- [{0} @ {1}] --------------------------------------------".format(project["id"], partition["index"]))

	offline = conf["offline"]

	if offline == "yes":
		log.info("Running Mutation assessor in local mode.")

		ma = MaLocal(conf["ma_cache_path"])
	else:
		log.info("Running Mutation assessor using web services.")

		ma = MaService(project["assembly"], cache_path=os.path.join(conf["cache_path"], "ma.db"))

	results_path = os.path.join(partition["base_path"], "{0:08d}.ma".format(partition["index"]))

	if not os.path.exists(results_path) or conf.get("consequences_overwrite", True):

		log.info("Querying Mutation assessor for 'missense_variant' consequences ...")

		projdb = ProjectDb(project["db"])

		missense_variants = set()

		with open(partition["vep_path"], "r") as f:
			for line in f:
				fields = line.rstrip().split("\t")
				var_id = int(fields[0])
				ctypes = fields[3].split(",")
				if so.match(ctypes, so.NON_SYNONYMOUS):
					missense_variants.add(var_id)

		with open(results_path, "w") as mf:
			for var_id in missense_variants:
				var = projdb.get_variant(var_id)

				start, end, ref, alt = var_to_tab(var)

				r = ma.get(var.chr, var.strand, start, ref, alt, var_id)
				if r is not None:
					tsv.write_line(mf, var_id, r.uniprot, r.fi_score, null_value="-")

		projdb.close()

	else:
		log.warn("Skipping MA, results already exist.")
		log.debug("MA results: {0}".format(results_path))

	ma.close()

	# Send results to the next module
	partition["ma_path"] = results_path
	results_port.send(partition)

task.run()