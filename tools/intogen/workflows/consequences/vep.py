import os

from wok.task import task

from bgcore import tsv

from intogensm.vep.local import VepLocal
from intogensm.vep.service import VepService

@task.foreach()
def vep_run(partition):
	log = task.logger
	conf = task.conf

	results_port = task.ports("results")

	project_id = partition["project"]["id"]
	log.info("--- [{0} @ {1}] --------------------------------------------".format(project_id, partition["index"]))

	offline = conf["offline"]

	if offline == "yes":
		log.info("Running VEP in local mode.")

		vep = VepLocal(
				perl_path=conf["perl_bin"],
				lib_path=conf["perl_lib"],
				script_path=os.path.join(conf["ext_bin_path"], "variant_effect_predictor", "variant_effect_predictor.pl"),
				cache_path=os.path.join(conf["data_path"], "vep_cache"))
	else:
		log.info("Running VEP using web services.")

		vep = VepService(cache_path=os.path.join(conf["cache_path"], "vep.db"))

	results_path = os.path.join(partition["base_path"], "{0:08d}.vep".format(partition["index"]))

	if not os.path.exists(results_path) or conf.get("consequences_overwrite", True):
		# Run VEP
		vep.run(partition["bed_path"])

		log.info("Saving results ...")
		log.debug("VEP results: {0}".format(vep.results_path))

		# Save results

		with open(results_path, "w") as f:
			for r in vep.results():
				tsv.write_line(f, r.var_id, r.gene, r.transcript,
								",".join(r.consequences),
								r.protein_pos, r.aa_change, r.protein,
								r.sift, r.polyphen, null_value="-")

	else:
		log.warn("Skipping VEP, results already exist.")
		log.debug("VEP results: {0}".format(results_path))

	vep.close()

	# Send results to the next module
	partition["vep_path"] = results_path
	results_port.send(partition)

task.run()