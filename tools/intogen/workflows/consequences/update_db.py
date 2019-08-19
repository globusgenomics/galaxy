import os.path

from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.model import Variant, Consequence, AffectedGene, Gene
from intogensm.varxrefsdb import VarXrefsDb
from intogensm.transfic import TransFIC

from wok.task import task

@task.foreach()
def update_db(project):
	log = task.logger
	conf = task.conf

	projects_port = task.ports("projects_out")

	log.info("--- [{0}] --------------------------------------------".format(project["id"]))

	partitions = project["partitions"]

	if not os.path.exists(conf["vardb_path"]):
		log.warn("Database for variantion external references not found")
		log.debug("> {0}".format(conf["vardb_path"]))

	varxdb = VarXrefsDb(conf["vardb_path"])
	varxdb.open()

	projdb = ProjectDb(project["db"])

	updated_variants = set()

	plen = len(partitions)

	gene_xrefs = {}

	for part in partitions:
		log.info("Updating database with partition data ({0} out of {1}) ...".format(part["index"] + 1, plen))

		log.info("  VEP results ...")

		ctype = lambda v: v.split(",")

		with open(part["vep_path"], "r") as vf:
			for fields in tsv.lines(vf, (int, str, str, ctype, str, str, str, float, float), null_value="-"):
				(var_id, gene, transcript, consequences,
					protein_pos, aa_change, protein,
					sift_score, pph2_score) = fields

				var = projdb.get_variant(var_id)

				xrefs = varxdb.get_xrefs(var.chr, var.start, var.ref, var.alt, var.strand)

				if xrefs is not None:
					xrefs = ["{0}:{1}".format(source, xref) for source, xref in xrefs]
					if gene in gene_xrefs:
						gene_xrefs[gene].update(xrefs)
					else:
						gene_xrefs[gene] = set(xrefs)

					if len(xrefs) == 0:
						xrefs = None

				projdb.update_variant(Variant(id=var_id, xrefs=xrefs))

				projdb.add_consequence(Consequence(var=Variant(id=var_id),
							transcript=transcript, gene=gene, ctypes=consequences,
							protein_pos=protein_pos, aa_change=aa_change, protein=protein))

		log.info("  Variant genes functional impacts ...")

		with open(part["gfi_path"], "r") as f:
			types = (int, str, float, int, str)
			for var_id, gene, impact, coding_region, prot_changes in tsv.lines(f, types, null_value="-"):
				projdb.add_affected_gene(AffectedGene(var=Variant(id=var_id), gene_id=gene, impact=impact,
													  coding_region=coding_region, prot_changes=prot_changes))

		log.info("  Consequence transformed functional impacts ...")

		with open(part["tfi_path"], "r") as f:
			types = (int, str, str, float, float, int, float, float, int, float, float, int, int)
			for fields in tsv.lines(f, types, null_value="-"):
				(var_id, transcript, uniprot,
				 sift_score, sift_tfic, sift_class,
				 pph2_score, pph2_tfic, pph2_class,
				 ma_score, ma_tfic, ma_class, impact) = fields

				projdb.update_consequence(Consequence(var=Variant(id=var_id), transcript=transcript, uniprot=uniprot,
											sift_score=sift_score, sift_tfic=sift_tfic, sift_tfic_class=sift_class,
											pph2_score=pph2_score, pph2_tfic=pph2_tfic, pph2_tfic_class=pph2_class,
											ma_score=ma_score, ma_tfic=ma_tfic, ma_tfic_class=ma_class, impact=impact))

	log.info("Updating database with gene external variant references ...")

	for gene, xrefs in gene_xrefs.items():
		projdb.update_gene(Gene(id=gene, xrefs=xrefs))

	projdb.commit()
	projdb.close()

	varxdb.close()

	del project["partitions"]

	projects_port.send(project)

task.run()