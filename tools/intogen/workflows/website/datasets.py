import os
import gzip
import bz2

from datetime import datetime

from wok.task import task

from bgcore import tsv

from intogensm import VERSION
from intogensm.projdb import ProjectDb
from intogensm.sigdb import SigDb
from intogensm.model import Consequence
from intogensm.transfic import TransFIC
from intogensm.variants.utils import var_to_tab
from intogensm.paths import get_website_results_path
from intogensm.projres import ProjectResults
from intogensm import so

def open_dataset(project_id, base_path, datasets_path, name, mode, logger):
	name, ext = os.path.splitext(name)
	ext = ext.lower()
	if len(ext) == 0:
		ext = ".gz"
		name = "{0}.tsv{1}".format(name, ext)
	else:
		name = name + ext

	path = os.path.join(datasets_path, name)
	logger.debug("> {0}".format(os.path.relpath(path, base_path)))

	f = tsv.open(path, mode)

	tsv.write_param(f, "version", VERSION)
	tsv.write_param(f, "date", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
	tsv.write_param(f, "PROJECT_ID", project_id)

	return f

@task.foreach()
def datasets(project):
	log = task.logger
	conf = task.conf

	project_id = project["id"]

	log.info("--- [{0}] --------------------------------------------".format(project_id))

	project_path = project["path"]
	temp_path = project["temp_path"]

	datasets_path = get_website_results_path(project_path)
	if not os.path.exists(datasets_path):
		os.makedirs(datasets_path)

	sigdb = SigDb(conf["sigdb_path"])
	sigdb.open()

	projdb = ProjectDb(project["db"])

	gene_sym = projdb.get_gene_symbols()

	total_samples = projdb.get_total_affected_samples()

	log.info("Exporting variant genes ...")

	vf = open_dataset(project_id, project_path, datasets_path, "variant_gene", "w", log)
	tsv.write_param(vf, "SAMPLE_TOTAL", total_samples)
	tsv.write_line(vf, "VAR_ID", "CHR", "STRAND", "START", "ALLELE",
					"GENE_ID", "IMPACT", "IMPACT_CLASS",
					"SAMPLE_FREQ", "SAMPLE_PROP",
					"CODING_REGION", "PROTEIN_CHANGES", "INTOGEN_DRIVER", "XREFS")

	sf = open_dataset(project_id, project_path, datasets_path, "variant-samples", "w", log)
	tsv.write_line(sf, "VAR_ID", "CHR", "STRAND", "START", "ALLELE", "SAMPLE")

	count = 0
	for afg in projdb.affected_genes(join_variant=True, join_samples=True, join_xrefs=True, join_rec=True):
		var = afg.var
		rec = afg.rec

		start, end, ref, alt = var_to_tab(var)

		allele = "{0}/{1}".format(ref, alt)

		xrefs = [xref for xref in var.xrefs]
		if sigdb.exists_variant(var.chr, start):
			xrefs += ["I:1"]
		xrefs = ",".join(xrefs)

		intogen_driver = 1 if sigdb.exists_gene(afg.gene_id) else 0

		tsv.write_line(vf, var.id, var.chr, var.strand, start, allele,
						afg.gene_id, afg.impact, TransFIC.class_name(afg.impact),
						rec.sample_freq, rec.sample_prop,
						afg.coding_region, afg.prot_changes, intogen_driver, xrefs, null_value="\N")

		for sample in var.samples:
			tsv.write_line(sf, var.id, var.chr, var.strand, start, allele, sample.name, null_value="\N")

		count += 1

	vf.close()
	sf.close()

	log.info("  {0} variant genes".format(count))

	log.info("Exporting consequences ...")

	cf = open_dataset(project_id, project_path, datasets_path, "consequence", "w", log)
	tsv.write_line(cf, "VAR_ID", "CHR", "STRAND", "START", "ALLELE", "TRANSCRIPT_ID",
				   "CT", "GENE_ID", "SYMBOL", "UNIPROT_ID", "PROTEIN_ID", "PROTEIN_POS", "AA_CHANGE",
					"SIFT_SCORE", "SIFT_TRANSFIC", "SIFT_TRANSFIC_CLASS",
					"PPH2_SCORE", "PPH2_TRANSFIC", "PPH2_TRANSFIC_CLASS",
					"MA_SCORE", "MA_TRANSFIC", "MA_TRANSFIC_CLASS",
					"IMPACT", "IMPACT_CLASS")

	count = 0
	for csq in projdb.consequences(join_variant=True):
		var = csq.var
		start, end, ref, alt = var_to_tab(var)

		allele = "{0}/{1}".format(ref, alt)

		uniprot = protein = protein_pos = aa_change = None
		sift_score = sift_tfic = sift_tfic_class = None
		pph2_score = pph2_tfic = pph2_tfic_class = None
		ma_score = ma_tfic = ma_tfic_class = None

		if so.match(csq.ctypes, so.ONCODRIVEFM):
			uniprot, protein = csq.uniprot, csq.protein

		if so.match(csq.ctypes, so.NON_SYNONYMOUS):
			protein_pos, aa_change = csq.protein_pos, csq.aa_change
			sift_score, sift_tfic, sift_tfic_class = csq.sift_score, csq.sift_tfic, TransFIC.class_name(csq.sift_tfic_class)
			pph2_score, pph2_tfic, pph2_tfic_class = csq.pph2_score, csq.pph2_tfic, TransFIC.class_name(csq.pph2_tfic_class)
			ma_score, ma_tfic, ma_tfic_class = csq.ma_score, csq.ma_tfic, TransFIC.class_name(csq.ma_tfic_class)

		tsv.write_line(cf, var.id, var.chr, var.strand, start, allele, csq.transcript,
						",".join(csq.ctypes), csq.gene, gene_sym.get(csq.gene),
						uniprot, protein, protein_pos, aa_change,
						sift_score, sift_tfic, sift_tfic_class,
						pph2_score, pph2_tfic, pph2_tfic_class,
						ma_score, ma_tfic, ma_tfic_class,
						csq.impact, TransFIC.class_name(csq.impact), null_value="\N")
		count += 1

	cf.close()

	log.info("  {0} consequences".format(count))

	log.info("Exporting genes ...")

	gf = open_dataset(project_id, project_path, datasets_path, "gene", "w", log)
	tsv.write_param(gf, "SAMPLE_TOTAL", total_samples)
	tsv.write_line(gf, "GENE_ID", "FM_PVALUE", "FM_QVALUE", "FM_EXC_CAUSE",
				   "CLUST_ZSCORE", "CLUST_PVALUE", "CLUST_QVALUE", "CLUST_EXC_CAUSE", "CLUST_COORDS",
				   "SAMPLE_FREQ", "SAMPLE_PROP", "INTOGEN_DRIVER")


	for gene in projdb.genes(join_rec=True):
		rec = gene.rec

		if rec.sample_freq is None or rec.sample_freq == 0:
			continue

		intogen_driver = 1 if sigdb.exists_gene(gene.id) else 0

		tsv.write_line(gf, gene.id, gene.fm_pvalue, gene.fm_qvalue, gene.fm_exc_cause,
					   gene.clust_zscore, gene.clust_pvalue, gene.clust_qvalue, gene.clust_exc_cause, gene.clust_coords,
					   rec.sample_freq or 0, rec.sample_prop or 0,
					   intogen_driver, null_value="\N")

	gf.close()

	log.info("Exporting pathways ...")

	pf = open_dataset(project_id, project_path, datasets_path, "pathway", "w", log)
	tsv.write_param(pf, "SAMPLE_TOTAL", total_samples)
	tsv.write_line(pf, "PATHWAY_ID", "GENE_COUNT", "FM_ZSCORE", "FM_PVALUE", "FM_QVALUE",
				   "SAMPLE_FREQ", "SAMPLE_PROP", "GENE_FREQ", "GENE_TOTAL", "GENE_PROP")

	for pathway in projdb.pathways(join_rec=True):
		rec = pathway.rec

		if rec.sample_freq is None or rec.sample_freq == 0:
			continue

		tsv.write_line(pf, pathway.id, pathway.gene_count, pathway.fm_zscore, pathway.fm_pvalue, pathway.fm_qvalue,
						rec.sample_freq or 0, rec.sample_prop or 0, rec.gene_freq or 0, pathway.gene_count, rec.gene_prop or 0, null_value="\N")

	pf.close()

	skip_oncodrivefm = conf.get("skip_oncodrivefm", False, dtype=bool)

	if not skip_oncodrivefm:

		log.info("Exporting genes per sample functional impact ...")

		with open_dataset(project_id, project_path, datasets_path, "gene_sample-fimpact", "w", log) as f:
			tsv.write_line(f, "GENE_ID", "SAMPLE",
					   "SIFT_SCORE", "SIFT_TRANSFIC", "SIFT_TRANSFIC_CLASS",
					   "PPH2_SCORE", "PPH2_TRANSFIC", "PPH2_TRANSFIC_CLASS",
					   "MA_SCORE", "MA_TRANSFIC", "MA_TRANSFIC_CLASS")

			for fields in projdb.sample_gene_fimpacts():
				(gene, sample,
					sift_score, sift_tfic, sift_tfic_class,
					pph2_score, pph2_tfic, pph2_tfic_class,
					ma_score, ma_tfic, ma_tfic_class) = fields
				tsv.write_line(f, gene, sample,
						   sift_score, sift_tfic, TransFIC.class_name(sift_tfic_class),
						   pph2_score, pph2_tfic, TransFIC.class_name(pph2_tfic_class),
						   ma_score, ma_tfic, TransFIC.class_name(ma_tfic_class), null_value="\N")

	projdb.close()

	sigdb.close()

	log.info("Saving project configuration ...")

	projres = ProjectResults(project)

	with open_dataset(project_id, project_path, datasets_path, "project.tsv", "w", log) as f:
		names = ["ASSEMBLY", "SAMPLES_TOTAL"]
		values = [project["assembly"], total_samples]
		names, values = projres.get_annotations_to_save(conf, project["annotations"], names=names, values=values)
		tsv.write_line(f, *names)
		tsv.write_line(f, *values, null_value="\N")

	projects_port = task.ports("projects_out")
	projects_port.send(project)

task.run()