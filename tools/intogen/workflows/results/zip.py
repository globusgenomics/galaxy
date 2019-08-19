import os.path

from datetime import datetime

from wok.task import task

from bgcore import tsv

from intogensm.woktmp import make_temp_file, remove_temp

from intogensm.projdb import ProjectDb
from intogensm.sigdb import SigDb
from intogensm.transfic import TransFIC
from intogensm.variants.utils import var_to_tab
from intogensm.archive import Archive
from intogensm.projres import ProjectResults
from intogensm import so

def write_line(f, *v):
	tsv.write_line(f, *v, null_value="-")

class ArcFile(object):
	def __init__(self, task, arc, project_id, name, mode):
		self.task = task
		self.arc = arc
		self.project_id = project_id
		self.name = name
		self.mode = mode
		name, ext = os.path.splitext(name)
		if len(ext) == 0:
			self.arcname = "{0}/{1}.tsv".format(self.project_id, self.name)
		else:
			self.arcname = "{0}/{1}".format(self.project_id, self.name)

	def __enter__(self):
		self.tmp = make_temp_file(self.task, suffix=".tsv")
		task.logger.debug("    > {0}".format(self.tmp))
		self.f = open(self.tmp, self.mode)
		return self.f

	def __exit__(self, exc_type, exc_val, exc_tb):
		self.f.close()
		if exc_val is None:
			self.arc.add(self.tmp, arcname=self.arcname)
		remove_temp(self.task, self.tmp)
		return False

@task.foreach()
def pack_results(project):
	log = task.logger
	conf = task.conf

	project_id = project["id"]

	log.info("--- [{0}] --------------------------------------------".format(project_id))

	project_path = project["path"]
	temp_path = project["temp_path"]

	dest_path = os.path.join(project_path, "results.zip")

	sigdb = SigDb(conf["sigdb_path"])
	sigdb.open()

	projdb = ProjectDb(project["db"])

	projres = ProjectResults(project)

	gene_sym = projdb.get_gene_symbols()

	total_samples = projdb.get_total_affected_samples()

	log.info("Compressing files ...")

	arc = None
	try:
		arc = Archive(dest_path, mode="w", fmt="zip")

		log.info("  Variant genes ...")

		with ArcFile(task, arc, project_id, "variant_genes", "w") as vf:
			write_line(vf, "PROJECT_ID", "CHR", "STRAND", "START", "ALLELE",
							"GENE_ID", "SYMBOL", "VAR_IMPACT", "VAR_IMPACT_DESC",
							"SAMPLE_FREQ", "SAMPLE_TOTAL", "SAMPLE_PROP",
							"CODING_REGION", "PROTEIN_CHANGES", "INTOGEN_DRIVER", "XREFS")

			for afg in projdb.affected_genes(join_variant=True, join_xrefs=True, join_rec=True):
				var = afg.var
				rec = afg.rec

				start, end, ref, alt = var_to_tab(var)

				xrefs = [xref for xref in var.xrefs]
				if sigdb.exists_variant(var.chr, start):
					xrefs += ["I:1"]
				xrefs = ",".join(xrefs)

				intogen_driver = 1 if sigdb.exists_gene(afg.gene_id) else 0

				write_line(vf, project_id, var.chr, var.strand, start, "{0}/{1}".format(ref, alt),
								afg.gene_id, gene_sym.get(afg.gene_id),
								afg.impact, TransFIC.class_name(afg.impact),
								rec.sample_freq or 0, total_samples, rec.sample_prop or 0,
								afg.coding_region, afg.prot_changes, intogen_driver, xrefs)

		log.info("  Variant samples ...")

		with ArcFile(task, arc, project_id, "variant_samples", "w") as vf:
			write_line(vf, "PROJECT_ID", "CHR", "STRAND", "START", "ALLELE", "SAMPLES")

			for var in projdb.variants(join_samples=True):
				start, end, ref, alt = var_to_tab(var)
				write_line(vf, project_id, var.chr, var.strand, start, "{0}/{1}".format(ref, alt),
						   ",".join([s.name for s in var.samples]))

		log.info("  Consequences ...")

		with ArcFile(task, arc, project_id, "consequences", "w") as cf:
			write_line(cf, "PROJECT_ID", "CHR", "STRAND", "START", "ALLELE", "TRANSCRIPT_ID", "CT",
					   		"GENE_ID", "SYMBOL", "UNIPROT_ID", "PROTEIN_ID", "PROTEIN_POS", "AA_CHANGE",
							"SIFT_SCORE", "SIFT_TRANSFIC", "SIFT_TRANSFIC_CLASS",
							"PPH2_SCORE", "PPH2_TRANSFIC", "PPH2_TRANSFIC_CLASS",
							"MA_SCORE", "MA_TRANSFIC", "MA_TRANSFIC_CLASS",
							"IMPACT", "IMPACT_CLASS")

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

				write_line(cf, project_id, var.chr, var.strand, start, allele, csq.transcript,
							",".join(csq.ctypes), csq.gene, gene_sym.get(csq.gene),
							uniprot, protein, protein_pos, aa_change,
							sift_score, sift_tfic, sift_tfic_class,
							pph2_score, pph2_tfic, pph2_tfic_class,
							ma_score, ma_tfic, ma_tfic_class,
							csq.impact, TransFIC.class_name(csq.impact))

		log.info("  Genes ...")

		with ArcFile(task, arc, project_id, "genes", "w") as gf:
			write_line(gf, "PROJECT_ID", "GENE_ID", "SYMBOL", "FM_PVALUE", "FM_QVALUE", "FM_EXC_CAUSE",
							"SAMPLE_FREQ", "SAMPLE_TOTAL", "SAMPLE_PROP",
							"CLUST_ZSCORE", "CLUST_PVALUE", "CLUST_QVALUE", "CLUST_EXC_CAUSE", "CLUST_COORDS",
							"INTOGEN_DRIVER", "XREFS")

			for gene in projdb.genes(join_xrefs=True, join_rec=True):
				if gene.rec.sample_freq is not None and gene.rec.sample_freq > 0:
					intogen_driver = 1 if sigdb.exists_gene(gene.id) else 0
					write_line(gf, project_id, gene.id, gene.symbol, gene.fm_pvalue, gene.fm_qvalue, gene.fm_exc_cause,
									gene.rec.sample_freq, total_samples, gene.rec.sample_prop or 0,
									gene.clust_zscore, gene.clust_pvalue, gene.clust_qvalue, gene.clust_exc_cause,
									gene.clust_coords, intogen_driver, ",".join(gene.xrefs))

		log.info("  Pathways ...")

		with ArcFile(task, arc, project_id, "pathways", "w") as pf:
			write_line(pf, "PROJECT_ID", "PATHWAY_ID", "GENE_COUNT", "FM_ZSCORE", "FM_PVALUE", "FM_QVALUE",
							"SAMPLE_FREQ", "SAMPLE_TOTAL", "SAMPLE_PROP", "GENE_FREQ", "GENE_TOTAL", "GENE_PROP")

			for pathway in projdb.pathways(join_rec=True):
				if pathway.rec.sample_freq is not None and pathway.rec.sample_freq > 0:
					write_line(pf, project_id, pathway.id, pathway.gene_count, pathway.fm_zscore, pathway.fm_pvalue, pathway.fm_qvalue,
									pathway.rec.sample_freq or 0, total_samples, pathway.rec.sample_prop or 0,
									pathway.rec.gene_freq or 0, pathway.gene_count, pathway.rec.gene_prop or 0)

		skip_oncodrivefm = conf.get("skip_oncodrivefm", False, dtype=bool)

		if not skip_oncodrivefm:

			log.info("  Genes per sample functional impact ...")

			with ArcFile(task, arc, project_id, "fimpact.gitools.tdm", "w") as f:
				write_line(f, "SAMPLE", "GENE_ID",
						   "SIFT_SCORE", "SIFT_TRANSFIC", "SIFT_TRANSFIC_CLASS",
						   "PPH2_SCORE", "PPH2_TRANSFIC", "PPH2_TRANSFIC_CLASS",
						   "MA_SCORE", "MA_TRANSFIC", "MA_TRANSFIC_CLASS")
				for fields in projdb.sample_gene_fimpacts():
					(gene, sample,
						 sift_score, sift_tfic, sift_tfic_class,
						 pph2_score, pph2_tfic, pph2_tfic_class,
						 ma_score, ma_tfic, ma_tfic_class) = fields
					write_line(f, sample, gene,
							   sift_score, sift_tfic, TransFIC.class_name(sift_tfic_class),
							   pph2_score, pph2_tfic, TransFIC.class_name(pph2_tfic_class),
							   ma_score, ma_tfic, TransFIC.class_name(ma_tfic_class))

		log.info("Saving project configuration ...")

		with ArcFile(task, arc, project_id, "project", "w") as f:
			names = ["PROJECT_ID", "ASSEMBLY", "SAMPLES_TOTAL"]
			values = [project_id, project["assembly"], total_samples]
			names, values = projres.get_annotations_to_save(conf, project["annotations"], names=names, values=values)
			tsv.write_line(f, *names)
			tsv.write_line(f, *values, null_value="-")
	finally:
		if arc is not None:
			arc.close()
		projdb.close()
		sigdb.close()

task.run()