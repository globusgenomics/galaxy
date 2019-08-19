import os
import re
from datetime import datetime

from wok.task import task

from bgcore import tsv

from intogensm.sigdb import SigDb
from intogensm.paths import get_results_path, get_combination_path

def add_cancer_site(gene_sites, gene, code, name):
	if gene not in gene_sites:
		gene_sites[gene] = set([(code, name)])
	else:
		gene_sites[gene].add((code, name))

@task.begin()
def drivers():
	log = task.logger
	conf = task.conf

	db_path = get_results_path(conf, "drivers.db")
	db = SigDb(db_path)
	db.open()

	log.info("Variants ...")

	path = get_combination_path(conf, "recurrences", "variant_gene-global-all.tsv.gz")
	with tsv.open(path, "r") as f:
		types = (str, str, int, str)
		for fields in tsv.lines(f, types, columns=("CHR", "STRAND", "START", "ALLELE"), header=True):
			chr, strand, start, allele = fields[:4]
			db.add_variant(chr, start)

	log.info("Genes ...")

	gene_sites = {}

	gene_fm = set()
	gene_clust = set()

	#SPECIAL_THRESHOLD = ["C18", "C34"]
	SPECIAL_THRESHOLD = []

	log.info("  OncodriveFM ...")

	filename_re = re.compile(r"gene-cancer_site-(.+)\.tsv.gz")
	base_path = get_combination_path(conf, "oncodrivefm")
	for path in os.listdir(base_path):
		m = filename_re.match(path)
		if not m:
			continue

		cancer_site_code = m.group(1)

		if cancer_site_code in SPECIAL_THRESHOLD:
			threshold = 1e-6
		else:
			threshold = 0.01

		with tsv.open(os.path.join(base_path, path), "r") as f:
			params = tsv.params(f)
			cancer_site_name = params["group_long_name"]
			for fields in tsv.lines(f, (str, float), columns=("ID", "QVALUE"), header=True):
				gene, qvalue = fields
				if qvalue < threshold:
					add_cancer_site(gene_sites, gene, cancer_site_code, cancer_site_name)

					gene_fm.add(gene)

	log.info("  OncodriveCLUST ...")

	filename_re = re.compile(r"cancer_site-(.+)\.tsv.gz")
	base_path = get_combination_path(conf, "oncodriveclust")
	for path in os.listdir(base_path):
		m = filename_re.match(path)
		if not m:
			continue

		cancer_site_code = m.group(1)

		with tsv.open(os.path.join(base_path, path), "r") as f:
			params = tsv.params(f)
			cancer_site_name = params["group_long_name"]
			for fields in tsv.lines(f, (str, float), columns=("ID", "QVALUE"), header=True):
				gene, qvalue = fields
				if qvalue < 0.05:
					add_cancer_site(gene_sites, gene, cancer_site_code, cancer_site_name)

					gene_clust.add(gene)

	log.info("  Updating db ...")
	sig_genes = gene_fm | gene_clust
	for gene in sig_genes:
		db.add_gene(gene, gene in gene_fm, gene in gene_clust)

	log.info("Saving driver genes cancer sites dataset ...")
	path = get_results_path(conf, "gene-driver_cancer_sites.tsv")
	log.debug("> {}".format(path))
	with open(path, "w") as f:
		tsv.write_param(f, "date", datetime.now())
		tsv.write_line(f, "GENE_ID", "FM", "CLUST", "CANCER_SITES_COUNT", "CANCER_SITE_CODES", "CANCER_SITE_NAMES")
		for gene, sites in gene_sites.items():
			tsv.write_line(f, gene,
						   1 if gene in gene_fm else 0,
						   1 if gene in gene_clust else 0,
						   len(sites),
						   ", ".join(sorted([code for code, name in sites])),
						   ", ".join(sorted([name for code, name in sites])))

	db.commit()
	db.close()

task.run()