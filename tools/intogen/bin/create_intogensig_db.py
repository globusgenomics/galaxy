#!/bin/env python

import os
import re
import argparse
import logging

from bgcore import tsv

from intogensm.sigdb import SigDb

_LOG_LEVELS = {
	"debug" : logging.DEBUG,
	"info" : logging.INFO,
	"warn" : logging.WARN,
	"error" : logging.ERROR,
	"critical" : logging.CRITICAL,
	"notset" : logging.NOTSET }

def main():
	parser = argparse.ArgumentParser(
		description="Create a database with the significant features found in IntOGen Mutations")

	parser.add_argument("source_path", metavar="SOURCE",
						help="The folder with the results")

	parser.add_argument("--db", dest="db_path", metavar="PATH",
						help="Database path")

	parser.add_argument("-L", "--log-level", dest="log_level", metavar="LEVEL", default=None,
						choices=["debug", "info", "warn", "error", "critical", "notset"],
						help="Define log level: debug, info, warn, error, critical, notset")

	args = parser.parse_args()

	logging.basicConfig(
		format = "%(asctime)s %(name)s %(levelname)-5s : %(message)s",
		datefmt = "%Y-%m-%d %H:%M:%S")

	if args.log_level is None:
		args.log_level = "info"
	else:
		args.log_level = args.log_level.lower()

	log = logging.getLogger("intogen-sig")
	log.setLevel(_LOG_LEVELS[args.log_level])

	if args.db_path is None:
		args.db_path = "intogen-sig.db"

	db = SigDb(args.db_path)
	db.open()

	log.info("Variant genes ...")

	path = os.path.join(args.source_path, "combination", "recurrences", "variant_gene-global-all.tsv.gz")
	with tsv.open(path, "r") as f:
		types = (str, str, int, str)
		for fields in tsv.lines(f, types, columns=("CHR", "STRAND", "START", "ALLELE"), header=True):
			chr, strand, start, allele = fields[:4]
			db.add_variant(chr, start)

	log.info("Genes ...")

	gene_fm = set()
	gene_clust = set()

	SPECIAL_THRESHOLD = ["C18", "C34"]

	log.info("  OncodriveFM ...")

	filter = re.compile(r"gene-cancer_site-(.+)\.tsv.gz")
	base_path = os.path.join(args.source_path, "combination", "oncodrivefm")
	for path in os.listdir(base_path):
		m = filter.match(path)
		if not m:
			continue

		if m.group(1) in SPECIAL_THRESHOLD:
			threshold = 1e-6
		else:
			threshold = 0.01

		with tsv.open(os.path.join(base_path, path), "r") as f:
			for fields in tsv.lines(f, (str, float), columns=("ID", "QVALUE"), header=True):
				gene, qvalue = fields
				if qvalue < threshold:
					gene_fm.add(gene)

	log.info("  OncodriveCLUST ...")

	filter = re.compile(r"cancer_site-(.+)\.tsv.gz")
	base_path = os.path.join(args.source_path, "combination", "oncodriveclust")
	for path in os.listdir(base_path):
		m = filter.match(path)
		if not m:
			continue

		with tsv.open(os.path.join(base_path, path), "r") as f:
			for fields in tsv.lines(f, (str, float), columns=("ID", "QVALUE"), header=True):
				gene, qvalue = fields
				if qvalue < 0.05:
					gene_clust.add(gene)

	log.info("  Updating db ...")
	sig_genes = gene_fm | gene_clust
	for gene in sig_genes:
		db.add_gene(gene, gene in gene_fm, gene in gene_clust)

	db.commit()
	db.close()

if __name__ == "__main__":
	main()