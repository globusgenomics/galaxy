#!/bin/env python

import os
import re
import tarfile
import logging
import sqlite3
import time

from datetime import datetime
from optparse import OptionParser

from bgcore import tsv

from intogensm.varxrefsdb import VarXrefsDb

LOG_LEVEL = {
	"debug" : logging.DEBUG,
	"info" : logging.INFO,
	"warn" : logging.WARN,
	"error" : logging.ERROR,
	"critical" : logging.CRITICAL,
	"notset" : logging.NOTSET }

SUFFIXES = ['K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

def hsize(size):
	if size < 0:
		raise ValueError('number must be non-negative')

	decimals = 0
	multiple = 1000.0
	for suffix in SUFFIXES:
		size /= multiple
		if size < multiple:
			fmt = "{0:." + str(3) + "f} {1}"
			return fmt.format(size, suffix)
		decimals += 3

	raise ValueError('number too large')

def main():
	parser = OptionParser(usage = "usage: %prog [options] <Variants gvf.gz file> ...")

	parser.add_option("--db", dest="db_path",
		help="Database path")

	parser.add_option("-L", "--log-level", dest="log_level",
		default="info", choices=["debug", "info", "warn", "error", "critical", "notset"],
		help="Which log level: debug, info, warn, error, critical, notset")

	(options, args) = parser.parse_args()

	logging.basicConfig(
		level=LOG_LEVEL[options.log_level],
		format="%(asctime)s %(levelname)-5s : %(message)s")

	log = logging.getLogger("var_db")

	if len(args) < 1:
		log.error("At least one variants file is required")
		parser.print_help()
		exit(-1)

	if options.db_path is None:
		log.error("The database path should be specified")
		parser.print_help()
		exit(-1)

	db_path = options.db_path

	log.info("Opening database ...")

	db = VarXrefsDb(db_path)

	db.open()

	db.begin()

	total_count = 0
	total_start_time = time.time()

	src_var_count = {}
	src_ratio = {}

	chromosomes = set()
	chr_var_count = {}
	strands = set()
	
	try:
		partial_count = 0
		partial_start_time = time.time()
		for xref_path in args:
			log.info("Reading {0} ...".format(xref_path))

			if not os.path.isfile(xref_path):
				log.error("File not found: {0}".format(xref_path))
				exit(-1)

			mtime = datetime.fromtimestamp(os.path.getmtime(xref_path))

			f = tsv.open(xref_path, "r")

			src_count = 0
			src_start_time = time.time()

			line_num = 1

			# discard headers
			line = f.readline()
			while line.startswith("#"):
				line = f.readline()
				line_num += 1

			src_var_count[xref_path] = 0

			for line in f:
				try:
					fields = [x if len(x) > 0 else None for x in line.rstrip("\n").split("\t")]

					chr, source, type, start, end, _1, strand, _2, extra = fields

					start = int(start)
					end = int(end)

					ref = None
					alt = None
					xref = None
					try:
						for var in extra.split(";"):
							try:
								key, value = var.split("=")
								if key == "Dbxref":
									pos = value.index(":")
									xref = value[pos + 1:]
								elif key == "Reference_seq":
									ref = value
								elif key == "Variant_seq":
									alt = value
							except:
								continue
					except:
						pass

					if sum([1 if x is None else 0 for x in [chr, start, strand, ref, alt, source, xref]]) > 0:
						log.warn("Discarding incomplete variant: {0}".format(",".join([chr, str(start), strand, ref, alt, source, xref])))
						continue

					src_var_count[xref_path] += 1

					chromosomes.add(chr)
					if chr in chr_var_count:
						chr_var_count[chr] += 1
					else:
						chr_var_count[chr] = 1
					
					strands.add(strand)

					db.add(chr, start, ref, alt, source, xref, strand)

					total_count += 1
					src_count += 1

					partial_count += 1
					elapsed_time = time.time() - partial_start_time
					if elapsed_time >= 10.0:
						ratio = float(partial_count) / elapsed_time
						log.debug("  {0:.1f} variants/second, {1} variants, {2} total variants".format(ratio,
								hsize(src_count), hsize(total_count)))
						partial_count = 0
						partial_start_time = time.time()

				except Exception as ex:
					log.error("Error at line {0}:\n{1}".format(line_num, line.rstrip("\n")))
					import sys
					import traceback
					traceback.print_exc(file=sys.stdout)
					continue
				finally:
					line_num += 1

			elapsed_time = time.time() - src_start_time
			ratio = float(src_count) / elapsed_time
			src_ratio[xref_path] = ratio
			log.info("  {0:.1f} variants/second, {1} variants, {2} total variants".format(ratio,
					hsize(src_count), hsize(total_count)))

			f.close()

		db.commit()
	except KeyboardInterrupt:
		db.commit()
		log.warn("Interrupted by the user with Ctrl-C")
		exit(-1)
	except:
		db.rollback()
		raise
	finally:
		db.close()

	elapsed_time = time.time() - total_start_time
	total_ratio = float(total_count) / elapsed_time

	log.info("Statistics:")
	log.info("  Sources:")
	for xref_path in args:
		log.info("    {0}: {1} variants".format(os.path.basename(xref_path), src_var_count[xref_path]))
	total_size = 0

	log.info("  Chromosomes:")
	for chr in chromosomes:
		log.info("    {0:>2}: {1:>7} variants".format(chr, str(chr_var_count[chr])))
	log.info("  Strands: {0}".format(", ".join(strands)))

	log.info("  Total {0} variants ({1:.1f} variants/sec)".format(hsize(total_count), total_ratio))

if __name__ == "__main__":
	main()