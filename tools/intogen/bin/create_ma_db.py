#!/bin/env python

import os
import re
import tarfile
import logging
import sqlite3
import time

from datetime import datetime
from optparse import OptionParser

from intogensm.ma.db import MaDb

NAME_RE = re.compile(r"^MA\.hg19/MA\.chr(.+)\.txt$")

LOG_LEVEL = {
	"debug" : logging.DEBUG,
	"info" : logging.INFO,
	"warn" : logging.WARN,
	"error" : logging.ERROR,
	"critical" : logging.CRITICAL,
	"notset" : logging.NOTSET }

SUFFIXES = ['K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

def approximate_size(size):
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
	parser = OptionParser(usage = "usage: %prog [options] <MA tar file> <output db file>")

	parser.add_option("-s", "--size", dest="size", type="int", default=10000,
		help="Partition size")

	parser.add_option("-L", "--log-level", dest="log_level",
		default="info", choices=["debug", "info", "warn", "error", "critical", "notset"],
		help="Which log level: debug, info, warn, error, critical, notset")

	(options, args) = parser.parse_args()

	logging.basicConfig(
		level=LOG_LEVEL[options.log_level],
		format="%(asctime)s %(levelname)-5s : %(message)s")

	log = logging.getLogger("ma_db")

	if len(args) != 2:
		log.error("Wrong number of arguments")
		parser.print_help()
		exit(-1)

	ma_path = args[0]
	if not os.path.isfile(ma_path):
		log.error("File not found: {0}".format(ma_path))
		exit(-1)

	ma_time = datetime.fromtimestamp(os.path.getmtime(ma_path))

	db_path = args[1]

	size = options.size

	log.info("Opening database ...")

	db = MaDb(db_path, part_size=size, source_time=ma_time)

	db.open()

	log.info("Opening {0} ...".format(ma_path))

	tar = tarfile.open(ma_path, "r|*")

	total_count = 0
	total_start_time = time.time()

	chr_ratio = {}

	try:
		partial_count = 0
		partial_start_time = time.time()
		for member in tar:
			if not member.isfile():
				continue

			m = NAME_RE.match(member.name)
			if m is None:
				continue

			chr = m.group(1)
			try:
				chr = str(int(chr))
			except:
				pass

			log.info("Processing chromosome {0} ...".format(chr))

			f = tar.extractfile(member)

			chr_count = 0
			chr_start_time = time.time()

			f.readline() # discard header

			for line in f:
				fields = [x if len(x) > 0 else None for x in line.rstrip("\n").split("\t")]

				mutation, ref_genome_variant, gene, uniprot, info, uniprot_variant, func_impact, fi_score = fields

				assembly, chr, start, ref, alt = mutation.split(",")

				start = int(start)

				if fi_score is not None:
					fi_score = float(fi_score)

				db.add(chr, start, ref, alt, ref_genome_variant, gene, uniprot, info, uniprot_variant, func_impact, fi_score)

				total_count += 1
				chr_count += 1

				partial_count += 1
				elapsed_time = time.time() - partial_start_time
				if elapsed_time >= 10.0:
					ratio = float(partial_count) / elapsed_time
					log.debug("  {0:.1f} mutations/second, {1} chromosome mutations, {2} total mutations".format(ratio,
							approximate_size(chr_count), approximate_size(total_count)))
					partial_count = 0
					partial_start_time = time.time()

			elapsed_time = time.time() - chr_start_time
			ratio = float(chr_count) / elapsed_time
			chr_ratio[chr] = ratio
			log.info("  {0:.1f} mutations/second, {1} chromosome mutations, {2} total mutations".format(ratio,
					approximate_size(chr_count), approximate_size(total_count)))

			db.commit()
	except KeyboardInterrupt:
		db.rollback()
		log.warn("Interrupted by the user with Ctrl-C")

	elapsed_time = time.time() - total_start_time
	total_ratio = float(total_count) / elapsed_time

	chromosomes, chr_parts, part_size = db.stats()

	log.info("Statistics:")
	total_size = 0
	for chr in chromosomes:
		log.info("  Chromosome {0}".format(chr))
		chr_size = 0
		for part_name in chr_parts[chr]:
			size = part_size[part_name]
			chr_size += size
			total_size += size
			log.debug("    {0}: {1} mutations".format(part_name, approximate_size(size)))

		if chr in chr_ratio:
			ratio = "({0:.1f} mutations/sec)".format(chr_ratio[chr])
		else:
			ratio = ""

		log.info("    Total {0} mutations {1}".format(approximate_size(chr_size), ratio))

	log.info("  Total {0} mutations ({1:.1f} mutations/sec)".format(approximate_size(total_size), total_ratio))

	tar.close()
	db.close()

if __name__ == "__main__":
	main()