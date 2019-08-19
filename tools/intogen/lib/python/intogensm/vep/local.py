import os
import re
import tempfile
import subprocess

from model import VepException, VepResult

_SIFT_POLYPHEN_RE = re.compile(r"^.+\((.+)\)$")

class VepLocal(object):
	def __init__(self, perl_path, lib_path, script_path, cache_path):
		self.perl_path = perl_path
		self.lib_path = lib_path
		self.script_path = script_path
		self.cache_path = cache_path
		self.results_path = None

	def run(self, variants_path):
		"""
		Run the VEP script and save results in a temporary file.

		:param variants_path: File with variants. In BED format. http://www.ensembl.org/info/docs/variation/vep/vep_script.html#custom_formats
		:return: True if successfull or False otherwise
		"""

		if self.results_path is None:
			self.results_path = tempfile.mkstemp()[1]

		self.cmd = " ".join([
			self.perl_path,
			self.script_path,
			"-i", variants_path,
			"-o", self.results_path,
			"--no_progress --compress 'gunzip -c'",
			"--force_overwrite --format=ensembl",
			"--cache --offline --dir={0}".format(self.cache_path),
			"--no_intergenic",
			"--protein --sift=b --polyphen=b"])

		self.retcode = subprocess.call(self.cmd, shell=True, env={"PERL5LIB" : self.lib_path, "COLUMNS" : "20", "LINES" : "10"})

		if self.retcode != 0:
			raise VepException("Error while running VEP:\n{0}".format(self.cmd))

	def results(self):
		"""
		Iterator that parses the results temporary file and yields VepResult's
		"""

		with open(self.results_path, "r") as f:
			for line in f:
				if line.startswith("#"):
					continue

				fields = line.rstrip("\n").split("\t")

				var_id, location, allele, gene, feature, feature_type, consequences = fields[0:7]

				if feature_type != "Transcript":
					continue

				chr, start = location.split(":")

				consequences = sorted(consequences.split(","))

				protein_pos, aa_change = [x if x != "-" else None for x in fields[9:11]]

				extra = fields[13]

				protein, sift, polyphen = (None, None, None)

				if len(extra) > 0:
					for var in extra.split(";"):
						key, value = var.split("=")
						if key == "ENSP":
							protein = value
						elif key == "SIFT":
							m = _SIFT_POLYPHEN_RE.match(value)
							if m is not None:
								sift = float(m.group(1))
						elif key == "PolyPhen":
							m = _SIFT_POLYPHEN_RE.match(value)
							if m is not None:
								polyphen = float(m.group(1))

				yield VepResult(var_id=var_id, chr=chr, start=start, allele=allele,
								gene=gene, transcript=feature, consequences=consequences,
								protein_pos = protein_pos, aa_change=aa_change, protein=protein,
								sift=sift, polyphen=polyphen)

	def close(self):
		"""
		Removes temporary files
		"""

		if self.results_path is not None:
			os.remove(self.results_path)
			self.results_path = None
