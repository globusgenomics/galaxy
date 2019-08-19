import os
import json
import tempfile

from model import VepException, VepResult

from bgcore import tsv
from bgcore.request import Request

def _ctype(value):
	return value.split(",")

class VepService(object):
	HOST = "beta.rest.ensembl.org"

	VEP_STRAND = { "+" : "1", "-" : "-1", "1" : "1", "-1" : "-1" }

	def __init__(self, cache_path, max_retries=3, max_freq=3):
		self.cache_path = cache_path
		self.results_path = None

		self.__restful = Request(max_retries=max_retries, max_freq=max_freq)

	def __parse_response(self, var_id, chr, start, end, strand, alt, response):
		root = json.load(response)

		if not isinstance(root, dict):
			raise Exception("Unexpected result from VEP web service:\n{0}".format(json.dumps(root)))

		results = []
		found = set()

		tag = ":".join([chr, str(start), str(end), strand, alt])

		for data in root["data"]:
			#chromosome = data["location"]["name"];
			#start = data["location"]["start"];

			for trans in data["transcripts"]:
				gene = trans.get("gene_id");
				transcript = trans.get("transcript_id")

				tstart = trans.get("translation_start")
				tend = trans.get("translation_end")
				if tstart is not None and tend is not None and tstart != tend:
					protein_pos = "{0}-{1}".format(tstart, tend)
				elif tstart is not None:
					protein_pos = tstart
				elif tend is not None:
					protein_pos = tend
				else:
					protein_pos = None

				protein = trans.get("translation_stable_id")

				for allele in trans.get("alleles", []):
					consequences = allele.get("consequence_terms")
					#allele_string = allele["allele_string"]
					aa_change = allele.get("pep_allele_string")
					sift_score = allele.get("sift_score")
					polyphen_score = allele.get("polyphen_score")

					key = "{0}|{1}".format(tag, transcript)

					if key not in found:
						found.add(key)

						results += [VepResult(
										var_id=var_id, chr=chr, start=start, allele=allele,
										gene=gene, transcript=transcript, consequences=consequences,
										protein_pos = protein_pos, aa_change=aa_change, protein=protein,
										sift=sift_score, polyphen=polyphen_score)]

		return results

	def get(self, chr, start, end, strand, alt, var_id=None):
		strand = self.VEP_STRAND[strand]

		url = "http://{0}/vep/human/{1}:{2}-{3}:{4}/{5}/consequences".format(
			self.HOST, chr, start, end, strand, alt)

		response = self.__restful.get(url, headers={"Content-type" : "application/json"})
		if response is None:
			return None

		return self.__parse_response(var_id, chr, start, end, strand, alt, response)

	def run(self, variants_path):
		"""
		Run the VEP service and save results in a temporary file.

		:param variants_path: File with variants. In BED format. http://www.ensembl.org/info/docs/variation/vep/vep_script.html#custom_formats
		:return: True if successfull or False otherwise
		"""

		if self.results_path is None:
			self.results_path = tempfile.mkstemp()[1]

		with open(self.results_path, "w") as rf:
			with open(variants_path, "r") as vf:
				column_types = (str, int, int, str, str, int)
				for fields in tsv.lines(vf, column_types):
					chr, start, end, allele, strand, var_id = fields

					alt = allele[allele.find("/") + 1:]

					results = self.get(chr, start, end, strand, alt, var_id)
					if results is None:
						continue

					for r in results:
						rf.write(tsv.line_text(
							var_id, chr, start, allele,
							r.gene, r.transcript, ",".join(sorted(r.consequences)),
							r.protein_pos, r.aa_change, r.protein,
							r.sift, r.polyphen, null_value="-"))

	def results(self):
		"""
		Iterator that parses the results temporary file and yields VepResult's
		"""

		with open(self.results_path, "r") as f:
			column_types = (int, str, int, str, str, str, _ctype, str, str, str, float, float)
			for fields in tsv.lines(f, column_types, null_value="-"):
				var_id, chr, start, allele,	gene, transcript, consequences, protein_pos, aa_change, protein, sift, polyphen = fields

				yield VepResult(var_id=var_id, chr=chr, start=start, allele=allele,
					gene=gene, transcript=transcript, consequences=consequences,
					protein_pos = protein_pos, aa_change=aa_change, protein=protein,
					sift=sift, polyphen=polyphen)

	def close(self):
		"""
		Removes temporary files
		"""

		if self.results_path is not None:
			os.remove(self.results_path)
			self.results_path = None
