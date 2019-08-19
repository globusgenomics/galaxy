import re
import tempfile

from bgcore.request import Request

from intogensm.utils import cast_type

from model import MaResult

class MaService(object):
	HOST = "mutationassessor.org"

	ISSUE_UNKNOWN_ID_TYPE = re.compile(r"unknown ID type")
	ISSUE_REFERENCE_ALLELE = re.compile(r"reference allele: ([ACGT])")

	def __init__(self, assembly, cache_path=None, max_retries=3, max_freq=3):
		self.assembly = assembly
		self.cache_path = cache_path

		self.__restful = Request(max_retries=max_retries, max_freq=max_freq)

	def get(self, chr, strand, start, ref, alt, var_id=None):
		done = False
		while not done:
			response = self.__restful.get("http://{0}/".format(self.HOST),
				params={
					"cm" : "var",
					"var" : "{0},{1},{2},{3},{4}".format(self.assembly, chr, start, ref, alt),
					"frm" : "txt",
					"fts" : "all"
				})

			if response is None:
				return None

			hdr = response.readline().rstrip("\n").split("\t")
			fields = response.readline().rstrip("\n").split("\t")
			hlen = len(hdr)
			if hlen == 0 or hlen != len(fields):
				return None

			r = {}
			for i in range(hlen):
				r[hdr[i]] = fields[i] if len(fields[i]) > 0 else None

			mapping_issue = r["Mapping issue"]
			if mapping_issue is not None:
				if self.ISSUE_UNKNOWN_ID_TYPE.match(mapping_issue):
					return None

				m = self.ISSUE_REFERENCE_ALLELE.match(mapping_issue)
				if m is not None:
					ref = m.group(1)
					strand = {"+" : "-", "-" : "+"}[strand]
					continue

					#TODO check: raise Exception("Infinite mapping issue for reference allele")

			done = True

		uniprot=r["Uniprot"]
		fi_score=cast_type(r["FI score"], float)
		snps_pos=r["SNPs@position"]

		if uniprot is not None or fi_score is not None or snps_pos is not None:
			return MaResult(
				var_id=var_id, chr=chr, start=start, ref=ref, alt=alt,
				uniprot=uniprot, fi_score=fi_score, snps_pos=snps_pos)
		else:
			return None

	def close(self):
		pass