from collections import namedtuple

class MaException(Exception):
	pass

MaResultBase = namedtuple('MaResult', 'chr, start, ref, alt, uniprot, fi_score, snps_pos, var_id')

class MaResult(MaResultBase):

	def __new__(cls, chr, start, ref, alt, uniprot, fi_score, snps_pos=None, var_id=None):
		return super(MaResult, cls).__new__(cls, chr, start, ref, alt, uniprot, fi_score, snps_pos, var_id=None)