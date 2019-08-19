from collections import namedtuple

class VepException(Exception):
	pass

VepResultBase = namedtuple('VepResult', 'var_id, chr, start, allele, gene, transcript, consequences, protein_pos, aa_change, protein, sift, polyphen')

class VepResult(VepResultBase):

	def __new__(cls, var_id, chr, start, allele, gene, transcript, consequences, protein_pos, aa_change, protein=None, sift=None, polyphen=None):
		return super(VepResult, cls).__new__(cls, var_id, chr, start, allele, gene, transcript, consequences, protein_pos, aa_change, protein, sift, polyphen)