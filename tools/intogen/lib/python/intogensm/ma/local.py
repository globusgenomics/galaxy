import tempfile

from db import MaDb
from model import MaResult

class MaLocal(object):

	__CB = {
		"A" : "T",
		"T" : "A",
		"G" : "C",
		"C" : "G"
	}

	def __init__(self, path):
		self.__db = MaDb(path)
		self.__db.open()

	def __complement(self, seq):
		return "".join([self.__CB[base] if base in self.__CB else base for base in seq.upper()])

	def get(self, chr, strand, start, ref, alt, var_id=None):
		if ref == "-" or alt == "-" or len(ref) > 1 or len(alt) > 1:
			return None

		r = self.__db.get(chr, start, ref, alt, var_id=var_id)
		if r is None and strand == "-":
			r = self.__db.get(chr, start, self.__complement(ref), self.__complement(alt), var_id=var_id)

		return r

	def close(self):
		self.__db.close()