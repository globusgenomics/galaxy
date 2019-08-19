import re

from intogensm.model import Sample, Variant
from intogensm.chromosome import parse_chromosome

from base import TextParser, ParserException, SkipLine

_ALLELE_RE = re.compile(r"^(\-|[ACGT]+)$", re.IGNORECASE)

_COL_ALLELE1 = "Tumor_Seq_Allele1"
_COL_ALLELE2 = "Tumor_Seq_Allele2"

_COLUMNS = ["Chromosome",
			"Start_Position",
			"Strand",
			"Variant_Type",
			"Reference_Allele",
			_COL_ALLELE1,
			_COL_ALLELE2,
			"Tumor_Sample_Barcode"]

class MafParser(TextParser):

	name = "MAF"

	def __init__(self, f, fname, default_sample_id):
		TextParser.__init__(self, f, fname, default_sample_id)

		self.__format = None

		# Metadata and comments
		line = self._readline()
		while len(line) > 0 and line.startswith("#"):
			if line.startswith("##fileformat="):
				self.__format = line[13:]
			line = self._readline()

		if len(line) > 0: # Header
			column_indices = {}
			columns = line.rstrip().split("\t")
			self._col_size = len(columns)
			for i, name in enumerate(columns):
				if name in _COLUMNS:
					column_indices[name] = i
			try:
				self._col_name_indices = column_indices
				self._col_indices = [column_indices[name] for name in _COLUMNS]
			except KeyError as ex:
				raise ParserException("Header column not found: {0}".format(ex.args[0]), self._location())
		else:
			raise ParserException("Header not found", (self._fname))

	def next(self):
		TextParser.next(self)

		var = None
		while var is None:
			line = self._readline()

			if len(line) == 0:
				raise StopIteration()

			fields = line.rstrip("\n").split("\t")

			chr, start, strand, vtype, ref, alt1, alt2, sample = [
				fields[i] if i < self._col_size else None for i in self._col_indices]

			#print ">>>", chr, start, strand, vtype, ref, alt1, alt2, sample

			# Chromosome

			chr = parse_chromosome(chr)
			if chr is None:
				self._discard_line()
				continue

			# Start

			try:
				start = int(start)
			except:
				self._discard_line()
				continue

			# Strand

			if strand not in ["+"]:
				self._discard_line()
				continue

			# Ref & alt

			if ref is None or alt1 is None or alt2 is None:
				self._discard_line()
				continue

			try:
				for i, x in enumerate([ref, alt1, alt2]):
					if _ALLELE_RE.match(x) is None:
						self._discard_line()
						raise SkipLine()
			except SkipLine:
				continue

			alt = alt1

			if ref == "-":
				# [1   2]  -->  [1 2] 3
				#  . - .         N
				#  . T .         N T
				ref = "N"
				alt = "N" + alt if alt != "-" else "N"
			elif alt == "-":
				# 1 [2] 3  -->  [1 2] 3
				# .  T  .        N T
				# .  -  .        N
				start -= 1
				ref = "N" + ref
				alt = "N"

			ref_len = len(ref)

			vtype = Variant.SUBST if ref_len == len(alt) else Variant.INDEL

			if len(sample) == 0:
				sample = self._default_sample_id

			if alt1 != alt2:
				fields[self._col_name_indices[_COL_ALLELE1]] = fields[self._col_name_indices[_COL_ALLELE2]]
				self._queue_line("\t".join(fields))

			if ref == alt:
				continue

			var = Variant(type=vtype, chr=chr, start=start, ref=ref, alt=alt, strand=strand,
						  samples=[Sample(name=sample)])

			#from intogensm.variants.utils import var_to_tab
			#print "***", var
			#print "+++", var_to_tab(var)
		return var

