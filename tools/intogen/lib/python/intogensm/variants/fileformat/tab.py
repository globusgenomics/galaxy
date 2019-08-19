import re

from intogensm.model import Sample, Variant
from intogensm.chromosome import parse_chromosome

from base import TextParser, ParserException, SkipLine

_ALLELE_RE = re.compile(r"^(\-|\*|(?:\+|\-)?[ACGTN]+|(?:\-|[ACGTN]+)(?:\/(?:\-|[ACGTN]+))?)$", re.IGNORECASE)

class TabParser(TextParser):

	name = "Tabulated text"

	def __init__(self, f, fname, default_sample_id):
		TextParser.__init__(self, f, fname, default_sample_id)

		#f.readline() # Discard header

	def __skip_line(self, line):
		if len(line) == 0:
			return False

		line = line.lstrip()
		return line.startswith("#") or len(line) == 0

	def __read_fields(self):
		line = self._readline()
		while self.__skip_line(line):
			line = self._readline()

		if len(line) == 0:
			raise StopIteration()

		return line.rstrip("\n").split("\t")

	def next(self):
		TextParser.next(self)

		var = None
		while var is None:
			fields = self.__read_fields()

			if len(fields) < 5:
				self._discard_line()
				continue

			if len(fields) < 6:
				fields += [self._default_sample_id]

			chr, start, end, strand, allele, sample = fields[0:6]

			#print ">>>", chr, start, end, strand, allele, sample

			# Chromosome

			chr = parse_chromosome(chr)
			if chr is None:
				self._discard_line()
				continue

			# Start and end

			try:
				start = int(start)
			except:
				self._discard_line()
				continue

			try:
				end = int(end)
			except:
				self._discard_line()
				continue

			if start > end:
				start, end = end, start

			# Strand

			if len(strand) == 0 or strand == "1" or strand == "+1":
				strand = "+"
			elif strand == "-1":
				strand = "-"
			elif strand not in ["+", "-"]:
				self._discard_line()
				continue

			# Alleles

			alleles = allele.split(">")
			if len(alleles) != 2:
				self._discard_line()
				continue

			ref, alt = alleles

			# Check that are well formed
			try:
				for a in [ref, alt]:
					if _ALLELE_RE.match(a) is None:
						self._discard_line()
						raise SkipLine()
			except SkipLine:
				continue

			# Special cases
			if ref == "-" or alt == "-": # ->A, GCT>-
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
			elif ref == "*" and len(alt) > 1 and alt[0] in ["-", "+"]: # *>-ACG, *>+CG
				if alt[0] == "-":
					start -= 1
					ref = "N" + alt[1:]
					alt = "N"
				elif alt[0] == "+":
					ref = "N"
					alt = "N" + alt[1:]
			elif "/" in ref or "/" in alt: # A/A>-/GGT, C/C>A/T, C/C>G/G
				ref = ref.split("/")
				alt = alt.split("/")
				if len(ref) != 2 or len(ref) != len(alt):
					self._discard_line()
					continue

				if ref[0] == ref[1] and alt[0] == alt[1]:
					ref.pop()
					alt.pop()

				for i in range(len(ref)):
					allele = "{0}>{1}".format(ref[i], alt[i])
					self._queue_line("\t".join([chr, str(start), str(end), strand, allele, sample]))

				continue

			ref_len = len(ref)
			alt_len = len(alt)

			vtype = Variant.SUBST if ref_len == alt_len else Variant.INDEL

			# Sample

			if len(sample) == 0:
				sample = self._default_sample_id

			var = Variant(type=vtype, chr=chr, start=start, ref=ref, alt=alt, strand=strand,
						  samples=[Sample(name=sample)])


			#from intogensm.variants.utils import var_to_tab
			#print "***", var
			#print "+++", var_to_tab(var)

		return var
