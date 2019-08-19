import re

GENOMIC = "G"

CHR_RE = re.compile(r"^(?:chr)?([XY]|1[0-9]|2[0-2]|[1-9])$", re.IGNORECASE)
INTEGER_RE = re.compile(r"^[0-9]+$")
NUCLEOTID_RE = re.compile(r"^[ACGT]$", re.IGNORECASE)
ALLELE_RE = re.compile(r"^([ACGT])[>/]([ACGT])$", re.IGNORECASE)
STRAND_RE = re.compile(r"^(\+|\+1|\-|\-1|1)$")
IDENTIFIER_RE = re.compile(r"^[a-zA-Z0-9_]+$")

def init(sm):
	pass

def chromosome(sm):
	sm.chr = sm.token_groups[0]

def start(sm):
	sm.start = int(sm.token)

def end(sm):
	sm.end = int(sm.token)

def allele(sm):
	sm.ref, sm.alt = sm.token_groups

def ref(sm):
	sm.alt = sm.token

def alt(sm):
	sm.ref = sm.alt
	sm.alt = sm.token

def strand(sm):
	if sm.strand is not None:
		raise UnexpectedToken(self.token)
	sm.strand = sm.token

def identifier(sm):
	sm.identifier = sm.token


class PrematureEnd(Exception):
	pass

class UnexpectedToken(Exception):
	pass

class ParserStateMachine(object):

	initial_state = None
	end_states = []
	transitions = {}
	handlers = {}

	def __init__(self):
		self.__reset()
		
	def __reset(self, line=None):
		self.line = line
		self._index = 0
		self.token = None
		self.token_groups = None
		self._state = self.initial_state
		self._init()

	def _init(self):
		pass

	def _next_token(self):
		if self._index > len(self.line):
			return None

		try:
			pos = self.line.index("\t", self._index)
		except ValueError:
			pos = len(self.line)

		token = self.line[self._index:pos]
		self._index = pos + 1
		return token

	def parse(self, line):
		self.__reset(line)

		done = False
		while not done:
			# action
			self.handlers[self._state](self)

			# transition
			next_state = None
			self.token = self._next_token()
			if self.token is not None:
				for condition, state in self.transitions[self._state]:
					match = condition.match(self.token)
					if match:
						self.token_groups = match.groups()
						next_state = state
						break

			#print "{} ---[{}]--> {} | {}".format(self._state, self.token, next_state, repr(self))

			if next_state is None:
				done = True
				if self._state not in self.end_states:
					if self._index > len(self.line):
						raise PrematureEnd()
					else:
						raise UnexpectedToken(self.token)
			else:
				self._state = next_state

		token = self._next_token()
		while token is not None:
			self.remaining_fields += [token]
			token = self._next_token()


class MutationParser(ParserStateMachine):

	def _init(self):
		self.chr = None
		self.start = None
		self.end = None
		self.ref = None
		self.alt = None
		self.strand = None
		self.identifier = None
		self.remaining_fields = []

	def __repr__(self):
		sb = []
		if self.coord is None:
			sb += ["?"]
		else:
			sb += ["chr:{}".format(self.chr)]
		if self.start is not None:
			sb += ["start:{}".format(self.start)]
		if self.end is not None:
			sb += ["end:{}".format(self.end)]
		if self.ref is not None:
			sb += ["ref:{}".format(self.ref)]
		if self.alt is not None:
			sb += ["alt:{}".format(self.alt)]
		if self.strand is not None:
			sb += ["strand:{}".format(self.strand)]
		if self.identifier is not None:
			sb += ["id:{}".format(self.identifier)]
		if len(self.remaining_fields) > 0:
			sb += ["extra:{}".format(";".join(self.remaining_fields))]
		return " ".join(sb)

class DnaMutationParser(MutationParser):

	initial_state = "INIT"

	end_states = ["ALLELE", "ALT", "STRAND2", "IDENTIFIER"]

	transitions = {
		"INIT" : [
			(CHR_RE, "CHR")
		],
		"CHR" : [
			(INTEGER_RE, "START")
		],
		"START" : [
			(INTEGER_RE, "END"),
			(ALLELE_RE, "ALLELE"),
			(NUCLEOTID_RE, "REF"),
			(STRAND_RE, "STRAND1")
		],
		"END" : [
			(ALLELE_RE, "ALLELE"),
			(NUCLEOTID_RE, "REF")
		],
		"STRAND1" : [
			(ALLELE_RE, "ALLELE"),
			(NUCLEOTID_RE, "REF")
		],
		"ALLELE" : [
			(STRAND_RE, "STRAND2"),
			(IDENTIFIER_RE, "IDENTIFIER")
		],
		"REF" : [
			(NUCLEOTID_RE, "ALT")
		],
		"ALT" : [
			(STRAND_RE, "STRAND2"),
			(IDENTIFIER_RE, "IDENTIFIER")
		],
		"STRAND2" : [
			(IDENTIFIER_RE, "IDENTIFIER")
		],
		"IDENTIFIER" : []
	}

	handlers = {
		"INIT" : init,
		"CHR" : chromosome,
		"START" : start,
		"END" : end,
		"STRAND1" : strand,
		"ALLELE" : allele,
		"REF" : ref,
		"ALT" : alt,
		"STRAND2" : strand,
		"IDENTIFIER" : identifier
	}
