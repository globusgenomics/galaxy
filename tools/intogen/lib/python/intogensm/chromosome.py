import re

CHR_INDEX = {
	"1" : 1, "2" : 2, "3" : 3, "4" : 4, "5" : 5, "6" : 6, "7" : 7, "8" : 8, "9" : 9, "10" : 10,
	"11" : 11, "12" : 12, "13" : 13, "14" : 14, "15" : 15, "16" : 16, "17" : 17, "18" : 18, "19" : 19, "20" : 20,
	"21" : 21, "22" : 22, "23" : 23, "X" : 24, "Y" : 25
}

STRAND_INDEX = {
	"+" : 0, "+1" : 0, "1" : 0,
	"-" : 1, "-1" : 1
}

_CHR_RE = re.compile(r"^(?:chr)?([XY]|1[0-9]|2[0-2]|[1-9])$", re.IGNORECASE)

def valid_chromosome(chr):
	return _CHR_RE.match(chr) is not None

def parse_chromosome(chr):
	m = _CHR_RE.match(chr)
	if m is not None:
		return m.group(1).upper()
	return None

__CB = {
	"A" : "T",
	"T" : "A",
	"G" : "C",
	"C" : "G"
}

def complementary_sequence(seq):
	return "".join([__CB[base] if base in __CB else base for base in seq.upper()])