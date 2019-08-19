import __builtin__

import sys
import os.path
import gzip
import bz2

def split_path(path):
	"""
	Splits a path into a tuple with the dirname, the file name and the extension.
	It is aware of extensions such as tsv.gz, vcf.bz2, tar.gz, ...
	:param path: The path to split
	:return tuple(dirname, name, ext)
	"""

	dirname = os.path.dirname(path)
	basename = os.path.basename(path)
	name, ext = os.path.splitext(basename)
	if ext in [".gz", ".bz2"]:
		name, ext2 = os.path.splitext(name)
		ext = ext2 + ext

	return (dirname, name, ext)

def open(path, mode="r"):
	"""
	Opens a file or the standard input/output if path is '-'.
	It recognises .gz and .bz2 extensions and opens them transparently.
	:param path: File path
	:param mode: r for read and w for write
	:return: a file handle
	"""
	if path == "-" and mode in ["r", "w"]:
		return {
			"r" : sys.stdin,
			"w" : sys.stdout }[mode]

	name, ext = os.path.splitext(path)
	ext = ext.lower()
	if ext == ".gz":
		if mode in ["w", "r"]:
			mode += "b"
		return gzip.GzipFile(path, mode)
	elif ext == ".bz2":
		return bz2.BZ2File(path, mode)
	else:
		return __builtin__.open(path, mode)

def skip_comments_and_empty(f):
	"""
	Skip comments and empty lines
	:param f: the tsv file handler
	:return: the first not empty line with no comments
	"""

	line = f.readline()
	while line.startswith("#") or line == "\n":
		line = f.readline()

	return line

def header_from_line(line):
	"""
	Parse a tsv header
	:param line: the header text
	:return: (column_names list, column_indices dict)
	"""

	if len(line) == 0:
		return [], {}

	column_names = line.rstrip("\n").split("\t")
	column_indices = dict([(name, index) for index, name in enumerate(column_names)])

	return column_names, column_indices

def header(f, comments=True):
	"""
	Parse a tsv header
	:param f: the file handler
	:param comments: skip comments ?
	:return: (column_names list, column_indices dict)
	"""

	if comments:
		line = skip_comments_and_empty(f) # Discard comments
	else:
		line = f.readline()

	return header_from_line(line)

def line(line, types, null_value=None, columns=None):
	"""
	Parses a tsv line
	:param line: The line to parse
	:param types: fields types
	:param null_value: value to interpret as a None
	:param columns: which columns to return
	:return: selected columns
	"""
	assert columns is None or (len(types) == len(columns))

	fields = line.rstrip("\n").split("\t")
	flen = len(fields)
	res = list()
	for i, vtype in enumerate(types):
		if columns is not None:
			index = columns[i]
		else:
			index = i

		if index < flen:
			field = fields[index]
		else:
			field = null_value

		if i < flen and (null_value is None or field != null_value):
			field = vtype(field)
		else:
			field = None

		res += [field]
	return tuple(res)

def lines(f, types, null_value=None, columns=None, comments=True, header=False):
	"""
	Parses a tsv file.
	:param f: the file handler
	:param types: fields types
	:param null_value: value to interpret as a None
	:param columns: which columns to return
	:param comments: skip comments ?
	:param header: has the tsv header ?
	:return: selected columns for each tsv line
	"""

	if comments:
		l = skip_comments_and_empty(f) # Discard comments
	else:
		l = f.readline()

	if len(l) == 0:
		return

	if header:
		hdr = {}
		for index, name in enumerate(l.rstrip("\n").split("\t")):
			hdr[name] = index

		if columns is not None:
			cols = [0] * len(columns)
			for index, name_or_index in enumerate(columns):
				if isinstance(name_or_index, basestring):
					if name_or_index in hdr:
						cols[index] = hdr[name_or_index]
				else:
					cols[index] = name_or_index
			columns = cols
	else:
		if not l.startswith("#"):
			yield line(l, types, null_value, columns)

	for l in f:
		if not l.startswith("#"):
			yield line(l, types, null_value, columns)

def row(line, null_value=None, columns=None):
	"""
	Parses a tsv line
	:param line: The line to parse
	:param null_value: value to interpret as a None
	:param columns: which columns to return
	:return: selected columns
	"""

	fields = line.rstrip("\n").split("\t")
	flen = len(fields)

	return tuple([fields[index] if index < flen else null_value for index in columns or range(flen)])

def rows(f, null_value=None, columns=None, comments=True, header=False):
	"""
	Parses a tsv file.
	:param f: the file handler
	:param null_value: value to interpret as a None
	:param columns: which columns to return
	:param comments: skip comments ?
	:param header: has the tsv header ?
	:return: selected columns for each tsv line
	"""

	if comments:
		line = skip_comments_and_empty(f) # Discard comments
	else:
		line = f.readline()

	if len(line) == 0:
		return

	if header:
		hdr = dict([(name, index) for index, name in enumerate(line.rstrip("\n").split("\t"))])

		if columns is not None:
			cols = [0] * len(columns)
			for index, name_or_index in enumerate(columns):
				if isinstance(name_or_index, basestring):
					if name_or_index in hdr:
						cols[index] = hdr[name_or_index]
				else:
					cols[index] = name_or_index
			columns = cols
	else:
		if not line.startswith("#"):
			yield row(line, null_value, columns)

	for line in f:
		if not line.startswith("#"):
			yield row(line, null_value, columns)

def line_text(*fields, **kargs):
	"""
	Transform a list of values into a tsv line
	:param fields: list of values as parameters
	:param kargs: optional arguments: null_value: value to interpret as None
	:return: tsv line text
	"""
	null_value = kargs.get("null_value", "")
	return "\t".join([str(x) if x is not None else null_value for x in fields]) + "\n"

def write_line(f, *values, **kargs):
	"""
	Write a tsvline into f from a list of values
	:param f: file handler
	:param values: list of values
	:param kargs: optional arguments: null_value: value to interpret as None
	"""
	null_value = kargs.get("null_value", "")
	f.write(line_text(*values, **kargs))

def write_param(f, key, value):
	"""
	:param key: Parameter key
	:param value: Parameter value
	"""
	f.write("## {0}={1}\n".format(key, value))

def params(f):
	"""
	Read header parameters like ## key=value into a dictionary
	:param f: file handler
	:return: {key: value}
	"""

	params = {}
	pos = f.tell()
	line = f.readline()
	while line.startswith("#") or line == "\n":
		if line.startswith("##"):
			line = line[2:].rstrip("\n")
			try:
				idx = line.index("=")
				key = line[0:idx].lstrip()
				value = line[idx + 1:]
				params[key] = value
			except:
				pass
		pos = f.tell()
		line = f.readline()
	f.seek(pos)

	return params