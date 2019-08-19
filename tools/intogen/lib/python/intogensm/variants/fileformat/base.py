from collections import deque

class ParserException(Exception):
	def __init__(self, msg, loc=None):

		if loc is not None and len(loc) > 0:
			loc_len = len(loc)
			if loc_len == 1:
				msg = "{0}: {2}".format(loc[0], msg)
			elif loc_len == 2:
				msg = "{0} at line {1}: {2}".format(loc[0], loc[1], msg)
			elif loc_len == 3:
				msg = "{0} at line {1}, column {2}: {3}".format(loc[0], loc[1], loc[2], msg)

		Exception.__init__(self, msg)

class SkipLine(Exception):
	pass

class Parser(object):

	name = ""

	def __init__(self, f, fname, default_sample_id):
		self._f = f
		self._fname = fname
		self._default_sample_id = default_sample_id

	def __iter__(self):
		return self

	def next(self):
		pass

	def discarded_lines(self):
		return []

class TextParser(Parser):

	name = "Text"

	def __init__(self, f, fname, default_sample_id):
		Parser.__init__(self, f, fname, default_sample_id)

		self._line_num = 0
		self._column = 0
		self.__read_lines = []
		self.__discarded_lines = []

		self.__queued_lines = deque()

	def next(self):
		self._column = 0
		self.__read_lines = []
		self.__discarded_lines = []

	def read_lines(self):
		return self.__read_lines

	def discarded_lines(self):
		return self.__discarded_lines

	def _discard_line(self):
		self.__discarded_lines += [self._line_num]

	def get_line_num(self):
		return self._line_num

	def _location(self, column=None):
		if column is not None:
			return (self._fname, self._line_num, column)
		else:
			return (self._fname, self._line_num)

	def _queue_line(self, line):
		self.__queued_lines.append(line + "\n")

	def _readline(self):
		if len(self.__queued_lines) > 0:
			return self.__queued_lines.popleft()

		self._line_num += 1
		self._line_text = self._f.readline()
		self.__read_lines += [(self._line_num, self._line_text.rstrip("\n"))]

		return self._line_text