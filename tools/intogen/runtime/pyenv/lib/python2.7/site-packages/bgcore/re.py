re = __import__("re")

_PAT_TYPE = type(re.compile(""))

class ReContext(object):
	"""
	Allows to perform regex operations and save the resulting state (MatchObject) for further manipulation.
	Example::
		matcher = ReState()
		if matcher.match("^(.)B$", "B"):
			print matcher.group(1)
		elif matcher.match("^(X|Y)Z$", "XZ"):
			print matcher.group(1), matcher.start(1)
	"""

	def __init__(self, text=None):
		self.text = text
		self.result = None

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_val, exc_tb):
		pass

	def match(self, pattern, text=None, flags=0):
		text = text or self.text
		if text is None:
			raise ValueError("Missing the text to match. Use either the ReContext constructor or the method argument")

		if isinstance(pattern, basestring):
			self.result = re.match(pattern, text, flags)
		elif isinstance(pattern, _PAT_TYPE):
			if flags != 0:
				raise TypeError("When pattern is a compiled regex flags are not used")
			self.result = pattern.match(text)
		else:
			raise TypeError("Expecting a string or a compiled regex for pattern but found: {0}".format(type(pattern)))
		return self.result

	def search(self, pattern, text=None, flags=0):
		text = text or self.text
		if text is None:
			raise ValueError("Missing the text for search. Use either the ReContext constructor or the method argument")

		if isinstance(pattern, basestring):
			self.result = re.search(pattern, text, flags)
		elif isinstance(pattern, _PAT_TYPE):
			if flags != 0:
				raise TypeError("When pattern is a compiled regex flags are not used")
			self.result = pattern.search(text)
		else:
			raise TypeError("Expecting a string or a compiled regex for pattern but found: {0}".format(type(pattern)))
		return self.result

	def expand(self, template):
		return self.result.expand(template)

	def group(self, index):
		return self.result.group(index)

	def groups(self):
		return self.result.groups()

	def groupdict(self, default=None):
		return self.result.groupsdict(default)

	def start(self, group=0):
		return self.result.start(group)

	def end(self, group=0):
		return self.result.end(group)

	def span(self, group=0):
		return self.result.span(group)
