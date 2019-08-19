class LabelFilter(object):
	"""
	Allows to filter labels according to include and exclude lists
	"""

	def __init__(self, include=None, exclude=None):
		if include is None:
			include = set()

		self.include = include

		if exclude is None:
			exclude = set()

		self.exclude = exclude

	def valid(self, label):
		"""
		Check whether the label passes the filter.
		:param label: a label
		:return: True if label not in exclude and (include list is empty or label is in the include list),
		         False otherwise.
		"""
		return (len(self.include) == 0 or label in self.include) and label not in self.exclude

	@property
	def include_count(self):
		return len(self.include)

	@property
	def exclude_count(self):
		return len(self.exclude)

	def load(self, path, clear=True):
		if clear:
			self.include = set()
			self.exclude = set()

		f = open(path) if isinstance(path, basestring) else path

		for line in f:
			line = line.rstrip()
			if len(line) == 0 or line.startswith("#"):
				continue

			if line.startswith("-"):
				line = line[1:]
				if len(line) > 0:
					self.exclude.add(line)
			elif line.startswith("+"):
				line = line[1:]
				if len(line) > 0:
					self.include.add(line)
			else:
				self.include.add(line)

		if isinstance(path, basestring):
			f.close()