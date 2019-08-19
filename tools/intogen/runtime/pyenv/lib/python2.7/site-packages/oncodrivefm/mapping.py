class MatrixMapping(object):
	"""
	Represents mappings between labeled groups and matrix rows.
	"""

	def __init__(self, matrix, group_row_names):
		"""
		Initialization of the matrix mapping.
		:param matrix: A labeled data matrix
		:param group_row_names: {group_name : [row_names]}
		"""

		self.num_groups = len(group_row_names)
		self.group_names = [None] * self.num_groups
		self.group_name_index = {}
		self.matrix_row_names = matrix.row_names
		self.group_row_indices = {}

		for group_index, (group_name, row_names) in enumerate(group_row_names.items()):
			self.group_names[group_index] = group_name
			self.group_name_index[group_name] = group_index
			row_indices = [matrix.row_name_index[name] for name in row_names if name in matrix.row_name_index]
			self.group_row_indices[group_index] = row_indices