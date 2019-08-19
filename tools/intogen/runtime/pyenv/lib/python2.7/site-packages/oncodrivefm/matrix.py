import random
import numpy as np

class Matrix(object):
	def __init__(self, num_rows, num_cols, num_slices,
				 initial_value=np.nan, dtype=float,
				 row_names=None, col_names=None, slice_names=None):

		self.num_rows = num_rows
		self.num_cols = num_cols
		self.num_slices = num_slices

		self.data = np.empty((num_slices, num_rows, num_cols), dtype=dtype)
		self.data[:] = initial_value

		self.set_row_names(row_names)
		self.set_col_names(col_names)
		self.set_slice_names(slice_names)

	def set_row_names(self, names):
		self.row_name_index = {}
		self.row_names = names
		if names is not None:
			for i, name in enumerate(names):
				self.row_name_index[name] = i

	def set_col_names(self, names):
		self.col_name_index = {}
		self.col_names = names
		if names is not None:
			for i, name in enumerate(names):
				self.col_name_index[name] = i

	def set_slice_names(self, names):
		self.slice_name_index = {}
		self.slice_names = names
		if names is not None:
			for i, name in enumerate(names):
				self.slice_name_index[name] = i

	def _indices_all(self, row, col, slice):
		if isinstance(row, basestring):
			row = self.row_name_index[row]
		if isinstance(col, basestring):
			col = self.col_name_index[col]
		if isinstance(slice, basestring):
			slice = self.slice_name_index[slice]
		return (row, col, slice)

	def _indices_row_col(self, row, col):
		if isinstance(row, basestring):
			row = self.row_name_index[row]
		if isinstance(col, basestring):
			col = self.col_name_index[col]
		return (row, col)

	def _slice_index(self, slice):
		if isinstance(slice, basestring):
			slice = self.slice_name_index[slice]
		return slice

	def get(self, row, col, slice=0):
		row, col, slice = self._indices_all(row, col, slice)
		return self.data[slice, row, col]

	def set(self, row, col, slice, value):
		row, col, slice = self._indices_all(row, col, slice)
		self.data[slice, row, col] = value

	def set_slices(self, row, col, values):
		row, col = self._indices_row_col(row, col)
		for slice, value in enumerate(values):
			self.data[slice, row, col] = value if value is not None else np.nan

	def masked_invalid(self, slice=0):
		return np.ma.masked_invalid(self.data[self._slice_index(slice)])

	def sample(self, size, slice=0):
		flat_data = self.data[self._slice_index(slice)].flat
		flat_data = flat_data[np.logical_not(np.isnan(flat_data))]
		return random.sample(flat_data, size)

	def __repr__(self):
		sb = ["Matrix("]
		sb += ["num_rows=", str(self.num_rows), ", "]
		sb += ["num_cols=", str(self.num_cols), ", "]
		sb += ["num_slices=", str(self.num_slices)]
		sb += ["):\n"]
		for slice, name in enumerate(self.slice_names):
			sb += ["\nslice ", name, ":\n"]
			sb += ["\t", "\t".join(self.col_names), "\n"]
			for row in xrange(min(self.num_rows, 100)):
				sb += [self.row_names[row]]
				for col in xrange(min(self.num_cols, 100)):
					sb += ["\t", str(self.data[slice, row, col])]
				sb += ["\n"]
		return "".join(sb)

if __name__ == "__main__":
	import sys
	from . import tdm
	print tdm.load_matrix(sys.argv[1])