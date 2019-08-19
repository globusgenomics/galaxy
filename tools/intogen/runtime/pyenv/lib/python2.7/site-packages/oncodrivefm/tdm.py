from matrix import Matrix

def load_matrix(path):
	# First step is to know which rows/columns/slices there are.
	with open(path) as f:
		first_data_line=0
		line = f.readline()
		while len(line) > 0 and line.startswith("#"):
			first_data_line += 1
			line = f.readline()

		# If the file is empty
		if len(line) == 0:
			return Matrix(0, 0, 0)

		# Parse header
		first_data_line += 1
		hdr = line.rstrip("\n").split("\t")
		if len(hdr) < 3:
			raise Exception("Malformed header, at least 3 columns are required")

		slice_names = hdr[2:]

		row_indices = {}
		col_indices = {}

		for line in f:
			if line.startswith("#"):
				continue

			r = line.rstrip("\n").split("\t")
			if len(r) < 3:
				raise Exception("Malformed row, at least 3 columns are required: {0}".format(line.rstrip("\n")))

			col_name, row_name = r[0:2]

			if row_name not in row_indices:
				row_indices[row_name] = len(row_indices)
			if col_name not in col_indices:
				col_indices[col_name] = len(col_indices)

	row_names = [None] * len(row_indices)
	for name, index in row_indices.items():
		row_names[index] = name

	col_names = [None] * len(col_indices)
	for name, index in col_indices.items():
		col_names[index] = name

	# Then initialize the Matrix and load the data
	mat = Matrix(num_rows=len(row_names), num_cols=len(col_names), num_slices=len(slice_names),
		row_names=row_names, col_names=col_names, slice_names=slice_names)

	with open(path) as f:
		# Discard headers
		while first_data_line > 0:
			f.readline()
			first_data_line -= 1

		for line in f:
			r = line.rstrip("\n").split("\t")
			col_name, row_name = r[0:2]
			fields = r[2:]
			d = [None] * len(fields)
			for i, x in enumerate(fields):
				try:
					d[i] = float(x)
				except:
					pass

				mat.set_slices(row_name, col_name, d)

	return mat
