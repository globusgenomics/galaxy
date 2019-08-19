def read_header(file_path, delimiter="\t"):
	hdr = {}
	f = open(file_path,"r")
	line = f.readline()
	if line is None:
		return hdr
	columns = line.rstrip().split(delimiter)
	indices = {}
	for i, col_name in enumerate(columns):
		indices[col_name] = i
	hdr["names"] = columns
	hdr["indices"] = indices
	f.close()
	return hdr

def index_file(file_path, key_columns, delimiter="\t"):
	index = {}
	f = open(file_path, "r")
	f.readline() # discard header
	file_pos = f.tell()
	line_num = 1
	line = f.readline()
	while line:
		fields = line.rstrip().split(delimiter)
		key = tuple([fields[i] for i in key_columns])
		if key in index:
			raise Exception("Duplicated key {0} at line {1} of file {2}".format(str(key), line_num, file_path))
		index[key] = (file_pos, line_num)
		file_pos = f.tell()
		line_num += 1
		line = f.readline()
	f.close()
	return index

def join_files(left_file, right_file, out_file, delimiter="\t", fields_filter=None):
	left_header = read_header(left_file, delimiter)
	right_header = read_header(right_file, delimiter)
	left_names = set(left_header["names"])
	right_names = set(right_header["names"])

	key_columns = left_names.intersection(right_names)
	left_indices = left_header["indices"]
	right_indices = right_header["indices"]
	left_key_columns = sorted([left_indices[name] for name in key_columns])
	right_key_columns = sorted([right_indices[name] for name in key_columns])
	right_data_columns = sorted([right_indices[name] for name in right_indices if name not in key_columns])
	right_index = index_file(right_file, right_key_columns, delimiter)

	empty_right_fields = ["\N"] * len(right_data_columns)

	left_f = open(left_file, "r")
	left_f.readline() # discard header
	right_f = open(right_file, "r")
	out_f = open(out_file, "w")
	out_f.write(delimiter.join(left_header["names"]))
	out_f.write(delimiter)
	out_f.write(delimiter.join([right_header["names"][i] for i in right_data_columns]))
	out_f.write("\n")
	for line in left_f:
		fields = line.rstrip().split(delimiter)
		key = tuple([fields[i] for i in left_key_columns])
		if key in right_index:
			file_pos, line_num = right_index[key]
			right_f.seek(file_pos)
			right_fields = right_f.readline().rstrip().split(delimiter)
			right_fields = [right_fields[i] for i in right_data_columns]
		else:
			right_fields = empty_right_fields
		if fields_filter is not None:
			fields_filter(fields, right_fields)
		out_f.write(delimiter.join(fields))
		out_f.write(delimiter)
		out_f.write(delimiter.join(right_fields))
		out_f.write("\n")
	out_f.close()
	right_f.close()
	left_f.close()
