import os
import time

def split_ext(file_name):
	name, ext = os.path.splitext(file_name)
	if ext.lower() in [".gz", ".bz2"]:
		name2, ext2 = os.path.splitext(name)
		if ext2.lower() == "tar":
			name = name2
			ext = ext2 + ext
	return (name, ext)

def ensure_path_exists(path):
	try:
		os.makedirs(path)
	except Exception as ex:
		time.sleep(0.5)
		if not os.path.exists(path):
			raise ex