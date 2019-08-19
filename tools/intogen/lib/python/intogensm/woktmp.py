import os
import shutil
import tempfile

def make_temp_file(task, suffix=""):
	return tempfile.mkstemp(
		prefix=task.id + ".",
		suffix=suffix,
		dir=task.conf.get("local_temp_path"))[1]

def make_temp_dir(task, suffix=""):
	return tempfile.mkdtemp(
		prefix=task.id + ".",
		suffix=suffix,
		dir=task.conf.get("local_temp_path"))

def remove_temp(task, *files):
	if task.conf.get("local_temp_remove", True, dtype=bool):
		if isinstance(files, basestring):
			files = [files]

		for file in files:
			if os.path.isdir(file):
				task.logger.debug("Removing temporary folder {0} ...".format(file))
				shutil.rmtree(file)
			else:
				task.logger.debug("Removing temporary file {0} ...".format(file))
				os.remove(file)