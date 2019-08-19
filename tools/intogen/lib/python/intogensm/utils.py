import os
import re
import json

from wok.config.data import Data

def normalize_id(name):
	return re.sub(r"[^a-zA-Z0-9-]", "_", name)

def cast_type(value, type):
	if value is None:
		return None
	return type(value)

def get_project_conf(conf, project, key, default=None, dtype=None):
	value = conf.get(key, default=default, dtype=dtype)
	if not Data.is_element(project):
		project = Data.create(project)
	return project.get(key, default=value, dtype=dtype)

# Scanning of projects ----------------------------------------------------------

PROJ_DEF_RE = re.compile(r"^(project.conf|.*\.smproject)$")

def match_id(value, patterns):
	for pattern in patterns:
		if pattern.match(value):
			return True
	return False

def list_projects(log, path, base_path=None):
	log.debug("  {0}".format(os.path.relpath(path, base_path or path)))

	if os.path.isfile(path):
		if PROJ_DEF_RE.match(os.path.basename(path)):
			if base_path is None:
				base_path = os.path.dirname(path)

			try:
				with open(path, "r") as f:
					project = json.load(f)
					if isinstance(project, dict):
						yield (path, project)
					else:
						log.warn("Discarding malformed file {0}".format(os.path.relpath(path, base_path)))
			except:
				log.warn("Unexpected error reading file {0}".format(os.path.relpath(path, base_path)))

	elif os.path.isdir(path):
		if base_path is None:
			base_path = path
		for entry in os.listdir(path):
			entry_path = os.path.join(path, entry)
			for p in list_projects(log, entry_path, base_path=base_path):
				yield p
