import os
import json
import shutil
import subprocess

import wok.logger
from wok.config.data import Data

from intogensm.pathutils import ensure_path_exists

class ProjectResults(object):
	def __init__(self, project=None, path=None):
		assert (path is None and project is not None) or (path is not None and project is None)

		if path is not None:
			self.project = self.load_def(path=path)
			self.project["path"] = path
		elif project is not None:
			self.project = project

	# Definition -------------------------------------------------------------------------------------------------------

	def load_def(self, path=None):
		if path is None:
			path = self.project["path"]

		path = os.path.join(path, "project.conf")
		if not os.path.exists(path):
			return {}

		with open(path) as f:
			project_def = json.load(f)

		return project_def

	def save_def(self, path=None):
		if path is None:
			path = self.project["path"]

		annotations = self.project.get("annotations")
		if annotations is None:
			annotations = {}

		p = {
			"id" : self.project["id"],
			"assembly" : self.project["assembly"],
			"files" : [os.path.relpath(fpath, path) for fpath in self.project["files"]],
			"db" : os.path.relpath(self.project["db"], path),
			"annotations" : annotations
		}

		for key in self.project:
			if key not in ["id", "assembly", "files", "db", "annotations"]:
				p[key] = self.project[key]

		with open(os.path.join(path, "project.conf"), "w") as f:
			json.dump(p, f, indent=4, sort_keys=True)

		temp_path = self.project["temp_path"]
		if os.path.exists(temp_path):
			# for debuging purposes
			with open(os.path.join(temp_path, "project.conf"), "w") as f:
				json.dump(self.project, f, indent=4, sort_keys=True)

	def get_annotations_to_save(self, conf, annotations, names=None, values=None):
		assert (names is None and values is None) or (names is not None and values is not None and len(names) == len(values))

		if names is None:
			names = []
			values = []
		else:
			names = [name for name in names]
			values = [value for value in values]

		ann_keys = conf.get("project.annotations")
		if ann_keys is None:
			ann_keys = []
		elif Data.is_list(ann_keys):
			ann_keys = ann_keys.to_native()
		else:
			log = logger.get_logger(__name__)
			log.warn("Wrong type for 'project.annotations', expecting a list but found:\n{0}".format(repr(ann_keys)))
			ann_keys = []

		for ann_key in ann_keys:
			default = None
			if isinstance(ann_key, basestring):
				key = name = ann_key
			elif isinstance(ann_key, list) and len(ann_key) == 2:
				key = ann_key[0]
				name = ann_key[1]

			value = annotations[key] if key in annotations else default

			names += [name]
			values += [value]

		return names, values

	# Quality control --------------------------------------------------------------------------------------------------

	def __load_qc_data(self, qc_path, qc, key):
		try:
			with open(os.path.join(qc_path, "{}.json".format(key)), "r") as f:
				qc[key] = json.load(f)
		except:
			pass

	def load_quality_control(self):
		quality_control = {}

		qc_path = os.path.join(self.project["path"], "quality_control")
		if not os.path.exists(qc_path):
			return quality_control

		self.__load_qc_data(qc_path, quality_control, "variants")
		self.__load_qc_data(qc_path, quality_control, "oncodrivefm")
		self.__load_qc_data(qc_path, quality_control, "oncodriveclust")

		return quality_control

	def save_quality_control(self, key, data):
		qc_path = os.path.join(self.project["path"], "quality_control")
		ensure_path_exists(qc_path)
		with open(os.path.join(qc_path, "{}.json".format(key)), "w") as f:
			json.dump(data, f, indent=True, sort_keys=True)

	# Clean up ---------------------------------------------------------------------------------------------------------

	def clean(self, conf, compress_db=False):
		if not "temp_path" in self.project:
			return

		temp_path = self.project["temp_path"]
		if conf.get("shared_temp_remove", True, dtype=bool) and os.path.exists(temp_path):
			shutil.rmtree(temp_path)

		if compress_db:
			projdb_path = self.project["db"]
			subprocess.call(" ".join(["gzip", projdb_path]), shell=True)
			self.project["db"] = "{0}.gz".format(projdb_path)
			if os.path.exists(projdb_path):
				os.remove(projdb_path)
