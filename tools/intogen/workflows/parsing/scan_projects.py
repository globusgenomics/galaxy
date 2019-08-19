import os.path
import shutil
import json
import re
import os

from wok.task import task
from wok.config.data import Data

from intogensm.exceptions import InternalError
from intogensm.utils import match_id, list_projects
from intogensm.paths import get_project_path, get_temp_path, create_combination_folders
from intogensm.projres import ProjectResults

def make_project(log, conf, project, base_path=None):
	project_id = project["id"]

	project_path = get_project_path(conf, project_id)
	if not os.path.exists(project_path):
		os.makedirs(project_path)

	temp_path = get_temp_path(conf, project_id)
	if not os.path.exists(temp_path):
		os.makedirs(temp_path)

	projdb_path = os.path.join(project_path, "project.db")

	if "annotations" in project:
		annotations = project["annotations"]
		if not isinstance(annotations, dict):
			log.warn("Overriding project annotations field with an empty dictionary")
			project["annotations"] = annotations = {}
	else:
		project["annotations"] = annotations = {}

	for key in project.keys():
		if key not in ["id", "assembly", "files", "annotations", "oncodriveclust", "oncodrivefm"]:
			value = project[key]
			del project[key]
			annotations[key] = value

	project["path"] = project_path
	project["temp_path"] = temp_path
	project["db"] = projdb_path

	if "assembly" not in project:
		project["assembly"] = "hg19"

	files = project["files"]

	# make absolute paths if necessary
	if base_path is not None:
		for i, file in enumerate(files):
			if not os.path.isabs(file):
				files[i] = os.path.join(base_path, file)

	missing_files = []

	# copy variants files into project path
	for i, file in enumerate(files):
		file_base_name = os.path.basename(file)
		if not os.path.dirname(file).startswith(project_path):
			dst_name = "{0:02d}-{1}".format(i, file_base_name)
			dst_path = os.path.join(project_path, dst_name)
			log.info("  {0} --> {1}".format(file_base_name, dst_name))
			if os.path.exists(file):
				shutil.copy(file, dst_path)
				files[i] = dst_path
			else:
				missing_files += [file]
		else:
			log.info("  {0}".format(file_base_name))

	if len(missing_files) > 0:
		raise Exception("Project {0} references some missing files:\n{1}".format(project_id, "\n".join(missing_files)))

	# save project.conf
	projres = ProjectResults(project)
	projres.save_def()

def single_project(projects_port):
	log = task.logger
	conf = task.conf

	log.info("Single project ...")

	files = conf["files"]
	if isinstance(files, basestring):
		files = [files]
	elif not Data.is_list(files):
		raise InternalError("Unexpected configuration type for 'files', only string and list are allowed")

	files = [str(f) for f in files] #unicode is not json serializable

	project_id = conf["project.id"]
	assembly = conf.get("assembly", "hg19")
	annotations = conf.get("annotations")

	project = {
		"id" : project_id,
		"files" : files,
		"assembly" : assembly
	}

	if annotations is not None and isinstance(annotations, dict):
		project["annotations"] = annotations

	make_project(log, conf, project)

	projects_port.send(project)

def multiple_projects(projects_port):
	log = task.logger
	conf = task.conf

	log.info("Scanning multiple projects ...")

	scan_paths = conf.get("scan_paths")
	if scan_paths is None or not Data.is_list(scan_paths):
		raise InternalError("Unexpected scan_paths type: {0}".format(type(scan_paths)))

	# Prepare includes and excludes
	
	includes = conf.get("includes")
	if includes is not None and not Data.is_list(includes):
		raise InternalError("Unexpected includes type: {0}".format(type(includes)))
	if includes is None or len(includes) == 0:
		includes = ["^.*$"]

	excludes = conf.get("excludes")
	if excludes is not None and not Data.is_list(excludes):
		raise InternalError("Unexpected excludes type: {0}".format(type(excludes)))
	if excludes is None:
		excludes = []

	log.debug("scan_paths:\n{0}".format("\n  ".join(conf["scan_paths"])))
	log.debug("Includes:\n{0}".format("\n  ".join(includes)))
	log.debug("Excludes:\n{0}".format("\n  ".join(excludes)))
	
	# compile regular expressions

	includes = [re.compile(inc) for inc in includes]
	excludes = [re.compile(exc) for exc in excludes]

	# scan paths

	project_ids = set()

	for scan_path in scan_paths:
		for path, project in list_projects(log, scan_path):
			if "id" not in project:
				log.warn("Discarding project that doesn't have 'id': {0}".format(path))
				continue
			if "files" not in project:
				log.warn("Discarding project that doesn't have 'files': {0}".format(path))
				continue

			project_id = project["id"]
			if "name" in project:
				project_name = ": " + project["name"]
			else:
				project_name = ""

			if match_id(project_id, includes) and not match_id(project_id, excludes):
				if project_id in project_ids:
					raise Exception("Duplicated project id at {0}".format(path))

				log.info("Included {0}{1}".format(project_id, project_name))

				make_project(log, conf, project, base_path=os.path.dirname(path))

				projects_port.send(project)
			else:
				log.info("Excluded {0}{1}".format(project_id, project_name))

@task.source()
def scan_projects(projects_port):
	log = task.logger
	conf = task.conf

	log.info("Creating combination folders ...")

	create_combination_folders(conf)

	if "files" in conf:
		single_project(projects_port)
	elif "scan_paths" in conf:
		multiple_projects(projects_port)
	else:
		raise InternalError("Either 'files' or 'scan_paths' configuration parameter is required")

task.run()