import os
import os.path
import shutil
import json

from datetime import datetime

from jinja2 import Environment, ChoiceLoader, PackageLoader, FileSystemLoader

from intogensm.version import VERSION
from intogensm.command import RunCommand
from intogensm.projres import ProjectResults

def copyfile(src, dst):
	dirname = os.path.dirname(dst)
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	#stime = os.stat(src).st_mtime if os.path.exists(src) else None
	#dtime = os.stat(dst).st_mtime if os.path.exists(dst) else None
	#updated = stime - dtime > 1 if stime is not None and dtime is not None else True
	#if updated:
	shutil.copy2(src, dst)

def copytree(src, dst):
	if not os.path.exists(dst):
		os.makedirs(dst)
	for item in os.listdir(src):
		s = os.path.join(src, item)
		d = os.path.join(dst, item)
		if os.path.isdir(s):
			copytree(s, d)
		else:
			copyfile(s, d)

def tojson(data):
	return json.dumps(data)

class QcReportCommand(RunCommand):

	cmd_name = "qc-report"

	#def add_arguments(self, parser):
	#	parser.add_argument("-w", "--workspace", dest="workspace", metavar="NAME", default="default",
	#					  help = "Define the workspace name.")

	def execute(self):
		workspace_path = os.path.join(self.results_path, self.args.workspace)
		projects_path = os.path.join(workspace_path, "projects")
		qc_path = os.path.join(workspace_path, "quality_control")

		self.log.info("Preparing output path at {} ...".format(os.path.relpath(qc_path, self.runtime_path)))

		self.create_output_path(qc_path)

		loader = ChoiceLoader([
			PackageLoader(__name__, 'qc_templates'),
			FileSystemLoader(os.path.join(self.root_path, "web", "templates"))
		])

		env = Environment(loader=loader, autoescape=False)
		env.globals.update(
			version=VERSION,
			creation_date=str(datetime.now())
		)

		env.filters["tojson"] = tojson

		self.log.info("Looking for projects ...")

		projects = []
		quality_controls = []

		for project_path in sorted(os.listdir(projects_path)):
			project_path = os.path.join(projects_path, project_path)
			proj_res = ProjectResults(path=project_path)
			project = proj_res.load_def()
			quality = proj_res.load_quality_control()
			projects += [project]
			quality_controls += [quality]

		self.log.info("Processing {} projects ...".format(len(projects)))

		for i, project in enumerate(projects):
			self.log.info("Processing report for project {} ...".format(project["id"]))

			quality = quality_controls[i]

			t = env.get_template("project.html")

			context=dict(
				index=i,
				project=project,
				quality=quality,
				projects=projects)

			t.stream(context).dump(os.path.join(qc_path, "project_{}.html".format(project["id"])))

		self.log.info("Generating index.html ...")

		t = env.get_template("index.html")

		t.stream(projects=projects).dump(os.path.join(qc_path, "index.html"))

	def create_output_path(self, path):
		if not os.path.exists(path):
			os.makedirs(path)

		self.log.debug("Copying static files ...")

		web_static_path = os.path.join(self.root_path, "web", "static")
		static_path = os.path.join(path, "static")

		copytree(os.path.join(os.path.dirname(__file__), "qc_static"), static_path)

		files = [
			"favicon.ico",
			"css/quality.css",
			"js/jqplot",
			"js/quality-plots.js"]

		for file in files:
			src_path = os.path.join(web_static_path, file)
			dst_path = os.path.join(static_path, file)
			if os.path.isdir(src_path):
				copytree(src_path, dst_path)
			elif os.path.isfile(src_path):
				copyfile(src_path, dst_path)
			else:
				raise Exception("File not found: {}".format(src_path))
