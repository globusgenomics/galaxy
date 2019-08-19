import os
import tempfile
import shutil

from datetime import datetime
from threading import Lock

from wok.core.runstates import READY, PAUSED, WAITING, RUNNING, FINISHED, FAILED

from intogensm.server.db import DT_FORMAT, fetchone_dict

from intogensm.archive import Archive

from intogensm.onexus import Onexus

class Project(object):
	def __init__(self,
			proj_id, user_id, path, temp_path, instance_name, state=None, created_time=None,
			data_file=None, assembly=None):
			
		self.id = proj_id
		self.user_id = user_id
		self.path = path
		self.temp_path = temp_path
		self.instance_name = instance_name
		self.state = state
		self.created_time = created_time or datetime.now()
		self.data_file = data_file
		self.assembly = assembly
		
		self.update_progress()
	
	def __modules_count_by_state(self, instance):
		if instance is None:
			counts = {}
		else:
			counts = instance.modules_count_by_state()
		
		total = 0
		for state in [READY, PAUSED, WAITING, RUNNING, FINISHED, FAILED]:
			if state.title not in counts:
				counts[state.title] = 0
			else:
				total += counts[state.title]
		counts["total"] = total
		
		return counts
		
	def update_progress(self, instance=None):
		counts = self.__modules_count_by_state(instance)
		progress = []

		total = float(counts["total"])
		total_percent = 0
		for label, state in [("success", FINISHED), ("danger", FAILED), ("warning", RUNNING)]:
			count = counts[state.title]
			if total == 0:
				percent = 0
			else:
				percent = round(100.0 * count / total)
			total_percent += percent
			progress += [{
				"label" : label,
				"count" : count,
				"tooltip" : "{0}: {1}".format(state.title, count),
				"percent" : percent }]
		
		self.progress = progress
		
		"""		
		info_tooltip = []
		info_count = 0
		for state in [READY, PAUSED, WAITING]:
			count = counts[state.title]
			info_tooltip += ["{0}: {1}".format(state.title, count)]
			info_count += count

		progress["info"] = {
			"count" : info_count,
			"tooltip" : ", ".join(info_tooltip),
			"percent" : 100 - total_percent }
		"""

	def update_from_instance(self, instance):
		self.has_instance = instance is not None
		
		if self.has_instance:
			self.state = instance.state.title
		else:
			self.state = ""

		self.update_progress(instance)

	@property
	def results_available(self):
		zip_path = os.path.join(self.path, "results.zip")
		return os.path.exists(zip_path)

	@property
	def website_available(self):
		website_path = os.path.join(self.path, "website", "results", "project.tsv")
		return os.path.exists(website_path)

class ProjectsManager(object):

	def __init__(self, onexus_projects=None):
		self.onexus_projects = onexus_projects

	def init_db(self, db, conn):
		c = conn.cursor()
		c.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='projects'")
		res = c.fetchone()
		c.close()
		count, = res
		if count == 0:
			self.__create_tables(conn)
	
	def __create_tables(self, conn):
		c = conn.cursor()
		c.execute("""
		CREATE TABLE projects (
			user_id        TEXT,
			proj_id        TEXT,
			instance_name  TEXT,
			created_time   TEXT,
			path           TEXT,
			temp_path      TEXT,

			data_file      TEXT,
			assembly       TEXT,
			
			PRIMARY KEY (user_id, proj_id) ON CONFLICT REPLACE
		)""")
		c.close()
		
	def __get_projects(self, conn, user_id, proj_id=None):
		c = conn.cursor()
		sql = ["""
			SELECT user_id, proj_id, instance_name, created_time,
					path, temp_path, data_file, assembly
			FROM projects WHERE user_id = ?"""]
		
		if proj_id is None:
			params = (user_id, )
		else:
			sql += ["AND proj_id = ?"]
			params = (user_id, proj_id)
			
		c.execute(" ".join(sql), params)
		
		projects = []
		
		r = fetchone_dict(c)
		while r is not None:
			projects += [Project(
							proj_id=r["proj_id"], user_id=r["user_id"],
							path=r["path"], temp_path=r["temp_path"],
							instance_name=r["instance_name"],
							created_time=datetime.strptime(r["created_time"], DT_FORMAT),
							data_file=r["data_file"], assembly=r["assembly"])]
			
			r = fetchone_dict(c)
		
		return projects
	
	def exists_project(self, conn, user_id, proj_id):
		c = conn.cursor()
		c.execute("SELECT COUNT(*) FROM projects WHERE user_id = ? AND proj_id = ?", (user_id, proj_id))
		res = c.fetchone()
		c.close()
		if res is None:
			return False
		return res[0] > 0

	def get_project(self, conn, user_id, proj_id):
		projects = self.__get_projects(conn, user_id, proj_id)
		if len(projects) == 0:
			return None
		return projects[0]

	def get_projects(self, conn, user_id):
		return self.__get_projects(conn, user_id)

	def get_projects_count(self, conn, user_id):
		c = conn.cursor()
		c.execute("SELECT COUNT(*) FROM projects WHERE user_id = ?", (user_id, ))
		res = c.fetchone()
		c.close()
		if res is None:
			return 0
		return res[0]

	def add_project(self, conn, user_id, project):
		c = conn.cursor()

		c.execute("""
			INSERT INTO projects(user_id, proj_id, instance_name, created_time,
		                 path, temp_path, data_file, assembly)
		    VALUES (?, ?, ?, ?, ?, ?, ?, ?)""", (
				user_id, project.id, project.instance_name, project.created_time.strftime(DT_FORMAT),
				project.path, project.temp_path, project.data_file, project.assembly))
		
		c.close()
		
	def remove_project(self, conn, user_id, proj_id):
		project = self.get_project(conn, user_id, proj_id)
		if project is None:
			return
		
		c = conn.cursor()
		c.execute("DELETE FROM projects WHERE user_id = ? AND proj_id = ?", (user_id, proj_id))
		c.close()

		if self.onexus_projects is not None and os.path.exists(os.path.join(project.path, "website")):
			onexus = Onexus(self.onexus_projects)
			onexus.remove_project(user_id, proj_id)

		if os.path.exists(project.path):
			shutil.rmtree(project.path)

		if os.path.exists(project.temp_path):
			shutil.rmtree(project.temp_path)

	def close(self):
		pass
