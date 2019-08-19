import os
import fcntl
import logging

def touchopen(filename, *args, **kwargs):
	# Open the file in R/W and create if it doesn't exist. *Don't* pass O_TRUNC
	fd = os.open(filename, os.O_RDWR | os.O_CREAT)
	# Encapsulate the low-level file descriptor in a python file object
	return os.fdopen(fd, *args, **kwargs)

class Onexus(object):
	def __init__(self, projects_list, logger=None):
		self.projects_list = projects_list

		if logger is None:
			logging.basicConfig(
				format = "%(asctime)s %(name)s %(levelname)-5s : %(message)s",
				datefmt = "%Y-%m-%d %H:%M:%S")
			self._log = logging.getLogger("onexus")
		else:
			self._log = logger

	def add_project(self, user_id, project_id, website_path):
		key = "private\\://{0}/{1}".format(user_id, project_id)
		uri = "private\\://{0}/{1}={2},{0}/{1}".format(user_id, project_id, website_path)

		self._log.info("Adding website project {0} ...".format(key))

		self._log.debug("URI: {0}".format(uri))

		self._log.debug("Opening the projects index ...")

		with touchopen(self.projects_list, "r+") as f:

			self._log.debug("Locking ...")

			fcntl.lockf(f, fcntl.LOCK_EX)

			lines = []

			self._log.debug("Reading repositories ...")
			for line in f:
				line = line.rstrip("\n")

				self._log.debug("  {0}".format(line))

				repo_key = None
				pos = line.find("=")
				if pos >= 0:
					repo_key = line[:pos]

				if repo_key is None or not line.startswith(key):
					lines += [line]

			lines += [uri]

			self._log.debug("Writing repositories ...")

			f.seek(0)
			for line in lines:
				f.write(line + "\n")
			f.truncate()

	def remove_project(self, user_id, project_id):
		key = "private\\://{0}/{1}".format(user_id, project_id)

		self._log.info("Removing website project {0} ...".format(key))

		self._log.debug("Opening the projects index ...")

		with touchopen(self.projects_list, "r+") as f:

			self._log.debug("Locking ...")

			fcntl.lockf(f, fcntl.LOCK_EX)

			lines = []

			self._log.debug("Reading repositories ...")
			for line in f:
				line = line.rstrip("\n")

				self._log.debug("  {0}".format(line))

				repo_key = None
				pos = line.find("=")
				if pos >= 0:
					repo_key = line[:pos]

				if repo_key is None or not line.startswith(key):
					lines += [line]

			self._log.debug("Writing repositories ...")

			f.seek(0)
			for line in lines:
				f.write(line + "\n")
			f.truncate()