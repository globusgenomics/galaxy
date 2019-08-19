import os
import sqlite3

from datetime import datetime

from intogensm import VERSION as SOFT_VERSION

DB_VERSION="2"

DT_FORMAT = "%Y-%m-%d %H:%M:%S"

class DbVersionError(Exception):
	pass
	
def fetchone_dict(cursor):
	res = cursor.fetchone()
	if res is None:
		return None
	
	r = {}
	for i, v in enumerate(res):
		r[cursor.description[i][0]] = v
	
	return r

class DbManager(object):
	def __init__(self, path):
		self.path = path
		if not os.path.exists(path):
			os.makedirs(path)
		
		self.db_file = os.path.join(path, "web.db")
		if not os.path.exists(self.db_file):
			self.__create_db()
		
		db_version = self.version(self.open_conn())
		if db_version != DB_VERSION:
			sb = []
			for i in range(int(db_version), int(DB_VERSION)):
				sb += ["{0:03}.sql".format(i)]
			raise DbVersionError("Database schema is too old, required {0} but found {1}. Please run the following SQL files under bin/launcher/db : {2}".format(DB_VERSION, db_version, ",".join(sb)))
	
	def __create_db(self):
		conn = self.open_conn()
		c = conn.cursor()
		
		c.execute("""
		CREATE TABLE db_meta (
		    db_version     TEXT,
		    soft_version   TEXT,
		    created_time   TEXT
		)""")
		
		c.execute("INSERT INTO db_meta VALUES (?, ?, ?)",
					(DB_VERSION, SOFT_VERSION, datetime.now().strftime(DT_FORMAT)))
		
		conn.commit()		
		c.close()
		conn.close()
	
	def open_conn(self):
		return sqlite3.connect(self.db_file)
	
	def version(self, conn):
		c = conn.cursor()
		c.execute("SELECT db_version FROM db_meta")
		res = c.fetchone()
		c.close()
		if res is None:
			return "0"
		return res[0]

	def close(self):
		pass

