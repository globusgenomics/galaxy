import sqlite3

from datetime import datetime
from chromosome import CHR_INDEX

DB_VERSION = "2"

DATE_FMT = "%Y-%m-%d %H:%M:%S"

class SigDb(object):
	def __init__(self, path):
		self.path = path
		self.__conn = None

	def open(self):
		self.__conn = sqlite3.connect(self.path)
		#self.__conn.row_factory = sqlite3.Row

		c = self.__conn.cursor()

		try:
			c.execute("SELECT db_version FROM meta")
			res = c.fetchone()
			db_version, = res
		except:
			self.create()
			db_version = DB_VERSION

		if db_version != DB_VERSION:
			self.drop()
			self.create()
			db_version = DB_VERSION

		self.db_version = db_version

		c.execute("SELECT creation_time FROM meta")
		res = c.fetchone()
		if res is None:
			raise Exception("unexpected error: meta table is empty")

		creation_time, = res

		self.creation_time = datetime.strptime(creation_time, DATE_FMT)

		c.close()

	def create(self):
		c = self.__conn.cursor()

		c.execute("CREATE TABLE meta (db_version, creation_time INTEGER)")
		c.execute("INSERT INTO meta VALUES (?,?)",
				  (DB_VERSION, datetime.now().strftime(DATE_FMT)))

		c.executescript("""
			CREATE TABLE variants (
				pos INTEGER PRIMARY KEY
			);

			CREATE TABLE genes (
				gene        TEXT PRIMARY KEY,
				fm          INTEGER,
				clust       INTEGER
			);
		""")

		c.close()

		self.__conn.commit()

	def drop(self):
		c = self.__conn.cursor()
		c.execute("DROP TABLE IF EXISTS meta")
		c.execute("DROP TABLE IF EXISTS variants")
		c.execute("DROP TABLE IF EXISTS genes")
		c.close()
		self.__conn.commit()

	def __get_pos(self, chr, start):
		return CHR_INDEX[chr] << 33 | start

	def add_variant(self, chr, start):

		try:
			pos = self.__get_pos(chr, start)
		except:
			return

		c = self.__conn.cursor()

		try:
			c.execute("INSERT INTO variants (pos) VALUES (?)", (pos, ))
		except sqlite3.IntegrityError:
			pass

		c.close()

	def exists_variant(self, chr, start):
		try:
			pos = self.__get_pos(chr, start)
		except:
			return

		c = self.__conn.cursor()

		exists = False
		try:
			c.execute("SELECT count(*) FROM variants WHERE pos=?", (pos, ))
			exists = c.fetchone()[0] != 0
		except sqlite3.IntegrityError:
			pass

		c.close()

		return exists

	def add_gene(self, gene, fm=False, clust=False):
		c = self.__conn.cursor()

		try:
			c.execute("INSERT INTO genes (gene, fm, clust) VALUES (?,?,?)",
					  (gene, 1 if fm else 0, 1 if clust else 0))
		except sqlite3.IntegrityError:
			c.execute("UPDATE genes SET fm=?, clust=? WHERE gene=?",
					  (1 if fm else 0, 1 if clust else 0, gene))

		c.close()

	def exists_gene(self, gene):
		c = self.__conn.cursor()

		exists = False
		try:
			c.execute("SELECT count(*) FROM genes WHERE gene=?", (gene, ))
			exists = c.fetchone()[0] != 0
		except sqlite3.IntegrityError:
			pass

		c.close()

		return exists

	def genes(self):
		for row in self.__conn.execute("SELECT gene, fm, clust FROM genes"):
			yield dict(gene=row["gene"], fm=row["fm"] != 0, clust=row["clust"] != 0)

	def commit(self):
		self.__conn.commit()

	def rollback(self):
		self.__conn.rollback()

	def close(self):
		if self.__conn is not None:
			self.__conn.close()
			self.__conn = None
