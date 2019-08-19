import sqlite3

from datetime import datetime
from chromosome import CHR_INDEX, STRAND_INDEX, complementary_sequence

DB_VERSION = "4"

DATE_FMT = "%Y-%m-%d %H:%M:%S"

class VarXrefsDb(object):
	def __init__(self, path):
		self.path = path

		self.__id_cache = {}

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

		source_cache = {}
		for source_id, name in c.execute("SELECT source_id, name FROM source"):
			source_cache[name] = source_id
		self.__id_cache["source"] = source_cache

		c.close()

	def create(self):
		c = self.__conn.cursor()

		c.execute("CREATE TABLE meta (db_version, creation_time INTEGER)")
		c.execute("INSERT INTO meta VALUES (?,?)",
			(DB_VERSION, datetime.now().strftime(DATE_FMT)))

		c.executescript("""
			CREATE TABLE IF NOT EXISTS source (
				source_id INTEGER PRIMARY KEY,
				name UNIQUE
			);

			CREATE TABLE xrefs (
				pos INTEGER, ref_alt TEXT,
				source_id INTEGER, xref TEXT
			);

			CREATE INDEX xrefs_pk ON xrefs (pos, ref_alt);

			PRAGMA synchronous=OFF;
		""")

		c.close()

		self.__conn.commit()

	def drop(self):
		c = self.__conn.cursor()
		c.execute("DROP TABLE IF EXISTS meta")
		c.execute("DROP TABLE IF EXISTS source")
		c.execute("DROP TABLE IF EXISTS xrefs")
		c.close()
		self.__conn.commit()

	def __get_id(self, c, var_name, value):
		if value is None:
			return None

		if var_name in self.__id_cache:
			var_cache = self.__id_cache[var_name]
			if value in var_cache:
				return var_cache[value]
		else:
			self.__id_cache[var_name] = {}

		c.execute("INSERT INTO {0} (name) VALUES (?)".format(var_name), (value,))
		var_id = c.lastrowid

		self.__id_cache[var_name][value] = var_id

		return var_id

	def __get_pos(self, chr_index, strand_index, start):
		return chr_index << 33 | strand_index << 32 | start

	def add(self, chr, start, ref, alt, source, xref, strand="+"):

		try:
			pos = self.__get_pos(CHR_INDEX[chr], STRAND_INDEX[strand], start)
			ref_alt = "{0}|{1}".format(ref, alt)
		except:
			return

		c = self.__conn.cursor()

		source_id = self.__get_id(c, "source", source)

		try:
			c.execute("INSERT INTO xrefs (pos, ref_alt, source_id, xref) VALUES (?,?,?,?)",
											(pos, ref_alt, source_id, xref))
		except sqlite3.IntegrityError:
			pass

		c.close()

	def __get_xrefs(self, chr_index, strand_index, start, ref, alt):

		pos = self.__get_pos(chr_index, strand_index, start)
		ref_alt = "{0}|{1}".format(ref, alt)

		c = self.__conn.cursor()

		c.execute("""
			SELECT s.name, x.xref
			FROM xrefs x
			LEFT JOIN source s USING (source_id)
			WHERE x.pos=? AND x.ref_alt=?
		""", (pos, ref_alt))

		xrefs = c.fetchall()

		c.close()

		return xrefs

	def get_xrefs(self, chr, start, ref, alt, strand="+"):
		chr = str(chr)
		strand = str(strand)
		if chr not in CHR_INDEX or strand not in STRAND_INDEX:
			return []

		chr_index = CHR_INDEX[chr]
		strand_index = STRAND_INDEX[strand]
		xrefs = self.__get_xrefs(chr_index, strand_index, start, ref, alt)
		if strand_index == 1 and len(xrefs) == 0:
			ref = complementary_sequence(ref)
			alt = complementary_sequence(alt)
			xrefs = self.__get_xrefs(chr_index, 0, start, ref, alt)

		return xrefs

	def begin(self):
		self.__conn.execute("BEGIN IMMEDIATE TRANSACTION")

	def commit(self):
		self.__conn.commit()

	def rollback(self):
		self.__conn.rollback()

	def close(self):
		if self.__conn is not None:
			self.__conn.close()
			self.__conn = None
