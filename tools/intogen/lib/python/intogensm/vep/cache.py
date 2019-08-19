class VepCache(object):
	def __init__(self, path):
		self.path = path

		self.__conn = None

	def open(self):
		pass

import sqlite3

from datetime import datetime

from model import VepResult

DB_VERSION = "1"

DATE_FMT = "%Y-%m-%d %H:%M:%S"

class VepCache(object):
	def __init__(self, path, part_size=None):
		self.path = path
		self.part_size = part_size

		self.__id_cache = {}
		self.__partitions = set()

	def open(self):
		self.__conn = sqlite3.connect(self.path)
		#self.__conn.row_factory = sqlite3.Row

		c = self.__conn.cursor()

		try:
			c.execute("SELECT db_version FROM meta")
			res = c.fetchone()
			db_version = res[0]
		except:
			self.create()
			db_version = DB_VERSION

		if db_version != DB_VERSION:
			self.drop()
			self.create()
			db_version = DB_VERSION

		self.db_version = db_version

		c.execute("SELECT creation_time, part_size FROM meta")
		res = c.fetchone()
		if res is None:
			raise Exception("unexpected error: meta table is empty")

		creation_time, part_size = res

		self.creation_time = datetime.strptime(creation_time, DATE_FMT)

		if self.part_size is not None and part_size != self.part_size:
			raise Exception("The database part_size ({0}) doesn't match with the asked one ({1})".format(part_size, self.part_size))

		self.part_size = part_size

		c.execute("SELECT name FROM partitions")
		self.__partitions = set([r[0] for r in c])

		c.close()

	def create(self):
		c = self.__conn.cursor()

		c.execute("CREATE TABLE meta (db_version, creation_time, part_size INTEGER)")
		c.execute("INSERT INTO meta VALUES (?,?,?)",
			(DB_VERSION, datetime.now().strftime(DATE_FMT), self.part_size))

		c.execute("CREATE TABLE partitions (chr, part INTEGER, name, PRIMARY KEY (chr, part))")

		#c.execute("CREATE TABLE ref_genome_variant (ref_genome_variant_id INTEGER PRIMARY KEY, name UNIQUE)")
		#c.execute("CREATE TABLE gene (gene_id INTEGER PRIMARY KEY, name UNIQUE)")
		#c.execute("CREATE TABLE uniprot (uniprot_id INTEGER PRIMARY KEY, name UNIQUE)")
		#c.execute("CREATE TABLE info (info_id INTEGER PRIMARY KEY, name UNIQUE)")
		#c.execute("CREATE TABLE uniprot_variant (uniprot_variant_id INTEGER PRIMARY KEY, name UNIQUE)")
		#c.execute("CREATE TABLE func_impact (func_impact_id INTEGER PRIMARY KEY, name UNIQUE)")

		c.close()

		self.__conn.commit()

	def drop(self):
		c = self.__conn.cursor()
		c2 = self.__conn.cursor()
		try:
			c.execute("SELECT name FROM partitions")
			for r in c:
				c2.execute("DROP TABLE IF EXISTS {0}".format(r[0]))
		finally:
			c2.close()
		c.execute("DROP TABLE IF EXISTS meta")
		#c.execute("DROP TABLE IF EXISTS uniprot")
		c.execute("DROP TABLE IF EXISTS partitions")
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

	def add(self, chr, start, ref, alt, ref_genome_variant, gene, uniprot, info, uniprot_variant, func_impact, fi_score):

		if uniprot is None and fi_score is None:
			return

		part = start / self.part_size

		part_name = "_{0}_{1}".format(chr, part)

		c = self.__conn.cursor()

		if part_name not in self.__partitions:
			c.execute("INSERT INTO partitions VALUES (?,?,?)", (chr, part, part_name))

			c.execute("""
				CREATE TABLE {0} (
					start INTEGER, ref_alt,
					uniprot_id INTEGER, fi_score REAL,
					PRIMARY KEY (start, ref_alt)
			)""".format(part_name))

			"""
			CREATE TABLE {0} (
					start INTEGER, ref_alt, ref_genome_variant,
					gene_id INTEGER, uniprot_id INTEGER, info_id INTEGER, uniprot_variant,
					func_impact_id INTEGER, fi_score REAL,
					PRIMARY KEY (start, ref_alt)
			"""

			self.__partitions.add(part_name)

		#ref_genome_variant_id = self.__get_id(c, "ref_genome_variant", ref_genome_variant)
		#gene_id = self.__get_id(c, "gene", gene)
		uniprot_id = self.__get_id(c, "uniprot", uniprot)
		#info_id = self.__get_id(c, "info", info)
		#uniprot_variant_id = self.__get_id(c, "uniprot_variant", uniprot_variant)
		#func_impact_id = self.__get_id(c, "func_impact", func_impact)

		"""
		c.execute("INSERT INTO {0} VALUES (?,?,?,?,?,?,?,?,?)".format(part_name),
									(start, "{0}/{1}".format(ref, alt), ref_genome_variant,
			 						gene_id, uniprot_id, info_id, uniprot_variant,
									func_impact_id, fi_score))
		"""
		c.execute("INSERT INTO {0} VALUES (?,?,?,?)".format(part_name),
			(start, "{0}/{1}".format(ref, alt),
			 uniprot_id, fi_score))

		c.close()

	def get(self, chr, start, ref, alt):
		part = start / self.part_size

		part_name = "_{0}_{1}".format(chr, part)

		if part_name not in self.__partitions:
			return None

		c = self.__conn.cursor()

		c.execute("""
			SELECT u.name, p.fi_score FROM {0} p
			LEFT JOIN uniprot u USING (uniprot_id)
			WHERE start=? AND ref_alt=?
		""".format(part_name), (start, "{0}/{1}".format(ref, alt)))

		res = c.fetchone()
		if res is None:
			return None

		uniprot, fi_score = res
		c.close()

		return MaResult(chr=chr, start=start, ref=ref, alt=alt, uniprot=uniprot, fi_score=fi_score)

	def stats(self):
		c = self.__conn.cursor()
		part_size = {}
		chr_parts = {}
		chr_list = []
		last_chr = None

		c.execute("SELECT chr, name FROM partitions ORDER BY chr*1,chr,part")

		for r in c.fetchall():
			chr, part_name = r

			if chr != last_chr:
				last_chr = chr
				chr_list += [chr]

			if chr in chr_parts:
				chr_parts[chr] += [part_name]
			else:
				chr_parts[chr] = [part_name]

			c.execute("SELECT COUNT(*) FROM {0}".format(part_name))
			size = c.fetchone()[0]
			part_size[part_name] = size

		c.close()

		return chr_list, chr_parts, part_size

	def commit(self):
		self.__conn.commit()

	def close(self):
		if self.__conn is not None:
			self.__conn.close()
			self.__conn = None
