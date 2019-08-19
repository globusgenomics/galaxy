import sys
import sqlite3

from intogensm.model import Sample, Variant, Consequence, AffectedGene, AffectedGeneRec, Gene, GeneRec, Pathway, PathwayRec

class ProjectDb(object):

	# consequences filters
	CSQ_CTYPES = "ctypes"
	CSQ_VAR = "var"
	CSQ_GENE = "gene"

	# oncodrive gene excludes causes
	GENE_EXC_NOT_FOUND = "N"
	GENE_EXC_FILTER = "F"
	GENE_EXC_THRESHOLD = "T"
	NO_GENE_EXC = ""

	def __init__(self, path):
		self.path = path
		self.__conn = sqlite3.connect(path, timeout=2**21-1)
		self.__conn.row_factory = sqlite3.Row

	def create(self):
		c = self.__conn.cursor()

		# sources

		c.executescript("""
		CREATE TABLE sources (
			source_id      INTEGER PRIMARY KEY,
			name           TEXT UNIQUE
		);

		CREATE TABLE source_lines (
			source_id       INTEGER,
			line_num        INTEGER,
			text            TEXT,

			PRIMARY KEY (source_id, line_num)
		);

		CREATE TABLE source_variants (
			source_id      INTEGER,
			line_num       INTEGER,
			var_id         INTEGER,

			PRIMARY KEY (source_id, line_num, var_id)
		);
		""")

		# samples

		c.execute("""
		CREATE TABLE samples (
			sample_id      INTEGER PRIMARY KEY,
			name           TEXT UNIQUE
		)""")

		# variants

		c.executescript("""
		CREATE TABLE variants (
			var_id         INTEGER PRIMARY KEY,
			type           TEXT,
			chr            TEXT,
			strand         TEXT,
			prev_start     INTEGER,
			start          INTEGER,
			ref            TEXT,
			alt            TEXT,

			UNIQUE (type, chr, strand, start, ref, alt)
		);

		CREATE INDEX variant_position ON variants (chr, start, ref, alt);

		CREATE TABLE variant_samples (
			var_id         INTEGER,
			sample_id      INTEGER,

			PRIMARY KEY (var_id, sample_id)
		);

		CREATE TABLE variant_xrefs (
			var_id			INTEGER,
			xref			TEXT
		);
		""")
		
		# affected_genes
		
		c.executescript("""
		CREATE TABLE affected_genes (
			afg_id			INTEGER PRIMARY KEY,
			var_id			INTEGER,
			gene_id			TEXT,
			impact			INTEGER,
			coding_region	INTEGER,
			prot_changes	TEXT,

			UNIQUE (var_id, gene_id)
		);
		""")
		
		# consequences

		c.executescript("""
		CREATE TABLE consequences (
			csq_id			INTEGER PRIMARY KEY,
			var_id			INTEGER,
			transcript_id	TEXT,
			gene_id			TEXT,
			ext_id			TEXT,
			uniprot_id		TEXT,
			protein_id		TEXT,
			protein_pos		TEXT,
			aa_change		TEXT,
			sift_score		REAL,
			sift_tfic		REAL,
			sift_tfic_class	INTEGER,
			pph2_score		REAL,
			pph2_tfic		REAL,
			pph2_tfic_class	INTEGER,
			ma_score		REAL,
			ma_tfic			REAL,
			ma_tfic_class	INTEGER,
			impact			INTEGER,

			UNIQUE (var_id, transcript_id)
		);

		CREATE INDEX csq_gene_id ON consequences (gene_id);

		CREATE TABLE consequence_types (
			csq_id          INTEGER,
			so_id           INTEGER,

			UNIQUE (csq_id, so_id)
		);

		CREATE TABLE so_terms (
			so_id           INTEGER PRIMARY KEY,
			name            TEXT UNIQUE
		);
		""")

		# features

		c.executescript("""
		CREATE TABLE genes (
			gene_id			TEXT PRIMARY KEY,
			symbol			TEXT,
			fm_pvalue		REAL,
			fm_qvalue		REAL,
			fm_exc_cause	TEXT DEFAULT '{0}',
			clust_coords	TEXT,
			clust_zscore	REAL,
			clust_pvalue	REAL,
			clust_qvalue	REAL,
			clust_exc_cause	TEXT DEFAULT '{0}'
		);""".format(self.GENE_EXC_NOT_FOUND))

		c.executescript("""
		CREATE TABLE gene_xrefs (
			gene_id			TEXT,
			xref			TEXT,
			
			PRIMARY KEY (gene_id, xref)
		);

		CREATE TABLE pathways (
			pathway_id		TEXT PRIMARY KEY,
			desc			TEXT,
			gene_count		INTEGER,
			fm_zscore		REAL,
			fm_pvalue		REAL,
			fm_qvalue		REAL
		);

		CREATE TABLE pathway_genes (
			pathway_id		INTEGER,
			gene_id			INTEGER,

			PRIMARY KEY (pathway_id, gene_id)
		);

		CREATE INDEX pathway_genes_by_gene ON pathway_genes (gene_id);
		""")

		# oncodrivefm

		c.executescript("""
		CREATE TABLE sample_gene_fimpact (
			sample_id		INTEGER,
			gene_id			TEXT,
			sift_score		REAL,
			sift_tfic		REAL,
			sift_tfic_class	INTEGER,
			pph2_score		REAL,
			pph2_tfic		REAL,
			pph2_tfic_class INTEGER,
			ma_score		REAL,
			ma_tfic			REAL,
			ma_tfic_class	INTEGER,

			PRIMARY KEY (sample_id, gene_id)
		);
		""")

		# recurrences

		c.executescript("""
		CREATE TABLE aff_gene_rec (
			afg_id			INTEGER PRIMARY KEY,
			sample_freq		INTEGER,
			sample_prop		REAL
		);

		CREATE TABLE gene_rec (
			gene_id			TEXT PRIMARY KEY,
			sample_freq		INTEGER,
			sample_prop		REAL
		);

		CREATE TABLE pathway_rec (
			pathway_id		TEXT PRIMARY KEY,
			sample_freq		INTEGER,
			sample_prop		REAL,
			gene_freq		INTEGER,
			gene_prop		REAL
		);
		""")

		c.close()

		self.__conn.commit()

	def close(self):
		self.__conn.close()

	def commit(self):
		self.__conn.commit()

	def rollback(self):
		self.__conn.rollback()

	# sources ------------------------------------

	def add_source(self, source):
		c = self.__conn.cursor()
		c.execute("INSERT INTO sources (name) VALUES (?)", (source,))
		source_id = c.lastrowid
		c.close()
		return source_id

	def add_source_line(self, source_id, line_num, text):
		c = self.__conn.cursor()
		try:
			c.execute("INSERT INTO source_lines (source_id, line_num, text) VALUES (?,?,?)", (source_id, line_num, text))
		except sqlite3.IntegrityError as ex:
			pass
		c.close()

	# samples ------------------------------------

	def add_sample(self, sample):
		c = self.__conn.cursor()
		try:
			c.execute("INSERT INTO samples (name) VALUES (?)", (sample.name,))
			sample_id = c.lastrowid
		except sqlite3.IntegrityError as ex:
			c.execute("SELECT sample_id FROM samples WHERE name=?", (sample.name,))
			res = c.fetchone()
			if res is None:
				raise ex
			sample_id = res[0]
		c.close()
		return sample_id

	def get_sample_by_name(self, name):
		c = self.__conn.cursor()
		c.execute("SELECT sample_id FROM samples WHERE name=?", (name,))
		res = c.fetchone()
		c.close()
		if res is None:
			return None
		sample = Sample(id=res[0], name=name)
		return sample

	# variants -----------------------------------

	def add_variant(self, var, source_id=None, line_num=None):
		c = self.__conn.cursor()
		try:
			c.execute("INSERT INTO variants (type, chr, strand, start, ref, alt) VALUES (?, ?, ?, ?, ?, ?)",
				(var.type, var.chr, var.strand, var.start, var.ref, var.alt))
			var_id = c.lastrowid
		except sqlite3.IntegrityError as ex:
			c.execute("SELECT var_id FROM variants WHERE type=? AND chr=? AND strand=? AND start=? AND ref=? AND alt=?",
				(var.type, var.chr, var.strand, var.start, var.ref, var.alt))
			res = c.fetchone()
			if res is None:
				raise ex
			var_id = res[0]

		if var.samples is not None:
			for sample in var.samples:
				sample_id = self.add_sample(sample)
				try:
					c.execute("INSERT INTO variant_samples (var_id, sample_id) VALUES (?, ?)", (var_id, sample_id))
				except sqlite3.IntegrityError as ex:
					pass

		if source_id is not None and line_num is not None:
			c.execute("INSERT INTO source_variants (source_id, line_num, var_id) VALUES (?,?,?)", (source_id, line_num, var_id))

		c.close()
		return var_id

	def update_variant(self, var):
		vars = [("chr", var.chr), ("strand", var.strand), ("start", var.start), ("ref", var.ref), ("alt", var.alt)]

		update = []
		params = []

		for t in vars:
			if t[1] is not None:
				update += ["{0}=?".format(t[0])]
				params += [t[1]]

		c = self.__conn.cursor()

		if len(update) > 0:
			if var.id is not None:
				where = "var_id=?"
				params += [var.id]
			elif var.chr is not None and var.start is not None and var.ref is not None and var.alt is not None:
				where = "chr=? AND start=? AND ref=? AND alt=?"
				params += [var.id, var.start, var.ref, var.alt]
			else:
				raise Exception("Either id or chr, start, ref and alt should be defined")

			sql = " ".join(["UPDATE variants SET",  ", ".join(update), "WHERE", where])
			c.execute(sql, tuple(params))

		if var.samples is not None or var.xrefs is not None:
			if var.id is not None:
				var_id = var.id
			elif var.chr is not None and var.start is not None and var.ref is not None and var.alt is not None:
				c.execute("SELECT var_id FROM variants WHERE chr=? AND start=? AND ref=? AND alt=?", (var.id, var.start, var.ref, var.alt))
				var_id, = c.fetchone()

		if var.samples is not None:
			c.execute("DELETE FROM variant_samples WHERE var_id=?", (var.id, ))
			for sample in var.samples:
				if sample.id is None:
					self.get_sample_by_name(sample.name)
				c.execute("INSERT INTO variant_samples (var_id, sample_id) VALUES (?,?)", (var_id, sample.id))
				
		if var.xrefs is not None:
			c.execute("DELETE FROM variant_xrefs WHERE var_id=?", (var.id, ))
			for xref in set(var.xrefs):
				c.execute("INSERT INTO variant_xrefs (var_id, xref) VALUES (?,?)", (var_id, xref))

		c.close()

	def update_variant_start(self, var_id, start):
		c = self.__conn.cursor()
		c.execute("UPDATE variants SET prev_start=start, start=? WHERE var_id=?", (start, var_id))
		c.close()

	#DEPRECATED Variants need to be tracked
	def remove_variant(self, var_id):
		c = self.__conn.cursor()
		c.execute("DELETE FROM variant_samples WHERE var_id=?", (var_id,))
		c.execute("DELETE FROM samples WHERE sample_id NOT IN (SELECT DISTINCT sample_id FROM variant_samples)")
		#c.execute("DELETE FROM variant_sources WHERE var_id=?", (var_id,))
		#c.execute("DELETE FROM sources WHERE source_id NOT IN (SELECT DISTINCT source_id FROM variant_sources)")
		c.execute("DELETE FROM variants WHERE var_id=?", (var_id,))
		c.close()

	def get_variant_samples(self, var_id):
		c = self.__conn.cursor()
		c.execute("SELECT sample_id, name FROM variant_samples LEFT JOIN samples USING (sample_id) WHERE var_id=? ORDER BY name", (var_id,))
		samples = [Sample(id=r[0], name=r[1]) for r in c]
		c.close()
		return samples

	def get_variant_xrefs(self, var_id):
		c = self.__conn.cursor()
		c.execute("SELECT xref FROM variant_xrefs WHERE var_id=? ORDER BY xref", (var_id,))
		xrefs = [r[0] for r in c]
		c.close()
		return xrefs

	def get_variant_genes(self, var_id):
		c = self.__conn.cursor()
		c.execute("SELECT DISTINCT gene_id FROM affected_genes WHERE var_id=?", (var_id,))
		genes = [r[0] for r in c]
		c.close()
		return genes

	def get_variant(self, var_id, join_samples=False, join_xrefs=False):
		c = self.__conn.cursor()
		c.execute("SELECT * FROM variants WHERE var_id=?", (var_id,))
		r = c.fetchone()
		c.close()
		if r is None:
			return None

		samples = xrefs = None

		if join_samples:
			samples = self.get_variant_samples(r["var_id"])

		if join_xrefs:
			xrefs = self.get_variant_xrefs(r["var_id"])

		return Variant(id=r["var_id"], type=r["type"], chr=r["chr"], strand=r["strand"],
			start=r["start"], ref=r["ref"], alt=r["alt"], samples=samples, xrefs=xrefs)

	def variants(self, join_samples=False, join_xrefs=False, order_by="position"):
		c = self.__conn.cursor()
			
		sql = ["SELECT * FROM variants v"]
		sql += ["WHERE v.start IS NOT NULL"]

		if order_by == "position":
			sql += ["ORDER BY chr*1, chr, strand, start"]
		elif order_by == "id":
			sql += ["ORDER BY var_id"]

		samples = xrefs = None

		for r in c.execute(" ".join(sql)):
			if join_samples:
				samples = self.get_variant_samples(r["var_id"])

			if join_xrefs:
				xrefs = self.get_variant_xrefs(r["var_id"])

			yield Variant(
					id=r["var_id"],
					type=r["type"], chr=r["chr"], start=r["start"],
					ref=r["ref"], alt=r["alt"], strand=r["strand"],
					samples=samples, xrefs=xrefs)

		c.close()

	def count_variants(self):
		c = self.__conn.cursor()

		# source counts

		c.execute("SELECT name, count(line_num) FROM source_lines JOIN sources USING (source_id) GROUP BY name")
		source_counts = {}
		for source, count in c:
			source_counts[source] = count

		# parser counts

		c.execute("SELECT name, count(var_id) FROM source_variants JOIN sources USING (source_id) GROUP BY name")
		parser_counts = {}
		for source, count in c:
			parser_counts[source] = count

		# liftover counts

		c.execute("""
			SELECT name, count(var_id) FROM variants
			JOIN source_variants USING (var_id)
			JOIN sources USING (source_id)
			WHERE start is not NULL
			GROUP BY name""")
		liftover_counts = {}
		for source, count in c:
			liftover_counts[source] = count

		# vep counts

		c.execute("""
			SELECT name, count(distinct var_id) FROM consequences
			JOIN source_variants USING (var_id)
			JOIN sources USING (source_id)
			GROUP BY name""")
		vep_counts = {}
		for source, count in c:
			vep_counts[source] = count

		c.execute("SELECT name FROM sources")
		sources = [s for s, in c.fetchall()]
		c.close()

		counts = []
		for source in sources:
			source_passed = source_counts.get(source, 0)
			parser_passed = parser_counts.get(source, 0)
			parser_lost = max(0, source_passed - parser_passed)
			liftover_passed = liftover_counts.get(source, 0)
			liftover_lost = max(0, parser_passed - liftover_passed)
			vep_passed = vep_counts.get(source, 0)
			vep_lost = max(0, liftover_passed - vep_passed)
			total_lost = parser_lost + liftover_lost + vep_lost
			counts += [(source, dict(
				source_passed=source_passed,
				parser_passed=parser_passed,
				parser_lost=parser_lost,
				#parser_lost_pct=float(parser_lost) * 100.0 / total_lost,
				liftover_passed=liftover_passed,
				liftover_lost=liftover_lost,
				#liftover_lost_pct=float(liftover_lost) * 100.0 / total_lost,
				vep_passed=vep_passed,
				vep_lost=vep_lost,
				#vep_lost_pct=float(vep_lost) * 100.0 / total_lost,
				total_lost=total_lost))]
				#total_lost_pct=float(total_lost) * 100.0 / source_passed))]
		return counts

	def add_affected_gene(self, afg):
		c = self.__conn.cursor()
		try:
			c.execute(
				"INSERT INTO affected_genes (var_id, gene_id, impact, coding_region, prot_changes) VALUES (?,?,?,?,?)",
				(afg.var.id, afg.gene_id, afg.impact, afg.coding_region, afg.prot_changes))
			afg_id = c.lastrowid
		except sqlite3.IntegrityError as ex:
			c.execute("SELECT afg_id FROM affected_genes WHERE var_id=? AND gene_id=?", (afg.var.id, afg.gene_id))
			res = c.fetchone()
			if res is None:
				raise ex
			afg_id = res[0]
			c.execute("UPDATE affected_genes SET impact=?, coding_region=?, prot_changes=? WHERE afg_id=?",
					  (afg.impact, afg.coding_region, afg.prot_changes, afg_id))

		c.close()
		return afg_id

	def affected_genes(self, join_variant=False, join_samples=False, join_xrefs=False, join_rec=False):
		c = self.__conn.cursor()

		sql = ["SELECT * FROM affected_genes ag"]
		if join_rec:
			sql += ["LEFT JOIN aff_gene_rec USING (afg_id)"]

		sql += ["ORDER BY var_id, gene_id"]

		var = rec = None

		for r in c.execute(" ".join(sql)):
			if join_variant:
				var = self.get_variant(r["var_id"], join_samples=join_samples, join_xrefs=join_xrefs)
			else:
				var = Variant(id=r["var_id"])

			if join_rec:
				rec = AffectedGeneRec(afg_id=r["afg_id"], sample_freq=r["sample_freq"], sample_prop=r["sample_prop"])

			yield AffectedGene(
				id=r["afg_id"], var=var, gene_id=r["gene_id"], impact=r["impact"],
				coding_region=r["coding_region"], prot_changes=r["prot_changes"], rec=rec)

		c.close()

	def get_total_affected_samples(self):
		c = self.__conn.cursor()
		c.execute("""SELECT COUNT(DISTINCT vs.sample_id) FROM affected_genes
		 				LEFT JOIN variant_samples vs USING (var_id)""")
		res = c.fetchone()
		c.close()
		if res is None:
			return 0
		return res[0]

	# consequences ------------------------------------------------------

	def add_consequence(self, csq):
		c = self.__conn.cursor()

		try:
			c.execute("""
				INSERT INTO consequences (var_id, transcript_id, gene_id,
											ext_id, uniprot_id, protein_id, protein_pos, aa_change,
											sift_score, sift_tfic, sift_tfic_class,
											pph2_score, pph2_tfic, pph2_tfic_class,
											ma_score, ma_tfic, ma_tfic_class, impact)
					VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
				csq.var.id, csq.transcript, csq.gene,
				csq.extid, csq.uniprot, csq.protein, csq.protein_pos, csq.aa_change,
				csq.sift_score, csq.sift_tfic, csq.sift_tfic_class,
				csq.pph2_score, csq.pph2_tfic, csq.pph2_tfic_class,
				csq.ma_score, csq.ma_tfic, csq.ma_tfic_class, csq.impact))
			csq_id = c.lastrowid
		except sqlite3.IntegrityError as ex:
			"""
			raise Exception("Consequence already exists: ({0})".format(repr(csq))
			"""
			return

		so_ids = []
		for ctype in csq.ctypes:
			c.execute("SELECT so_id FROM so_terms WHERE name=?", (ctype,))
			res = c.fetchone()
			if res is None:
				c.execute("INSERT INTO so_terms (name) VALUES (?)", (ctype,))
				so_id = c.lastrowid
			else:
				so_id = res[0]
			so_ids += [so_id]

		for so_id in so_ids:
			c.execute("INSERT INTO consequence_types (csq_id, so_id) VALUES (?, ?)", (csq_id, so_id))

		c.close()

	def update_consequence(self, csq):

		vars = [("gene_id", csq.gene), ("ext_id", csq.extid), ("uniprot_id", csq.uniprot),
				("protein_id", csq.protein), ("protein_pos", csq.protein_pos), ("aa_change", csq.aa_change),
				("sift_score", csq.sift_score), ("sift_tfic", csq.sift_tfic), ("sift_tfic_class", csq.sift_tfic_class),
				("pph2_score", csq.pph2_score), ("pph2_tfic", csq.pph2_tfic), ("pph2_tfic_class", csq.pph2_tfic_class),
				("ma_score", csq.ma_score), ("ma_tfic", csq.ma_tfic), ("ma_tfic_class", csq.ma_tfic_class),
				("impact", csq.impact)]

		update = []
		params = []

		for t in vars:
			if t[1] is not None:
				update += ["{0}=?".format(t[0])]
				params += [t[1]]

		if len(update) == 0:
			return

		if csq.id is not None:
			where = "csq_id=?"
			params += [csq.id]
		elif csq.var.id is not None and csq.transcript is not None:
			where = "var_id=? AND transcript_id=?"
			params += [csq.var.id, csq.transcript]
		elif csq.var.id is not None:
			where = "var_id=?"
			params += [csq.var.id]
		else:
			raise Exception("Either csq_id, only var_id or var_id and transcript should be defined")

		sql = " ".join(["UPDATE consequences SET",  ", ".join(update), "WHERE", where])
		c = self.__conn.cursor()
		c.execute(sql, tuple(params))
		c.close()

	def consequences(self, join_variant=False, join_samples=False, join_xrefs=False, join_ctypes=True, filters=None):
		"""
		Retrieve consequences.
		:param join_variant: Whether to join the variants table. Attribute var will have all the variant information.
		:param join_samples: Whether to join the variant samples table. Attribute samples will contain a list of Samples.
		:param join_xrefs: Whether to join external references for variants table. Attribute xrefs will contain a list of xrefs.
		:param join_ctypes: Whether to join consequence types table. Attribute ctypes will contain a list of SO terms.
		:param filters: Dictionary with filters:
		                **ctypes**: Value = list of SO terms. Retrieve only those consequences with any of its ctypes in values.
		                **var**: Value = var_id. Retrieve the consequences for the given variant.
		                **gene**: Value = gene_id. Retrieve the consequences for the given gene.
		:return: Iterator for the retrieved consequences.
		"""

		sql = ["SELECT * FROM consequences c"]
		if join_variant:
			sql += ["LEFT JOIN variants v USING (var_id)"]

		params = []
		if filters is not None:
			filters_sql = []
			for filter, value in filters.items():
				if filter == self.CSQ_CTYPES:
					filters_sql += ["""
					csq_id IN (
						SELECT DISTINCT csq_id FROM consequence_types ct
						JOIN so_terms so USING (so_id)
						WHERE so.name IN ({0}))""".format(",".join(["?"] * len(value)))]
					params += list(value)
					# WHERE so.name IN ({0}))""".format(",".join(["'{0}'".format(term) for term in value]))]
				elif filter == self.CSQ_VAR:
					filters_sql += ["var_id = ?"]
					params += [value]
				elif filter == self.CSQ_GENE:
					filters_sql += ["gene_id = ?"]
					params += [value]

			if len(filters_sql) > 0:
				sql += ["WHERE ", " AND ".join(filters_sql)]

		for r in self.__conn.execute(" ".join(sql), tuple(params)):
			so_terms = None

			if join_ctypes:
				sql = "SELECT so.name FROM consequence_types t LEFT JOIN so_terms so USING (so_id) WHERE t.csq_id=?"
				so_terms = [n[0] for n in self.__conn.execute(sql, (r["csq_id"],))]

			samples = xrefs = None

			if join_samples:
				samples = self.get_variant_samples(r["var_id"])

			if join_xrefs:
				xrefs = self.get_variant_xrefs(r["var_id"])

			if join_variant:
				var = Variant(
					id=r["var_id"],
					type=r["type"], chr=r["chr"], start=r["start"],
					ref=r["ref"], alt=r["alt"], strand=r["strand"],
					samples=samples, xrefs=xrefs)
			else:
				var = Variant(id=r["var_id"], samples=samples, xrefs=xrefs)

			csq = Consequence(
				id=r["csq_id"], var=var, transcript=r["transcript_id"], gene=r["gene_id"], ctypes=so_terms,
				extid=r["ext_id"], protein=r["protein_id"], uniprot=r["uniprot_id"],
				protein_pos=r["protein_pos"], aa_change=r["aa_change"],
				sift_score=r["sift_score"], sift_tfic=r["sift_tfic"], sift_tfic_class=r["sift_tfic_class"],
				pph2_score=r["pph2_score"], pph2_tfic=r["pph2_tfic"], pph2_tfic_class=r["pph2_tfic_class"],
				ma_score=r["ma_score"], ma_tfic=r["ma_tfic"], ma_tfic_class=r["ma_tfic_class"], impact=r["impact"])

			yield csq

	# genes --------------------------------------------------------------

	def load_genes(self, path):
		c = self.__conn.cursor()
		with open(path, "r") as f:
			hdr = {}
			for i, name in enumerate(f.readline().rstrip("\n").split("\t")):
				hdr[name] = i

			ensg_index = hdr["Ensembl Gene ID"]
			symbol_index = hdr["HGNC symbol"]
			for line in f:
				fields = line.rstrip("\n").split("\t")
				gene_id = fields[ensg_index]
				symbol = fields[symbol_index]

				try:
					c.execute("INSERT INTO genes (gene_id, symbol) VALUES (?,?)", (gene_id, symbol))
				except sqlite3.IntegrityError as ex:
					pass

		c.close()

	def get_gene_symbols(self):
		c = self.__conn.cursor()
		gene_symbol = {}
		for gene_id, symbol in c.execute("SELECT gene_id, symbol FROM genes"):
			gene_symbol[gene_id] = symbol
		c.close()
		return gene_symbol

	def update_gene(self, gene):
		vars = [("symbol", gene.symbol), ("fm_pvalue", gene.fm_pvalue), ("fm_qvalue", gene.fm_qvalue),
				("fm_exc_cause", gene.fm_exc_cause), ("clust_coords", gene.clust_coords),
				("clust_zscore", gene.clust_zscore), ("clust_pvalue", gene.clust_pvalue),
				("clust_qvalue", gene.clust_qvalue), ("clust_exc_cause", gene.clust_exc_cause)]

		update = []
		params = []

		for t in vars:
			if t[1] is not None:
				update += ["{0}=?".format(t[0])]
				params += [t[1]]

		if gene.id is not None:
			where = "gene_id=?"
			params += [gene.id]
		else:
			raise Exception("id should be defined")

		c = self.__conn.cursor()

		if len(update) > 0:
			sql = " ".join(["UPDATE genes SET",  ", ".join(update), "WHERE", where])
			print sql
			c.execute(sql, tuple(params))

		if gene.xrefs is not None:
			c.execute("DELETE FROM gene_xrefs WHERE gene_id=?", (gene.id, ))
			for xref in gene.xrefs:
				c.execute("INSERT INTO gene_xrefs (gene_id, xref) VALUES (?,?)", (gene.id, xref))

		c.close()

	def get_gene_xrefs(self, gene_id):
		c = self.__conn.cursor()
		c.execute("SELECT xref FROM gene_xrefs WHERE gene_id=? ORDER BY xref", (gene_id,))
		xrefs = [r[0] for r in c]
		c.close()
		return xrefs

	def genes(self, join_xrefs=False, join_rec=False):
		c = self.__conn.cursor()

		xrefs = rec = None

		sql = ["SELECT * FROM genes"]
		if join_rec:
			sql += ["LEFT JOIN gene_rec USING (gene_id)"]

		for r in c.execute(" ".join(sql)):
			if join_xrefs:
				xrefs = self.get_gene_xrefs(r["gene_id"])

			if join_rec:
				rec = GeneRec(gene_id=r["gene_id"], sample_freq=r["sample_freq"], sample_prop=r["sample_prop"])

			yield Gene(
					id=r["gene_id"], symbol=r["symbol"], xrefs=xrefs,
					fm_pvalue=r["fm_pvalue"], fm_qvalue=r["fm_qvalue"], fm_exc_cause=r["fm_exc_cause"],
					clust_coords=r["clust_coords"], clust_zscore=r["clust_zscore"],
					clust_pvalue=r["clust_pvalue"], clust_qvalue=r["clust_qvalue"],
					clust_exc_cause=r["clust_exc_cause"], rec=rec)
		c.close()

	def delete_sample_gene_fimpact(self):
		c = self.__conn.cursor()
		c.execute("DELETE FROM sample_gene_fimpact")

	def add_sample_gene_fimpact(self, sample_id, gene_id, sift_score, sift_tfic, sift_tfic_class,
								pph2_score, pph2_tfic, pph2_tfic_class, ma_score, ma_tfic, ma_tfic_class):
		c = self.__conn.cursor()
		c.execute("""
		INSERT INTO sample_gene_fimpact (sample_id, gene_id, sift_score, sift_tfic, sift_tfic_class,
			 							pph2_score, pph2_tfic, pph2_tfic_class, ma_score, ma_tfic, ma_tfic_class)
		VALUES (?,?,?,?,?,?,?,?,?,?,?)""", (sample_id, gene_id, sift_score, sift_tfic, sift_tfic_class,
											pph2_score, pph2_tfic, pph2_tfic_class, ma_score, ma_tfic, ma_tfic_class))
		c.close()

	def sample_gene_fimpacts(self):
		c = self.__conn.cursor()
		c.execute("""
			SELECT gene_id, s.name AS sample_name,
					sift_score, sift_tfic, sift_tfic_class,
					pph2_score, pph2_tfic, pph2_tfic_class,
					ma_score, ma_tfic, ma_tfic_class
			FROM sample_gene_fimpact
			LEFT JOIN samples s USING (sample_id)
			ORDER BY sample_name""")

		for r in c:
			yield (r["gene_id"], r["sample_name"],
					r["sift_score"], r["sift_tfic"], r["sift_tfic_class"],
					r["pph2_score"], r["pph2_tfic"], r["pph2_tfic_class"],
					r["ma_score"], r["ma_tfic"], r["ma_tfic_class"])

		c.close()

	# pathways -----------------------------------------------------------

	def load_pathways(self, desc_path, genes_map_path):
		c = self.__conn.cursor()

		pathway_genes = {}
		with open(genes_map_path, "r") as f:
			for line in f:
				gene_id, pathway_id = line.rstrip("\n").split("\t")

				if pathway_id not in pathway_genes:
					pathway_genes[pathway_id] = 1
				else:
					pathway_genes[pathway_id] += 1

				c.execute("INSERT INTO pathway_genes (pathway_id, gene_id) VALUES (?,?)", (pathway_id, gene_id))

		with open(desc_path, "r") as f:
			for line in f:
				pathway_id, desc = line.rstrip("\n").split("\t")

				genes_count = pathway_genes.get(pathway_id, 0)

				c.execute("INSERT INTO pathways (pathway_id, desc, gene_count) VALUES (?,?,?)",
												(pathway_id, desc, genes_count))

		c.close()

	def update_pathway(self, pathway):
		vars = [("gene_count", pathway.gene_count), ("fm_zscore", pathway.fm_zscore), ("fm_pvalue", pathway.fm_pvalue),
				("fm_qvalue", pathway.fm_qvalue)]

		update = []
		params = []

		for t in vars:
			if t[1] is not None:
				update += ["{0}=?".format(t[0])]
				params += [t[1]]

		if pathway.id is not None:
			where = "pathway_id=?"
			params += [pathway.id]
		else:
			raise Exception("id should be defined")

		sql = " ".join(["UPDATE pathways SET",  ", ".join(update), "WHERE", where])
		c = self.__conn.cursor()
		c.execute(sql, tuple(params))
		c.close()

	def pathways(self, join_rec=False):
		c = self.__conn.cursor()
		sql = ["SELECT * FROM pathways"]
		if join_rec:
			sql += ["LEFT JOIN pathway_rec USING (pathway_id)"]

		rec = None

		for r in c.execute(" ".join(sql)):
			if join_rec:
				rec=PathwayRec(
						pathway_id=r["pathway_id"], sample_freq=r["sample_freq"], sample_prop=r["sample_prop"],
						gene_freq=r["gene_freq"], gene_prop=r["gene_prop"])

			yield Pathway(
				id=r["pathway_id"], gene_count=r["gene_count"],
				fm_zscore=r["fm_zscore"], fm_pvalue=r["fm_pvalue"], fm_qvalue=r["fm_qvalue"],
				rec=rec)

		c.close()

	# recurrences --------------------------------------------------------

	def compute_affected_genes_recurrences(self, total_samples):
		c = self.__conn.cursor()
		c.execute("DELETE FROM aff_gene_rec")
		c.execute("""
		INSERT INTO aff_gene_rec (afg_id, sample_freq, sample_prop)
			SELECT ag.afg_id,
				COUNT(DISTINCT vs.sample_id),
				COUNT(DISTINCT vs.sample_id) / {0}.0
			FROM affected_genes ag
			JOIN variant_samples vs USING (var_id)
			WHERE coding_region = 1
			GROUP BY ag.var_id, ag.gene_id
		""".format(total_samples))
		rowcount = c.rowcount
		c.close()
		return rowcount

	def compute_gene_recurrences(self, total_samples):
		c = self.__conn.cursor()
		c.execute("DELETE FROM gene_rec")
		c.execute("""
		INSERT INTO gene_rec (gene_id, sample_freq, sample_prop)
			SELECT ag.gene_id,
				COUNT(DISTINCT vs.sample_id),
				COUNT(DISTINCT vs.sample_id) / {0}.0
			FROM affected_genes ag
			JOIN variant_samples vs USING (var_id)
			GROUP BY ag.gene_id
		""".format(total_samples))
		rowcount = c.rowcount
		c.close()
		return rowcount

	def compute_pathway_recurrences(self, total_samples):
		c = self.__conn.cursor()

		c.executescript("""
		DROP TABLE IF EXISTS gene_samples;
		CREATE TABLE gene_samples (gene_id TEXT, sample_id INTEGER, PRIMARY KEY (gene_id, sample_id));
		INSERT INTO gene_samples
			SELECT DISTINCT gene_id, sample_id
			FROM affected_genes
			JOIN variant_samples USING (var_id);

		DROP TABLE IF EXISTS pathway_samples;
		CREATE TABLE pathway_samples (pathway_id TEXT PRIMARY KEY, sample_freq INTEGER, genes_sample_freq INTEGER);
		INSERT INTO pathway_samples
			SELECT p.pathway_id, COUNT(DISTINCT sample_id), COUNT(sample_id)
			FROM pathways p
			JOIN pathway_genes pg ON (p.pathway_id = pg.pathway_id)
			JOIN gene_samples gs ON (pg.gene_id = gs.gene_id)
			GROUP BY p.pathway_id;
		""")

		c.execute("DELETE FROM pathway_rec")
		c.execute("""
		INSERT INTO pathway_rec (pathway_id, sample_freq, sample_prop, gene_freq, gene_prop)
			SELECT
				p.pathway_id,
				ps.sample_freq,
				ps.sample_freq / {0}.0,
				ps.genes_sample_freq,
				(ps.sample_freq / {0}.0) / p.gene_count
			FROM pathways p
			JOIN pathway_samples ps USING (pathway_id);
		""".format(total_samples))
		rowcount = c.rowcount

		c.execute("UPDATE pathway_rec SET gene_freq=0.0 WHERE gene_freq IS NULL")

		c.close()
		return rowcount
