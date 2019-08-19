import os
import sqlite3

from wok.task import task

from intogensm.woktmp import make_temp_file, remove_temp
from bgcore import tsv

from intogensm.projdb import ProjectDb
from intogensm.utils import normalize_id
from intogensm.paths import get_combination_path
from intogensm.variants.utils import var_to_tab, complementary_sequence
from intogensm.transfic import TransFIC

def create_db(conn):
	conn.executescript("""
	CREATE TABLE variants (
		var_id			INTEGER PRIMARY KEY,
		chr				TEXT,
		strand			TEXT,
		start			INTEGER,
		ref				TEXT,
		alt				TEXT,

		xrefs			TEXT,

		UNIQUE (chr, strand, start, ref, alt)
	);

	CREATE TABLE variant_genes (
		var_id			INTEGER,
		gene_id			TEXT,

		impact			INTEGER,
		coding_region	INTEGER,
		prot_changes	TEXT,

		sample_freq		INTEGER,
		sample_prop		REAL,

		PRIMARY KEY (var_id, gene_id)
	);

	CREATE TABLE genes (
		gene_id			TEXT PRIMARY KEY,

		sample_freq		INTEGER,
		sample_prop		REAL
	);

	CREATE TABLE pathways (
		pathway_id		TEXT PRIMARY KEY,
		
		sample_freq		INTEGER,
		sample_prop		REAL
	);
	""")

	conn.commit()

@task.foreach()
def combination_recurrences(projects_set):
	log = task.logger
	conf = task.conf

	classifier, projects = projects_set

	classifier_id = classifier["id"]

	group_values = classifier["group_values"]
	short_values = classifier["group_short_values"]
	long_values = classifier["group_long_values"]

	group_name = classifier["group_name"]
	group_short_name = classifier["group_short_name"]
	group_long_name = classifier["group_long_name"]

	if len(group_values) == 0:
		group_file_prefix = classifier_id
	else:
		group_file_prefix = "{0}-{1}".format(classifier_id, group_short_name)

	group_file_prefix = normalize_id(group_file_prefix)

	log.info("--- [{0} ({1}) ({2}) ({3})] {4}".format(classifier["name"], group_long_name, group_short_name, group_name, "-" * 30))

	log.info("Creating database ...")

	db_path = make_temp_file(task, suffix="-{0}.db".format(group_file_prefix))
	log.debug("  > {0}".format(db_path))

	conn = sqlite3.connect(db_path)
	conn.row_factory = sqlite3.Row

	create_db(conn)

	log.info("Combining recurrences ...")

	c = conn.cursor()

	sample_total = 0

	project_ids = []
	for project in projects:
		project_ids += [project["id"]]

		log.info("  Project {0}:".format(project["id"]))

		projdb = ProjectDb(project["db"])

		project_sample_total = projdb.get_total_affected_samples()

		sample_total += project_sample_total

		log.info("    Total samples = {0}".format(project_sample_total))

		log.info("    Variant genes ...")

		count = 0
		for afg in projdb.affected_genes(join_variant=True, join_xrefs=True, join_rec=True):
			var = afg.var
			rec = afg.rec

			if rec.sample_freq is None:
				log.warn("Discarding variant gene without sample frequency: {0}".format(repr(afg)))
				continue

			start, end, ref, alt = var_to_tab(var)

			try:
				c.execute("INSERT INTO variants (chr, strand, start, ref, alt, xrefs) VALUES (?,?,?,?,?,?)",
						  (var.chr, var.strand, start, ref, alt, ",".join(var.xrefs)))
				var_id = c.lastrowid
			except sqlite3.IntegrityError:
				c.execute("SELECT var_id FROM variants WHERE chr=? AND strand=? AND start=? AND ref=? AND alt=?",
						  (var.chr, var.strand, start, ref, alt))
				r = c.fetchone()
				var_id = r[0]

			try:
				c.execute("INSERT INTO variant_genes (var_id, gene_id, impact, coding_region, prot_changes, sample_freq) VALUES (?,?,?,?,?,?)",
						  (var_id, afg.gene_id, afg.impact, afg.coding_region, afg.prot_changes, rec.sample_freq))
			except sqlite3.IntegrityError:
				c.execute("""
					UPDATE variant_genes
					SET sample_freq=sample_freq + ?
					WHERE var_id=? AND gene_id=?""",
						(rec.sample_freq, var_id, afg.gene_id))

			count += 1

		log.info("      {0} variant genes".format(count))

		log.info("    Genes ...")

		count = 0
		for gene in projdb.genes(join_xrefs=True, join_rec=True):
			rec = gene.rec

			if rec.sample_freq is None:
				continue

			c.execute("SELECT COUNT(*) FROM genes WHERE gene_id=?", (gene.id,))
			r = c.fetchone()
			if r[0] == 0:
				c.execute("INSERT INTO genes (gene_id, sample_freq) VALUES (?,?)",
					  (gene.id, rec.sample_freq))
			else:
				c.execute("UPDATE genes SET sample_freq=sample_freq + ? WHERE gene_id=?",
						  (rec.sample_freq, gene.id))
			count += 1

		log.info("      {0} genes".format(count))

		log.info("    Pathways ...")

		count = 0
		for pathway in projdb.pathways(join_rec=True):
			rec = pathway.rec

			if rec.sample_freq is None:
				continue

			c.execute("SELECT COUNT(*) FROM pathways WHERE pathway_id=?", (pathway.id,))
			r = c.fetchone()
			if r[0] == 0:
				c.execute("INSERT INTO pathways (pathway_id, sample_freq) VALUES (?,?)",
						  (pathway.id, rec.sample_freq))
			else:
				c.execute("UPDATE pathways SET sample_freq=sample_freq + ? WHERE pathway_id=?",
						  (rec.sample_freq, pathway.id))
			count += 1

		log.info("      {0} pathways".format(count))

		projdb.close()

	log.info("Calculating proportions with {0} samples in total among projects ...".format(sample_total))

	if sample_total > 0:
		c.execute("UPDATE variant_genes SET sample_prop=CAST(sample_freq AS REAL)/{0}.0".format(sample_total))
		c.execute("UPDATE genes SET sample_prop=CAST(sample_freq AS REAL)/{0}.0".format(sample_total))
		c.execute("UPDATE pathways SET sample_prop=CAST(sample_freq AS REAL)/{0}.0".format(sample_total))

	c.close()
	conn.commit()
	
	log.info("Saving results ...")
	
	c = conn.cursor()

	base_path = get_combination_path(conf, "recurrences")

	log.info("  Variant genes ...")

	with tsv.open(os.path.join(base_path, "variant_gene-{0}.tsv.gz".format(group_file_prefix)), "w") as f:
		tsv.write_param(f, "classifier", classifier["id"])
		tsv.write_param(f, "group_id", group_name)
		tsv.write_param(f, "group_short_name", group_short_name)
		tsv.write_param(f, "group_long_name", group_long_name)
		tsv.write_param(f, "projects", ",".join(project_ids))
		tsv.write_param(f, "SAMPLE_TOTAL", sample_total)
		tsv.write_line(f, "CHR", "STRAND", "START", "ALLELE", "GENE_ID", "IMPACT", "IMPACT_CLASS", "SAMPLE_FREQ", "SAMPLE_PROP", "PROT_CHANGES", "XREFS")
		for r in c.execute("SELECT * FROM variant_genes JOIN variants USING (var_id) ORDER BY chr*1, chr, strand, start, gene_id"):
			strand, ref, alt = r["strand"], r["ref"], r["alt"]
			allele = "{0}/{1}".format(ref, alt)
			tsv.write_line(f, r["chr"], strand, r["start"], allele,
						   r["gene_id"], r["impact"], TransFIC.class_name(r["impact"]),
						   r["sample_freq"], r["sample_prop"], r["prot_changes"], r["xrefs"], null_value="-")
			
	log.info("  Genes ...")

	with tsv.open(os.path.join(base_path, "gene-{0}.tsv.gz".format(group_file_prefix)), "w") as f:
		tsv.write_param(f, "classifier", classifier["id"])
		tsv.write_param(f, "group_id", group_name)
		tsv.write_param(f, "group_short_name", group_short_name)
		tsv.write_param(f, "group_long_name", group_long_name)
		tsv.write_param(f, "projects", ",".join(project_ids))
		tsv.write_param(f, "SAMPLE_TOTAL", sample_total)
		tsv.write_line(f, "GENE_ID", "SAMPLE_FREQ", "SAMPLE_PROP")
		for r in c.execute("SELECT * FROM genes ORDER BY gene_id"):
			tsv.write_line(f, r["gene_id"], r["sample_freq"], r["sample_prop"], null_value="-")

	log.info("  Pathways ...")

	with tsv.open(os.path.join(base_path, "pathway-{0}.tsv.gz".format(group_file_prefix)), "w") as f:
		tsv.write_param(f, "classifier", classifier["id"])
		tsv.write_param(f, "group_id", group_name)
		tsv.write_param(f, "group_short_name", group_short_name)
		tsv.write_param(f, "group_long_name", group_long_name)
		tsv.write_param(f, "projects", ",".join(project_ids))
		tsv.write_param(f, "SAMPLE_TOTAL", sample_total)
		tsv.write_line(f, "PATHWAY_ID", "SAMPLE_FREQ", "SAMPLE_PROP")
		for r in c.execute("SELECT * FROM pathways ORDER BY pathway_id"):
			tsv.write_line(f, r["pathway_id"], r["sample_freq"], r["sample_prop"], null_value="-")
			
	conn.close()

	remove_temp(task, db_path)

task.run()