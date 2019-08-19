import os.path
import gzip
import shutil

from wok.task import task

from intogensm.woktmp import make_temp_file, make_temp_dir, remove_temp

from intogensm.projdb import ProjectDb
from intogensm.archive import Archive
from intogensm.pathutils import split_ext
from intogensm.exceptions import InternalError
from intogensm.variants.fileformat.factory import create_variants_parser
from intogensm.paths import *

_SINGLE_CONTAINERS = [".gz", ".bz2"]
_SUPPORTED_CONTAINERS = [".zip", ".tar.gz", ".tar.bz2"] + _SINGLE_CONTAINERS
_SUPPORTED_EXTENSIONS = [".tab", ".vcf", ".maf"]

def archived_files(container_file):
	name, ext = split_ext(container_file)
	path = os.path.dirname(name)
	name = os.path.basename(name)
	if ext in _SUPPORTED_CONTAINERS:
		arc_ext = ext
		arc = Archive(container_file)
		for entry in arc.list():
			name, ext = os.path.splitext(entry.name)
			if arc_ext not in _SINGLE_CONTAINERS and ext.lower() not in _SUPPORTED_EXTENSIONS:
				continue

			path = os.path.dirname(name)
			name = os.path.basename(name)

			tmp = make_temp_dir(task)
			entry.extract(tmp)
			f = open(os.path.join(tmp, entry.name))
			yield (container_file, path, name, ext, f)
			f.close()
			shutil.rmtree(tmp)
	else:
		f = open(container_file, "r")
		yield (None, path, name, ext, f)
		f.close()

@task.foreach()
def scan_files(project):
	log = task.logger
	conf = task.conf

	projects_port, liftover_projects_port = task.ports("projects_out", "liftover_projects")

	project_id = project["id"]
	temp_path = project["temp_path"]
	project_path = project["path"]
	projdb_path = project["db"]
	assembly = project["assembly"]

	log.info("--- [{0}] --------------------------------------------".format(project_id))

	if assembly == "hg18":
		out_port = liftover_projects_port
	elif assembly == "hg19":
		out_port = projects_port
	else:
		raise Exception("Unexpected assembly: {0}".format(assembly))

	#if os.path.exists(projdb_path):
	#	log.warn("Variations database already created, skipping this step.")
	#	out_port.send(project)
	#	return

	if os.path.exists(projdb_path):
		os.remove(projdb_path)

	log.info("Creating variants database ...")

	projdb_tmp_path = make_temp_file(task, suffix=".db")

	log.debug(projdb_tmp_path)

	projdb = ProjectDb(projdb_tmp_path)
	projdb.create()

	data_path = conf["data_path"]

	log.info("Loading genes ...")

	projdb.load_genes(get_data_ensembl_genes_path(conf))

	log.info("Loading pathways ...")

	projdb.load_pathways(
		get_data_kegg_def_path(conf),
		get_data_kegg_ensg_map_path(conf))

	log.info("Parsing variants ...")

	for file in project["files"]:
		if not os.path.isabs(file):
			raise InternalError("Non absolute path found: {0}".format(file))

		if not os.path.exists(file):
			raise Exception("Input file not found: {0}".format(file))

		if not os.path.isfile(file):
			raise Exception("Not a file: {0}".format(file))

		for container_name, path, name, ext, f in archived_files(file):
			fname = os.path.join(path, name + ext)
			if container_name is not None:
				source_name = "{0}:{1}".format(os.path.basename(container_name), fname)
			else:
				source_name = name + ext

			log.info("=> {0} ...".format(source_name))

			sample_id = os.path.basename(name)

			if ext.lower() in _SUPPORTED_EXTENSIONS:
				parser_type = ext[1:]
			else:
				parser_type = "tab"

			parser = create_variants_parser(parser_type, f, source_name, sample_id)

			source_id = projdb.add_source(source_name)

			var_ids = set()
			for var in parser:
				for line_num, text in parser.read_lines():
					projdb.add_source_line(source_id, line_num, text)

				var_id = projdb.add_variant(var, source_id=source_id, line_num=parser.get_line_num())
				var_ids.add(var_id)

			for line_num, text in parser.read_lines():
				projdb.add_source_line(source_id, line_num, text)

			num_variants = len(var_ids)
			log.info("   {0} variants".format(num_variants))

			if num_variants == 0:
				raise Exception("No variants found in source '{}'. "
								"Please check the documentation for the expected input for '{}' format.".format(
								source_name, parser.name))

	projdb.commit()
	projdb.close()

	log.info("Copying variants database ...")

	log.debug("{0} -> {1}".format(projdb_tmp_path, projdb_path))

	shutil.copy(projdb_tmp_path, projdb_path)

	remove_temp(task, projdb_tmp_path)

	out_port.send(project)

task.run()