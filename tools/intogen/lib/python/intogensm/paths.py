import os

def get_results_path(conf, *path):
	return os.path.join(conf["results_path"], conf["workspace"], *path)

def get_projects_path(conf, *path):
	return get_results_path(conf, "projects", *path)

def get_project_path(conf, project_id):
	return get_projects_path(conf, project_id)

def get_temp_path(conf, project_id=None):
	if project_id is not None:
		return os.path.join(conf["temp_path"], conf["workspace"], project_id)
	else:
		return os.path.join(conf["temp_path"], conf["workspace"])

# Combination paths ------------------------------------------------------------

def get_combination_path(conf, *path):
	return os.path.join(conf["results_path"], conf["workspace"], "combination", *path)

def create_combination_folders(conf):
	combination_path = get_combination_path(conf)
	if not os.path.exists(combination_path):
		os.makedirs(combination_path)
	for name in ["recurrences", "oncodrivefm", "oncodriveclust"]:
		path = os.path.join(combination_path, name)
		if not os.path.exists(path):
			os.makedirs(path)
	return combination_path

# Data paths -------------------------------------------------------------------

def get_data_path(conf, *path):
	return os.path.join(conf["data_path"], *path)

def get_data_ensembl_path(conf, *path):
	return get_data_path(conf, "Ensembl", *path)

def get_data_ensembl_genes_path(conf):
	return get_data_ensembl_path(conf, "genes.tsv")

def get_data_ensembl_gene_transcripts_path(conf):
	return get_data_ensembl_path(conf, "gene_transcripts.tsv")

def get_data_kegg_path(conf, *path):
	return get_data_path(conf, "KEGG", *path)

def get_data_kegg_def_path(conf):
	return get_data_kegg_path(conf, "kegg.tsv")

def get_data_kegg_ensg_map_path(conf):
	return get_data_kegg_path(conf, "ensg_kegg.tsv")

def get_data_gene_filter_path(conf):
	return get_data_path(conf, "gene-filter.txt")

# website -----------------------------------------------------------------------

def get_website_path(project_path, *path):
	return os.path.join(project_path, "website", *path)

def get_website_results_path(project_path):
	return get_website_path(project_path, "results")

# Examples -----------------------------------------------------------------------

def get_examples_path(conf, *path):
	return os.path.join(conf["examples_path"], *path)

# Quality control ----------------------------------------------------------------

def get_quality_variants_path(project_path):
	return os.path.join(project_path, "quality_variants.json")

def get_quality_oncodrivefm_path(project_path):
	return os.path.join(project_path, "quality_oncodrivefm.json")

def get_quality_oncodriveclust_path(project_path):
	return os.path.join(project_path, "quality_oncodriveclust.json")