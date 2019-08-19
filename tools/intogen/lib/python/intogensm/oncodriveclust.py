import os

from bgcore import tsv
from bgcore.labelfilter import LabelFilter

from intogensm import so
from intogensm.projdb import ProjectDb
from intogensm.paths import get_data_ensembl_gene_transcripts_path
from intogensm.utils import get_project_conf
from intogensm.paths import get_data_gene_filter_path
from intogensm.constants.oncodriveclust import *

NON_SYN = 0
SYN = 1

def load_cds_len(conf):
	cds_len = {}
	with tsv.open(get_data_ensembl_gene_transcripts_path(conf), "r") as f:
		for gene, transcript, transcript_len in tsv.lines(f, (str, str, int), header=True):
			cds_len[transcript] = transcript_len
	return cds_len

def retrieve_data(projdb, cds_len):
	data = {}

	for csq in projdb.consequences(join_samples=True,
								   filters={ProjectDb.CSQ_CTYPES : so.PROTEIN_AFFECTING | so.SYNONYMOUS}):

		if csq.transcript not in cds_len:
			continue

		transcript_len = cds_len[csq.transcript]

		if so.match(csq.ctypes, so.PROTEIN_AFFECTING):
			cls = NON_SYN
		elif so.match(csq.ctypes, so.SYNONYMOUS):
			cls = SYN
		else:
			continue

		for sample in csq.var.samples:
			key = (cls, csq.gene, sample.name)
			if key not in data:
				data[key] = (csq.transcript, transcript_len, csq.protein_pos)
			else:
				transcript, tlen, protein_pos = data[key]
				if transcript_len > tlen:
					data[key] = (csq.transcript, transcript_len, csq.protein_pos)

	return data

def get_oncodriveclust_configuration(log, conf, project):
	log.info("OncodriveCLUST configuration:")

	mutations_threshold = get_project_conf(conf, project, "oncodriveclust.mutations_threshold", ONCODRIVECLUST_MUTATIONS_THRESHOLD)

	default_filter = get_data_gene_filter_path(conf)
	genes_filter_enabled = get_project_conf(conf, project, "oncodriveclust.genes_filter_enabled", ONCODRIVECLUST_GENES_FILTER_ENABLED)
	genes_filter = get_project_conf(conf, project, "oncodriveclust.genes_filter", default_filter)
	if genes_filter is None: # user can assign a null
		genes_filter_enabled = False
		genes_filter = default_filter

	log.info("  mutations_threshold = {0}".format(mutations_threshold))
	log.info("  genes_filter_enabled = {0}".format(genes_filter_enabled))
	log.info("  genes_filter = {0}".format(os.path.basename(genes_filter)))

	filt = LabelFilter()

	if genes_filter_enabled:
		log.info("Loading expression filter ...")
		log.debug("> {0}".format(genes_filter))
		filt.load(genes_filter)

	return (mutations_threshold, genes_filter_enabled, genes_filter, filt)