import os

from bgcore.labelfilter import LabelFilter

from intogensm import so
from intogensm.projdb import ProjectDb
from intogensm.utils import get_project_conf
from intogensm.paths import get_data_gene_filter_path
from intogensm.constants.oncodrivefm import *

def get_threshold(log, conf, project, key, default, num_samples):
	threshold = get_project_conf(conf, project, key, default)
	try:
		threshold = int(threshold)
	except:
		if isinstance(threshold, basestring) and threshold.endswith("%"):
			try:
				threshold = max(num_samples * int(threshold[0:-1]) / 100, 2)
			except:
				log.warn("Wrong threshold value for '{0}'. Using the default value '{1}' instead.".format(threshold[0:-1], default))
				threshold = default
		else:
			log.warn("Wrong threshold value for '{0}'. Using the default value '{1}' instead.".format(threshold, default))
			threshold = default
	return threshold

def load_expression_filter(log, conf, project):

	default_filter = get_data_gene_filter_path(conf)
	genes_filter_enabled = get_project_conf(conf, project,
											"oncodrivefm.genes.filter_enabled", ONCODRIVEFM_GENES_FILTER_ENABLED)

	genes_filter = get_project_conf(conf, project, "oncodrivefm.genes.filter", default_filter)
	if genes_filter is None: # user can assign a null
		genes_filter_enabled = False
		genes_filter = default_filter

	filt = LabelFilter()

	if genes_filter_enabled:
		log.info("Loading expression filter ...")
		log.debug("> {0}".format(genes_filter))
		filt.load(genes_filter)

	return genes_filter_enabled, genes_filter, filt

def get_oncodrivefm_configuration(log, conf, project, num_samples):
	log.info("OncodriveFM configuration:")

	num_cores = get_project_conf(conf, project, "oncodrivefm.num_cores", ONCODRIVEFM_NUM_CORES)
	estimator = get_project_conf(conf, project, "oncodrivefm.estimator", ONCODRIVEFM_ESTIMATOR)

	log.info("  num_cores = {0}".format(num_cores))
	log.info("  estimator = {0}".format(estimator))

	log.info("  Genes:")

	genes_num_samplings = get_project_conf(conf, project,
										   "oncodrivefm.genes.num_samplings", ONCODRIVEFM_GENES_NUM_SAMPLINGS)

	genes_threshold = get_threshold(log, conf, project,
											"oncodrivefm.genes.threshold", ONCODRIVEFM_GENES_THRESHOLD, num_samples)

	genes_filter_enabled, genes_filter, filt = load_expression_filter(log, conf, project)

	log.info("    num_samplings = {0}".format(genes_num_samplings))
	log.info("    threshold = {0}".format(genes_threshold))
	log.info("    filter enabled = {0}".format(genes_filter_enabled))
	log.info("    filter = {0}".format(os.path.basename(genes_filter) if genes_filter is not None else None))

	log.info("  Pathways:")

	pathways_num_samplings = get_project_conf(conf, project,
									"oncodrivefm.pathways.num_samplings", ONCODRIVEFM_PATHWAYS_NUM_SAMPLINGS)
	pathways_threshold = get_threshold(log, conf, project,
									"oncodrivefm.pathways.threshold", ONCODRIVEFM_PATHWAYS_THRESHOLD, num_samples)

	log.info("    num_samplings = {0}".format(pathways_num_samplings))
	log.info("    threshold = {0}".format(pathways_threshold))

	return (num_cores, estimator,
			genes_num_samplings, genes_threshold, genes_filter_enabled, genes_filter, filt,
			pathways_num_samplings, pathways_threshold)

_AGGREGATE = (max, max, min,
			  max, max, min,
			  max, max, min)

def retrieve_data(projdb):
	data = {}
	for csq in projdb.consequences(join_samples=True, join_ctypes=False,
								   filters={ProjectDb.CSQ_CTYPES : so.ONCODRIVEFM}):

		var = csq.var
		for sample in var.samples:
			key = (sample.id, csq.gene)
			if key not in data:
				data[key] = (csq.sift_score, csq.sift_tfic, csq.sift_tfic_class,
							 csq.pph2_score, csq.pph2_tfic, csq.pph2_tfic_class,
							 csq.ma_score, csq.ma_tfic, csq.ma_tfic_class)
			else:
				a = data[key]

				b = (csq.sift_score, csq.sift_tfic, csq.sift_tfic_class,
					 csq.pph2_score, csq.pph2_tfic, csq.pph2_tfic_class,
					 csq.ma_score, csq.ma_tfic, csq.ma_tfic_class)

				c = [0.0] * len(a)

				for i in range(len(a)):
					c[i] = _AGGREGATE[i](a[i], b[i])

				data[key] = tuple(c)

	return data
