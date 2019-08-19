import os

from wok.task import task

from bgcore.labelfilter import LabelFilter

from intogensm import so
from intogensm.projdb import ProjectDb
from intogensm.projres import ProjectResults
from intogensm.utils import get_project_conf
from intogensm.paths import get_data_gene_filter_path
from intogensm.oncodriveclust import load_cds_len, get_oncodriveclust_configuration, SYN, NON_SYN

@task.foreach()
def oncodriveclust(project):
	log = task.logger
	conf = task.conf

	log.info("--- [{0}] --------------------------------------------".format(project["id"]))

	source_genes = {}
	syn_genes = set()
	selected_genes = set()
	filter_genes = set()
	threshold_genes = set()

	source_samples = {}
	selected_samples = set()
	filter_samples = set()
	threshold_samples = set()

	selected_gene_sample_count = {} # number of samples for each selected gene
	filter_gene_sample_count = {} # number of samples per each gene passing the filter

	# get configuration

	samples_threshold, genes_filter_enabled, genes_filter, filt = get_oncodriveclust_configuration(log, conf, project)

	log.info("Retrieving gene alterations ...")

	projdb = ProjectDb(project["db"])

	data = set()

	for csq in projdb.consequences(join_samples=True):
								   #filters={ProjectDb.CSQ_CTYPES : so.PROTEIN_AFFECTING | so.SYNONYMOUS}):

		is_selected = so.match(csq.ctypes, so.PROTEIN_AFFECTING)
		is_synonymous = so.match(csq.ctypes, so.SYNONYMOUS)

		if csq.gene not in source_genes:
			source_genes[csq.gene] = gene_index = len(source_genes)

		if is_selected:
			selected_genes.add(gene_index)

		if is_synonymous:
			syn_genes.add(gene_index)

		for sample in csq.var.samples:
			if sample.name not in source_samples:
				source_samples[sample.name] = sample_index = len(source_samples)

			if is_selected:
				selected_samples.add(sample_index)
				data.add((csq.gene, sample_index))

	projdb.close()

	log.info("Counting selected, filtered and threshold ...")

	# calculate selected and filter counts

	data2 = set()

	for gene, sample_index in data:
		gene_index = source_genes[gene]
		if gene_index not in selected_gene_sample_count:
			selected_gene_sample_count[gene_index] = 1
		else:
			selected_gene_sample_count[gene_index] += 1

		if filt.valid(gene):
			data2.add((gene_index, sample_index))
			filter_genes.add(gene_index)
			filter_samples.add(sample_index)
			if gene_index not in filter_gene_sample_count:
				filter_gene_sample_count[gene_index] = 1
			else:
				filter_gene_sample_count[gene_index] += 1

	# calculate threshold counts

	for gene_index, sample_index in data2:
		if selected_gene_sample_count[gene_index] >= samples_threshold:
			threshold_genes.add(gene_index)
			threshold_samples.add(sample_index)

	log.info("Counting significant genes ...")

	# significance of q-values

	projdb = ProjectDb(project["db"])
	sig_thresholds = [0.0, 0.001, 0.005] + [i / 100.0 for i in range(1, 11)] + [1.0]
	sig_count = [0] * len(sig_thresholds)
	for gene in projdb.genes():
		if gene.id in source_genes and source_genes[gene.id] in threshold_genes:
			i = 0
			while i < len(sig_thresholds) and gene.fm_qvalue > sig_thresholds[i]:
				i += 1

			for j in range(i, len(sig_count)):
				sig_count[j] += 1

	projdb.close()

	source_genes_count = len(source_genes)
	syn_genes_count = len(syn_genes)
	selected_genes_count = len(selected_genes)
	filter_genes_count = len(filter_genes)
	threshold_genes_count = len(threshold_genes)

	source_samples_count = len(source_samples)
	selected_samples_count = len(selected_samples)
	filter_samples_count = len(filter_samples)
	threshold_samples_count = len(threshold_samples)

	sorted_filter_genes = sorted(filter_genes, reverse=True, key=lambda gi: filter_gene_sample_count[gi])

	qc_data = dict(
			source=dict(
				genes=sorted(source_genes.keys(), key=lambda k: source_genes[k]),
				genes_count=source_genes_count,
				genes_lost_count=max(0, source_genes_count - syn_genes_count - threshold_genes_count),
				samples=sorted(source_samples.keys(), key=lambda k: source_samples[k]),
				samples_count=source_samples_count),
				samples_lost_count=max(0, source_samples_count - threshold_samples_count),
			synonymous=dict(
				genes=sorted(syn_genes),
				genes_count=syn_genes_count,
				ratio=(float(syn_genes_count) / selected_genes_count) if selected_genes_count > 0 else 0),
			selected=dict(
				genes=sorted(selected_genes),
				genes_count=selected_genes_count,
				genes_lost=sorted(set(source_genes.values()) - syn_genes - selected_genes),
				genes_lost_count=max(0, source_genes_count - syn_genes_count - selected_genes_count),
				samples=sorted(selected_samples),
				samples_count=selected_samples_count,
				samples_lost=sorted(set(source_samples.values()) - selected_samples),
				samples_lost_count=max(0, source_samples_count - selected_samples_count)),
			filter=dict(
				genes=sorted_filter_genes,
				genes_count=filter_genes_count,
				genes_lost=sorted(selected_genes - filter_genes),
				genes_lost_count=max(0, selected_genes_count - filter_genes_count),
				genes_sample_count=[filter_gene_sample_count[gene_index] for gene_index in sorted_filter_genes],
				samples=sorted(filter_samples),
				samples_count=filter_samples_count,
				samples_lost=sorted(selected_samples - filter_samples),
				samples_lost_count=max(0, selected_samples_count - filter_samples_count)),
			threshold=dict(
				genes=sorted(threshold_genes),
				genes_count=threshold_genes_count,
				genes_lost=sorted(filter_genes - threshold_genes),
				genes_lost_count=max(0, filter_genes_count - threshold_genes_count),
				samples=sorted(threshold_samples),
				samples_count=threshold_samples_count,
				samples_threshold=samples_threshold,
				samples_lost=sorted(filter_samples - threshold_samples),
				samples_lost_count=max(0, filter_samples_count - threshold_samples_count)),
			results=dict(
				sig_thresholds=sig_thresholds[1:],
				sig_count=sig_count[1:])
			)

	project_results = ProjectResults(project)
	project_results.save_quality_control("oncodriveclust", qc_data)

task.run()