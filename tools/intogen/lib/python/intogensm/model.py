from collections import namedtuple

# Sample

SampleBase = namedtuple("Sample", "id, name")

class Sample(SampleBase):
	def __new__(cls, id=None, name=None):
		return super(Sample, cls).__new__(cls, id, name)

# Variant

VariantBase = namedtuple("Variant", """id, type, chr, start, ref, alt, strand, samples, xrefs""")

class Variant(VariantBase):
	SUBST = "SUBST"
	INDEL = "INDEL"
	COMPLEX = "COMPLEX"

	def __new__(cls, id=None, type=None, chr=None, start=None, ref=None, alt=None, strand=None, samples=None, xrefs=None):
		return super(Variant, cls).__new__(cls, id, type, chr, start, ref, alt, strand, samples, xrefs)

# Consequence

ConsequenceBase = namedtuple("ConsequenceBase", """id, var, transcript, gene, ctypes,
								extid, protein, uniprot, protein_pos, aa_change,
								sift_score, sift_tfic, sift_tfic_class,
								pph2_score, pph2_tfic, pph2_tfic_class,
								ma_score, ma_tfic, ma_tfic_class,
								impact, rec""")

class Consequence(ConsequenceBase):

	def __new__(cls, id=None, var=None, transcript=None, gene=None, ctypes=None,
				 extid=None, protein=None, uniprot=None, protein_pos=None, aa_change=None,
				 sift_score=None, sift_tfic=None, sift_tfic_class=None,
				 pph2_score=None, pph2_tfic=None, pph2_tfic_class=None,
				 ma_score=None, ma_tfic=None, ma_tfic_class=None, impact=None, rec=None):
		return super(Consequence, cls).__new__(cls,
			id, var, transcript, gene, ctypes,
			extid, protein, uniprot, protein_pos, aa_change,
			sift_score, sift_tfic, sift_tfic_class,
			pph2_score, pph2_tfic, pph2_tfic_class,
			ma_score, ma_tfic, ma_tfic_class, impact, rec)


# AffectedGene

AffectedGeneBase = namedtuple("AffectedGene", """id, var, gene_id, impact, coding_region, prot_changes, rec""")

class AffectedGene(AffectedGeneBase):
	def __new__(cls, id=None, var=None, gene_id=None, impact=None, coding_region=None, prot_changes=None, rec=None):
		return super(AffectedGene, cls).__new__(cls, id, var, gene_id, impact, coding_region, prot_changes, rec)


AffectedGeneRec = namedtuple("AffectedGeneRec", "afg_id, sample_freq, sample_prop")

# Gene

GeneBase = namedtuple('Gene', 'id, symbol, xrefs, fm_pvalue, fm_qvalue, fm_exc_cause,'
						' clust_coords, clust_zscore, clust_pvalue, clust_qvalue, clust_exc_cause, rec')

class Gene(GeneBase):
	def __new__(cls, id=None, symbol=None, xrefs=None, fm_pvalue=None, fm_qvalue=None, fm_exc_cause=None,
				clust_coords=None, clust_zscore=None, clust_pvalue=None, clust_qvalue=None, clust_exc_cause=None, rec=None):
		return super(Gene, cls).__new__(cls, id, symbol, xrefs, fm_pvalue, fm_qvalue, fm_exc_cause,
										clust_coords, clust_zscore, clust_pvalue, clust_qvalue, clust_exc_cause,
										rec)

GeneRec = namedtuple("GeneRec", "gene_id, sample_freq, sample_prop")

# Pathway

PathwayBase = namedtuple('Gene', 'id, gene_count, fm_zscore, fm_pvalue, fm_qvalue, rec')

class Pathway(PathwayBase):
	def __new__(cls, id=None, gene_count=None, fm_zscore=None, fm_pvalue=None, fm_qvalue=None, rec=None):
		return super(Pathway, cls).__new__(cls, id, gene_count, fm_zscore, fm_pvalue, fm_qvalue, rec)

PathwayRec = namedtuple("PathwayRec", "pathway_id, sample_freq, sample_prop, gene_freq, gene_prop")