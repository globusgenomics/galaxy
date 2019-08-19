import re
import os.path

from bgcore import tsv
from bgcore.re import ReContext

from intogensm.projdb import ProjectDb

from intogensm.transfic import TransFIC

from intogensm import so

from wok.task import task

SIMPLE_AA_CHANGE_RE = re.compile(r"^([\*A-Z])(?:/([\*A-Z]))?$", re.IGNORECASE)
COMPLEX_AA_CHANGE_RE = re.compile(r"^([\*\-A-Z]+)(?:/([\*\-A-Z]+))?$", re.IGNORECASE)

IMPACT_CLASSES = set([
	TransFIC.HIGH_IMPACT_CLASS,
	TransFIC.MEDIUM_IMPACT_CLASS,
	TransFIC.LOW_IMPACT_CLASS])

def update_attr(attrs, key, attr_name, value, update, new=lambda value: value):
	if value is not None:
		if key not in attrs:
			attrs[key] = {attr_name : new(value)}
		elif attr_name not in attrs[key]:
			attrs[key][attr_name] = new(value)
		else:
			prev_value = attrs[key][attr_name]
			attrs[key][attr_name] = update(prev_value, value)


@task.foreach()
def fimpact_run(partition):
	log = task.logger
	conf = task.conf

	results_port = task.ports("results")

	project = partition["project"]

	log.info("--- [{0} @ {1}] --------------------------------------------".format(project["id"], partition["index"]))

	log.info("Reading MA scores ...")

	ma_uniprot = {}
	ma_scores = {}

	with open(partition["ma_path"], "r") as f:
		for var_id, uniprot, fi_score in tsv.lines(f, (int, str, float), null_value="-"):
			ma_uniprot[var_id] = uniprot
			ma_scores[var_id] = fi_score

	log.info("Reading VEP results and calculating functional impact ...")

	tfic = TransFIC(data_path=os.path.join(conf["data_path"], "TransFIC"))

	tfi_path = os.path.join(partition["base_path"], "{0:08d}.tfi".format(partition["index"]))
	cf = open(tfi_path, "w")

	aff_gene_attrs = {}

	with open(partition["vep_path"], "r") as f:
		for fields in tsv.lines(f, (int, str, str, str, str, str, str, float, float), null_value="-"):
			(var_id, gene, transcript, ct,
				protein_pos, aa_change, protein,
				sift_score, pph2_score) = fields

			if ct is not None:
				ct = ct.split(",")
			else:
				ct = []

			# Invert sift score
			if sift_score is not None:
				sift_score = 1.0 - sift_score

			ma_score = None

			uniprot = ma_uniprot[var_id] if var_id in ma_uniprot else None

			sift_impact = pph2_impact = ma_impact = None # TransFIC.UNKNOWN_IMPACT_CLASS

			coding_region = so.match(ct, so.CODING_REGION)

			calculate_transfic = True

			ct_type = None
			if so.match(ct, so.NON_SYNONYMOUS):       # missense
				ct_type = TransFIC.CT_NON_SYNONYMOUS
				ma_score = ma_scores[var_id] if var_id in ma_scores else None
			elif so.match(ct, so.STOP):               # stop
				ct_type = TransFIC.CT_STOP
				sift_impact = pph2_impact = ma_impact = TransFIC.HIGH_IMPACT_CLASS
				sift_score = pph2_score = 1.0
				ma_score = 3.5
			elif so.match(ct, so.FRAMESHIFT):         # frameshift
				ct_type = TransFIC.CT_FRAMESHIFT
				sift_impact = pph2_impact = ma_impact = TransFIC.HIGH_IMPACT_CLASS
				sift_score = pph2_score = 1.0
				ma_score = 3.5
			elif so.match(ct, so.SPLICE):             # splice
				ct_type = "splice"
				sift_impact = pph2_impact = ma_impact = TransFIC.HIGH_IMPACT_CLASS if so.match(ct, so.SPLICE_JUNCTION) else TransFIC.UNKNOWN_IMPACT_CLASS
				calculate_transfic = False
			elif so.match(ct, so.SYNONYMOUS):         # synonymous
				ct_type = TransFIC.CT_SYNONYMOUS
				sift_impact = pph2_impact = ma_impact = TransFIC.NONE_IMPACT_CLASS
				sift_score = pph2_score = 0.0
				ma_score = -2
			else:
				sift_impact = pph2_impact = ma_impact = TransFIC.NONE_IMPACT_CLASS
				calculate_transfic = False

			if calculate_transfic:
				(sift_tfic, sift_class,
				 pph2_tfic, pph2_class,
				 ma_tfic, ma_class) = tfic.calculate("gosmf", gene, ct_type, sift_score, pph2_score, ma_score)

				# if the impact was not preassigned get it from the transFIC calculated class
				sift_impact = sift_class if sift_impact is None and sift_class in IMPACT_CLASSES else sift_impact
				pph2_impact = pph2_class if pph2_impact is None and pph2_class in IMPACT_CLASSES else pph2_impact
				ma_impact = ma_class if ma_impact is None and ma_class in IMPACT_CLASSES else ma_impact
			else:
				sift_tfic, sift_class, pph2_tfic, pph2_class, ma_tfic, ma_class = (None, None, None, None, None, None)

			aff_gene = (var_id, gene)

			# update aggregated impact for all the predictors
			update_attr(aff_gene_attrs, aff_gene, "sift_impact", sift_impact, update=TransFIC.higher_impact)
			update_attr(aff_gene_attrs, aff_gene, "pph2_impact", pph2_impact, update=TransFIC.higher_impact)
			update_attr(aff_gene_attrs, aff_gene, "ma_impact", ma_impact, update=TransFIC.higher_impact)

			# update whether the affected gene is a coding region or not
			update_attr(aff_gene_attrs, aff_gene, "coding_region", coding_region,
						update=lambda prev_value, value: prev_value or value)

			# aggregate protein changes per affected_gene
			# try to follow the convention http://www.hgvs.org/mutnomen/recs-prot.html
			prot_change = None
			if ct_type == TransFIC.CT_FRAMESHIFT:
				if protein_pos is None:
					prot_change = "fs"
				else:
					prot_change = "fs {0}".format(protein_pos)
				#log.debug("FRAMESHIFT: gene={}, protein={}, pos={}, change={}".format(gene, protein, protein_pos, aa_change))
			elif ct_type == "splice":
				prot_change = "r.spl?"
				#log.debug("SPLICE: gene={}, protein={}, pos={}, change={}".format(gene, protein, protein_pos, aa_change))
			elif protein_pos is not None and aa_change is not None:
				rc = ReContext()
				if rc.match(SIMPLE_AA_CHANGE_RE, aa_change):
					prot_change = "{ref}{pos}{alt}".format(pos=protein_pos, ref=rc.group(1), alt=rc.group(2) or "=")
				elif rc.match(COMPLEX_AA_CHANGE_RE, aa_change):
					prot_change = "{0} {1}".format(aa_change, protein_pos)
				else:
					log.warn("Unmatched aa change: gene={}, protein={}, pos={}, change={}, ct=[{}]".format(
													gene, protein, protein_pos, aa_change, ", ".join(ct)))

			if prot_change is not None:
				update_attr(aff_gene_attrs, aff_gene, "prot_changes", prot_change,
								 new=lambda value: set([value]),
								 update=lambda prev_value, value: prev_value | set([value]))

			impact = ma_impact or pph2_impact or sift_impact or TransFIC.UNKNOWN_IMPACT_CLASS

			tsv.write_line(cf, var_id, transcript, uniprot,
						   sift_score, sift_tfic, sift_class,
						   pph2_score, pph2_tfic, pph2_class,
						   ma_score, ma_tfic, ma_class,
						   impact, null_value="-")

	cf.close()

	log.info("Saving variant impacts ...")

	gfi_path = os.path.join(partition["base_path"], "{0:08d}.gfi".format(partition["index"]))
	vf = open(gfi_path, "w")
	for aff_gene, attrs in aff_gene_attrs.items():
		var_id, gene = aff_gene
		# get the impact by trust priority: ma, pph2, sift
		impact = attrs.get("ma_impact") or attrs.get("pph2_impact") or attrs.get("sift_impact") or TransFIC.UNKNOWN_IMPACT_CLASS
		coding_region = attrs.get("coding_region", False)
		coding_region = 1 if coding_region else 0
		prot_changes = attrs.get("prot_changes")
		prot_changes = ",".join(prot_changes) if prot_changes is not None else None

		tsv.write_line(vf, var_id, gene, impact, coding_region, prot_changes, null_value="-")
	vf.close()

	# Send results to the next module
	partition["tfi_path"] = tfi_path
	partition["gfi_path"] = gfi_path
	results_port.send(partition)

task.run()