import os.path

from math import log

class TransFIC(object):
	# Consequence group
	CT_SYNONYMOUS = 1
	CT_STOP = 2
	CT_FRAMESHIFT = 3
	CT_NON_SYNONYMOUS = 4

	# Impact class
	HIGH_IMPACT_CLASS = 1
	MEDIUM_IMPACT_CLASS = 2
	LOW_IMPACT_CLASS = 3
	UNKNOWN_IMPACT_CLASS = 4
	NONE_IMPACT_CLASS = 5

	CLASS_NAME = {
		HIGH_IMPACT_CLASS : "high",
		MEDIUM_IMPACT_CLASS : "medium",
		LOW_IMPACT_CLASS : "low",
		UNKNOWN_IMPACT_CLASS : "unknown",
		NONE_IMPACT_CLASS : "none"
	}

	@staticmethod
	def higher_impact(a, b):
		return min(a, b)

	@staticmethod
	def class_name(cls):
		if cls is None:
			return None
		return TransFIC.CLASS_NAME[cls]

	def __init__(self, data_path):
		self.data_path = data_path

		self.__partitions = {}

	def __load_partition(self, name):
		if name not in self.__partitions:
			self.__partitions[name] = {}

		part = self.__partitions[name]

		with open(os.path.join(self.data_path, name + ".genes"), "r") as f:
			f.readline() # discard header
			for line in f:
				fields = line.rstrip("\n").split("\t")
				gene = fields[0]
				#sift_mean, sift_sd, pph2_mean, pph2_sd, ma_mean, ma_sd = [float(x) for x in fields[1:7]]
				part[gene] = tuple([float(x) for x in fields[1:7]])

		return part

	def __get_partition(self, name):
		if name not in self.__partitions:
			self.__load_partition(name)

		return self.__partitions[name]

	def __calculate(self, ct_type, score, mean, sd, class_thresholds, is_ma=False):
		if score is None:
			return (None, None)

		if ct_type == self.CT_SYNONYMOUS:
			tfic = mean - 2.0 * sd
		elif ct_type in [self.CT_FRAMESHIFT, self.CT_STOP]:
			tfic = mean + 2.0 * sd
		elif ct_type == self.CT_NON_SYNONYMOUS:
			if not is_ma:
				if score == 0.0:
					score = 0.001
				elif score == 1.0:
					score = 0.999

				tfic = log(score / (1 - score)) - (mean / sd)
			else:
				tfic = (score - mean) / sd

		if tfic < class_thresholds[0]:
			tfic_class = self.LOW_IMPACT_CLASS
		elif tfic < class_thresholds[1]:
			tfic_class = self.MEDIUM_IMPACT_CLASS
		else:
			tfic_class = self.HIGH_IMPACT_CLASS

		return tfic, tfic_class

	def calculate(self, part_name, gene, ct, sift_score, pph2_score, ma_score):
		part = self.__get_partition(part_name)

		if gene not in part:
			return (None, None, None, None, None, None)

		sift_mean, sift_sd, pph2_mean, pph2_sd, ma_mean, ma_sd = part[gene]

		sift_tfic, sift_class = self.__calculate(ct, sift_score, sift_mean, sift_sd, (0.0, 2.0))
		pph2_tfic, pph2_class = self.__calculate(ct, pph2_score, pph2_mean, pph2_sd, (0.0, 1.5))
		ma_tfic, ma_class = self.__calculate(ct, ma_score, ma_mean, ma_sd, (1.0, 3.0), is_ma=True)

		return (sift_tfic, sift_class, pph2_tfic, pph2_class, ma_tfic, ma_class)