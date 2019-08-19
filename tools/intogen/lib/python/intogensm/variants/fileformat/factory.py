from tab import TabParser
from vcf import VcfParser
from maf import MafParser

__PARSERS = {
	"tab" : TabParser,
	"vcf" : VcfParser,
	"maf" : MafParser
}

def create_variants_parser(type, f, fname, default_sample_id):
	if type not in __PARSERS:
		raise Exception("Unknown mutations parser type: {0}".format(type))

	return __PARSERS[type](f, fname, default_sample_id)