###############################################################################
#                                                                             #
#    Copyright 2009-2011, Universitat Pompeu Fabra                            #
#                                                                             #
#    This file is part of the OncodriveFM interfaces package.                 #
#    Here you find all the functions available in the package.                #
#    This functions can be directly called or used via the different          #
#    interfaces found in the same package.                                    #
###############################################################################

import os
import sys
import re

from wok.config import ConfigBuilder

"""
#FIXME Moved to intogensm.utils
def norm_id(name):
	return re.sub(r"[^a-zA-Z0-9 -\.]", "_", name)

#FIXME Moved to intogensm.utils
def split_ext(file_name):
	name, ext = os.path.splitext(file_name)
	if ext.lower() in [".gz", ".bz2"]:
		name2, ext2 = os.path.splitext(name)
		if ext2.lower() == "tar":
			name = name2
			ext = ext2 + ext
	return (name, ext)
"""

def project_analysis(files, assembly="hg19", conf_builder=None):

	cb = ConfigBuilder(conf_builder)
	cb.add_value("files", files)
	cb.add_value("assembly", assembly)

	conf = cb.get_conf()

	flow_uri = "intogen-mutations:project"

	return cb, flow_uri

def projects_analysis(scan_paths, includes, excludes, conf_builder=None):

	cb = ConfigBuilder(conf_builder)
	cb.add_value("scan_paths", scan_paths)
	cb.add_value("includes", includes)
	cb.add_value("excludes", excludes)

	conf = cb.get_conf()

	flow_uri = "intogen-mutations:project"

	return cb, flow_uri


