#!/usr/bin/python
import time

from optparse import OptionParser
args=None
parser = OptionParser()

parser.add_option("--wait", dest="wait", help="wait")

options, args = parser.parse_args(args)

wait = int(options.wait)

time.sleep(wait)