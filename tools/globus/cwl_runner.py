#!/usr/bin/python

import subprocess
import time
import os
import sys
from shutil import copyfile
from optparse import OptionParser


TMP_DIR = "/scratch/galaxy/tmp"

WORKING_DIR = "/ephemeral/0/tmp/"

args=None
parser = OptionParser()

parser.add_option("--cwl-file", dest="cwl_file", help="CWL File")
parser.add_option("--inputs", dest="inputs", help="Inputs JSON file")
parser.add_option("--output", dest="output", help="Output")

options, args = parser.parse_args(args)

#cwl_file = options.cwl_file
cwl_file = "/home/galaxy/sbg_dockstore_tools/topmed-workflows/alignment/topmed-alignment.cwl"
print cwl_file
inputs = options.inputs
print inputs
output = options.output
print output

outdir = os.path.join(TMP_DIR, str(int(time.time()*10000)))
if not os.path.exists(outdir):
    os.makedirs(outdir)

command = 'cwl-runner --no-match-user --outdir {0} --tmpdir-prefix {3} --tmp-outdir-prefix {3} {1} {2}'.format(outdir, cwl_file, inputs, WORKING_DIR)

subprocess.call(command, shell=True)


output_list = os.listdir(outdir)
print output_list
if output_list == []:
    sys.exit("No output file.")
else:
    output_file = output_list[0]

source_path = os.path.join(outdir, output_file)

if os.path.exists(output):
    os.remove(output)

copyfile(source_path, output)




"""
{
  "output_name": "output1",
  "output": {
    "path": "/home/ubuntu/test/output1",
    "class": "File"
  },
    "bwa_index": {
        "class": "File",
        "path": "ark:/57799/b91414"
    },
    "dbsnp": {
        "class": "File",
        "path": "ark:/57799/b9wb0k"
    },
    "input_file": {
        "class": "File",
        "path": "ark:/57799/b9rm50"
    },
    "reference_genome": {
        "class": "File",
        "path": "ark:/57799/b9mt4f"
    }
}

hs38DH.fa.tar    hs38DH_fa_tar_bdbag.zip    ark:/57799/b91414
Homo_sapiens_assembly38.dbsnp138.vcf.gz    Homo_sapiens_assembly38_dbsnp138_bdbag.zip    ark:/57799/b9wb0k
NWD176325.0005.recab.cram    NWD176325_0005_recab_cram_bdbag.zip    ark:/57799/b9rm50
hs38DH.fa    hs38DH_fa_bdbag.zip    ark:/57799/b9mt4f

{"bwa_index": {"class": "File", "path": "ark:/57799/b91414" }, "dbsnp": {"class": "File", "path": "ark:/57799/b9wb0k"}, "input_file": {"class": "File", "path": "ark:/57799/b9rm50"}, "reference_genome": {"class": "File", "path": "ark:/57799/b9mt4f"}}

"""