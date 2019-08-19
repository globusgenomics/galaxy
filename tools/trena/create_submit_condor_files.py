#!/usr/bin/python
import sys, os, glob

input_dir = sys.argv[1]
print input_dir
results_dirs = glob.glob("%s/*.results" % input_dir)

for idir in results_dirs:
    sh_name = "%s.sh" % (idir)
    output_f = "%s/%s.bam" % (idir, os.path.basename(idir))
    input_bams = glob.glob("%s/*/*.bam" % (idir))
    sh_f = open(sh_name, "w")
    sh_f.write("#!/bin/sh\n")
    sh_f.write("PACKAGE_BASE=/mnt/galaxyTools/tools/samtools/1.2; export PACKAGE_BASE; . /mnt/galaxyTools/tools/samtools/1.2/env.sh; samtools merge -f %s %s\n" % (output_f, " ".join(input_bams)))
    sh_f.close()

    out_log = "%s.o" % (idir)
    err_log = "%s.e" % (idir)
    con_log = "%s.condor.log" % (idir)
    submit_f = open("%s.condor.submit" % idir, "w")
    submit_f.write("universe = vanilla\ngetenv = true\nexecutable = %s\n" % (sh_name))
    submit_f.write("output = %s\nerror = %s\nlog = %s\n" % (out_log, err_log, con_log) )
    submit_f.write("notification = NEVER\nrequest_cpus=8\nqueue\n")
