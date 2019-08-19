"""Wrapper script providing conversion from BAM to fastq, handling paired ends.

Requires:
  Picard (http://picard.sourceforge.net/)
    The SamToFastq.jar file needs to be linked from this directory or available
    in a standard directory like /usr/share/java/picard.
  pysam (http://code.google.com/p/pysam/)
"""
import os
import sys
import subprocess

import pysam

def main(in_bam, out_fastq, out_id, extra_file_dir):
    """out_fastq2 = check_for_paired(in_bam, out_id, extra_file_dir)"""
    out_fastq2 = ""
    picard_jar = find_picard_jar("SamToFastq")
    opts = [("INPUT", in_bam), ("FASTQ", out_fastq),
            ("QUIET", "true"), ("VERBOSITY", "WARNING"),
            ("VALIDATION_STRINGENCY", "SILENT"),
            ("TMP_DIR", "/glusterfs/galaxy-data/tmp")]
    if out_fastq2:
        opts.append(("SECOND_END_FASTQ", out_fastq2))
    opts = ["%s=%s" % (x, y) for x, y in opts]
    cl = ["java", "-jar", picard_jar] + opts
    subprocess.check_call(cl)

def find_picard_jar(name):
    test_dirs = [os.path.dirname(__file__), "/usr/share/java/picard", "/nfs/software/picard"]
    for d in test_dirs:
        f = os.path.join(d, "%s.jar" % name)
        if os.path.exists(f):
            return f
    raise ValueError("Could not find %s in %s" % (name, test_dirs))

def check_for_paired(in_bam, out_id, extra_file_dir):
    if is_paired(in_bam):
        return os.path.join(extra_file_dir, "%s_%s_%s_%s_%s" %
                            ('primary', out_id, 'pair2', 'visible', 'fastqsanger'))
    else:
        return None

def is_paired(in_bam):
    samfile = pysam.Samfile(in_bam, "rb")
    read = samfile.fetch(until_eof=True).next()
    samfile.close()
    return read.is_paired

if __name__ == "__main__":
    main(*sys.argv[1:])
