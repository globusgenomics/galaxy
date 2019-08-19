#!/usr/bin/python
import sys, os

coverage_file = sys.argv[1]
vcf_file = sys.argv[2]
depth_cut = sys.argv[3]
sample_name = sys.argv[4]
coverage = {}

poor_cov_t = 60
low_cov_t = 200

lcr = []
pcr = []

c_fh = open(coverage_file, "r")
header = c_fh.readline()
for line in c_fh:
    line = line.rstrip("\n")
    values = line.split("\t")

    if values[0] in coverage:
        coverage[values[0]].append([values[1], values[2], values[6]])
    else:
        coverage[values[0]] = [[values[1], values[2], values[6]]]

    if float(values[6]) <= float(poor_cov_t):
        pcr.append("%s:%s-%s" % (values[0], values[1], values[2]))
    elif float(values[6]) <= float(low_cov_t):
        lcr.append("%s:%s-%s" % (values[0], values[1], values[2]))

#print coverage
c_fh.close()

v_fh = open(vcf_file, "r")
for line in v_fh:
    line = line.rstrip("\n")
    if "Sample1" in line:
        line = line.replace("Sample1", sample_name)
        print line
    elif line.startswith("#") and not line.startswith("##source"):
        print line
    elif line.startswith("##source=Mutect2"):
        print line
        print "##PCR_threshold_value=<60>"
        print "##LCR_threshold_value=<200>"
        print "##PCR=%s" % (",".join(pcr))
        print "##LCR=%s" % (",".join(lcr))
        print "##INFO=<ID=TargetCov,Number=1,Type=String,Description=\"Pass or fail depth coverage in region given Threshold %s\"" % depth_cut
        print "##INFO=<ID=TargetReg,Number=1,Type=String,Description=\"Region from target covered by variant.\""
    elif "##contig" in line and "_" in line:
        continue
    else:
        values = line.split("\t")
        for region in coverage[values[0]]:
            if region[0] <= values[1] <= region[1]:
                result = None
                if float(region[2]) >= float(depth_cut):
                    result = "PASS"
                else:
                    result = "FAIL"
                values[7] += ";TargetCov=%s;TargetReg=%s:%s-%s" % (result, values[0], region[0], region[1])
                print "%s" % ("\t".join(values))
