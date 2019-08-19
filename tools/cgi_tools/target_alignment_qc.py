#!/usr/bin/python
import os, sys

# calculate the TMB given a filtered VCF and a targeted BED file
target_bed = sys.argv[1]
base_count = 0
region_length = {}
base2target = {}
fh = open(target_bed, "r")
for line in fh:
    line = line.rstrip("\n")
    values = line.split("\t")
    base_count += int(values[2]) - int(values[1]) + 1
    #region_count += 1
    #region_length[values[4]] = int(values[2]) - int(values[1]) + 1
    for x in range(int(values[1]), int(values[2])+1):
        base2target[values[0] + "_" + str(x)] = values[4]
        if values[4] in region_length:
            region_length[values[4]] += 1
        else:
            region_length[values[4]] = 1
        #print values[0] + "_" + str(x) + "\t" + values[4]
fh.close()

region_count = len(region_length)
#for region in region_length:
#   print "%s\t%s" % (region, region_length[region])

#print base2target
per_base_coverage = sys.argv[2]
per_target_coverage = sys.argv[3]
cov_threshold = 60

base_fh = open(per_base_coverage, "r")
base_covs_per_chr = {}
base_covs_per_region = {}
below_threshold_targets = {}
above_threshold_targets = {}

total_base_cov = 0
header_base = base_fh.readline()
histogram_cov_count = {}
histogram_cov_count_list = []
for line in base_fh:
    line = line.rstrip("\n")
    values = line.split("\t")
    location = values[0] + "_" + values[1]
    target_name = base2target[location]
    histogram_cov_count_list.append(int(values[3]))
    if int(values[3]) in histogram_cov_count:
        histogram_cov_count[int(values[3])] += 1
    else:
        histogram_cov_count[int(values[3])] = 1

    if values[0] in base_covs_per_chr:
        base_covs_per_chr[values[0]] += int(values[3])
    else:
        base_covs_per_chr[values[0]] = int(values[3])
    if target_name in base_covs_per_region:
        base_covs_per_region[target_name] += int(values[3])
    else:
        base_covs_per_region[target_name] = int(values[3])

    if int(values[3]) < cov_threshold:
        if target_name in below_threshold_targets:
            below_threshold_targets[target_name] += 1
        else:
            below_threshold_targets[target_name] = 1
    else:
        if target_name in above_threshold_targets:
            above_threshold_targets[target_name] += 1
        else:
            above_threshold_targets[target_name] = 1

    total_base_cov += int(values[3])

#print above_threshold_targets
below_threshold_size = 0
below_threshold_bases = 0
percent_target_coverage = {}
for region in below_threshold_targets:
    below_threshold_size += region_length[region]
    below_threshold_bases += below_threshold_targets[region]
    percent_target_coverage[region] = float(below_threshold_targets[region])/float(region_length[region])

for region in above_threshold_targets:
    percent_target_coverage[region] = float(above_threshold_targets[region])/float(region_length[region])

#for region in percent_target_coverage:
#    print "%s\t%s" % (region, percent_target_coverage[region])

count_above_percent = {}
for i in [0,10,20,30,40,50,60,70,80,90,100]:
    #print i
    for region in percent_target_coverage:
        #print "%s\t%s" % (region, percent_target_coverage[region])
        if percent_target_coverage[region]  *100 >= i:
            if str(i) in count_above_percent:
                count_above_percent[str(i)] += 1
            else:
                count_above_percent[str(i)] = 1

#print count_above_percent
below_threshold_regions = len(below_threshold_targets)
average_coverage = float(total_base_cov) / float(base_count)

### fraction of target covered
lis = range(0,102)
it = iter(lis)
fractions_of_targets_covered = []
for x in it:
    y = next(it)
    record = [x,y]
    count = 0
    for region in percent_target_coverage:
        if percent_target_coverage[region]*100 >= x and percent_target_coverage[region]*100 < y+1:
            count += 1
    record.append(count)
    fractions_of_targets_covered.append(record)

################
print "Summary"
print "#########"
print "Number target regions\t%s" % str(region_count)
print "Total length of target regions\t%s" % str(base_count)
print "Average coverage\t%s" % average_coverage
print "Number of target regions with coverage below 60\t%s" % str(below_threshold_regions)
print "Total length of regions of targets with coverage below 60\t%s" % str(below_threshold_size)
print "Total bases in regions with coverage below 60\t%s" % str(below_threshold_bases)
print "\n"
print "Fractions of targets with coverage at least 60"
print "################################################"
print "Number of targeted regions for which\tCount\tPercentage"
print ">100 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["100"]), str(100*(float(count_above_percent["100"])/float(region_count))))
print ">90 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["90"]), str(100*(float(count_above_percent["90"])/float(region_count))))
print ">80 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["80"]), str(100*(float(count_above_percent["80"])/float(region_count)))) 
print ">70 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["70"]), str(100*(float(count_above_percent["70"])/float(region_count))))
print ">60 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["60"]), str(100*(float(count_above_percent["60"])/float(region_count))))
print ">50 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["50"]), str(100*(float(count_above_percent["50"])/float(region_count))))
print ">40 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["40"]), str(100*(float(count_above_percent["40"])/float(region_count))))
print ">30 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["30"]), str(100*(float(count_above_percent["30"])/float(region_count))))
print ">20 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["20"]), str(100*(float(count_above_percent["20"])/float(region_count))))
print ">10 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["10"]), str(100*(float(count_above_percent["10"])/float(region_count))))
print ">0 percent of the targeted region has coverage at least %s\t%s\t%s" % (str(cov_threshold), str(count_above_percent["0"]), str(100*(float(count_above_percent["0"])/float(region_count))))
print ""
print "fraction of target covered\tcount: Fractions of targets with coverage at least 60"
for record in fractions_of_targets_covered:
    print "%s-%s\t%s" % (record[0], record[1], record[2])
print ""
#print "Coverage of target region positions"
#print "####################################"
#print "Coverage\tnumber of positions: Coverage distribution"
#hist_bins = 1000
#slices = np.linspace(sorted(histogram_cov_count_list)[0], sorted(histogram_cov_count_list)[-1], bins+1, True).astype(np.int)
#counts = np.diff(slices)
#histogram_cov_count_list
