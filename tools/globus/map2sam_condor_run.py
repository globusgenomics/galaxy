#!/usr/bin/env python

"""
Usage:

map2sam_condor_run.py myinput.cga $reference $output.cga $output_dir

"""


import json
import os
import sys

from condor_run import CondorJob, CondorQueueItem


if __name__ == '__main__':
    input_path = sys.argv[1]
    reference_file = sys.argv[2]
    output_json_path = sys.argv[3]
    output_dir = sys.argv[4]

    selected_lanes = "01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16"
    selected_parts = "01,02,03,04,05,06,07,08,09,10,11,12"

    if len(sys.argv) > 5:
        selected_lanes = sys.argv[5]
        print 'Filtering to only process lanes %s'%(selected_lanes)
    if len(sys.argv) > 6:
        selected_parts = sys.argv[6]
        print 'Filtering to only process parts %s'%(selected_parts)

    lane_numbers = map(int, selected_lanes.split(","))
    part_numbers = map(int, selected_parts.split(","))

    try:
        os.mkdir(output_dir)
    except Exception as e:
        print "Unable to create directory: %s"%(e,)

    items = []

    out_paths = []

    input_dict = json.load(open(input_path, "r"))
    parts = input_dict["parts"]
    for name in parts:
        part = parts[name]
        basename = os.path.basename(part["reads"])
        #if basename == "reads_GS21184-FS3-L05_012.tsv.bz2":
        if basename:
            # cgatools map2sam --reads=$reads --mappings=$mappings  --reference=$reference > $output

            # parse lane number and part number out of:
            # the reads file name reads_GS21910-FS3-L04_004.tsv.bz2
            basename = os.path.basename(part["reads"])
            sans_extension = basename.split(".")[0]
            lane_part_segment = sans_extension.split("-")[-1]
            lane_string, part_string = lane_part_segment.split("_")
            lane_number = int(lane_string.replace("L", ""))
            part_number = int(part_string)

            if part_number in part_numbers and lane_number in lane_numbers:
                #  $ cgapipe.sh output.bam output.fastq --reads=reads.tsv \
                #               --mappings=mappings.tsv --reference=ref.crr
                bam_path = os.path.join(output_dir, "%s.bam" % name)
                out_paths.append(bam_path)
                fastq_path = os.path.join(output_dir, "%s.fastq" % name)
                out_paths.append(fastq_path)
                args = [bam_path, fastq_path,
                        "--reads=%s" % part["reads"],
                        "--mappings=%s" % part["mapping"],
                        "--reference=%s" % reference_file]
                out_path = os.path.join(output_dir, "%s.out"%(name,))
                out_paths.append(out_path)
                items.append(CondorQueueItem(args, out_path))

    print "Spawning job with %s items"%(len(items),)
    job = CondorJob("/nfs/software/bin/cgapipe.sh", items=items)
    print job.condor_job

    cluster_id = job.submit()
    print "Job submitted with cluster_id %s" % cluster_id
    result = job.wait()

    if job.logdata:
        sys.stdout.write(job.logdata)
    if job.outdata:
        sys.stdout.write(job.outdata)
    if job.errdata:
        # galaxy looks for output to stderr to determine "failure",
        # not the exit code, so we send stderr to stdout unless we got
        # a failure exit code:
        if result == 0:
            sys.stdout.write(job.errdata)
        else:
            sys.stderr.write(job.errdata)

    json.dump({"title":"multiple sam files",
               "files":out_paths,}, open(output_json_path, "w"))

