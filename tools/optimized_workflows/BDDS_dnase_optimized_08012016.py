# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the Trena workflow
See below for options
"""

import requests, glob, time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))


def file_type(filename):
    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
        }

    max_len = max(len(x) for x in magic_dict)

    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "txt"

def get_ncores():
    return multiprocessing.cpu_count()

def get_linuxRAM():
    # return the available memory on the instance in MB
    totalMemory = os.popen("free -m").readlines()[1].split()[1]
    return int(totalMemory)

def get_readLength(input_fastq):
    if file_type(input_fastq) == "gz":
        fh_zip = gzip.open(input_fastq)
        fh_zip.readline()
        return len(fh_zip.readline().rstrip()) 
    else:
        return len(os.popen("head -n 2 %s" % input_fastq).readlines()[1].rstrip())

def get_sampleName(input_fastq):
    if input_fastq.endswith("fastq.gz"):
        name = os.path.basename(input_fastq)
        return (".").join(name.split(".")[:-2])
    else:
        name = os.path.basename(input_fastq)
        return (".").join(name.split(".")[:-1])

def run_cmd_parallel( cmd, wd_tmpdir, descriptor):
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    p = subprocess.Popen(args=cmd, stderr=stderr, shell=True)
    return p

def run_cmd ( cmd , wd_tmpdir, descriptor):
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    proc = subprocess.Popen( args=cmd, stderr=stderr, shell=True )

    exit_code = proc.wait()

    if exit_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
    stderr.flush()
    stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    stderr.close()

def _get_input_urls(sample_id, config_file, input_dir):
    urls = []
    fh = open(config_file, "r")
    for line in fh:
        line = line.rstrip("\n")
        name, url = line.split("\t")
        if name == sample_id:
            urls.append(url)
            filename = "%s/%s" % (input_dir, os.path.basename(url))
            r = requests.get(url)
            with open(filename, 'wb') as test:
                test.write(r.content)
    return urls

def create_wf_meta(options, output_dir, tmp_dir, input_dir):

    # get ncores
    ncpu = get_ncores()/4
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()/4

    # get the reference datasets
    reference_db = options.snap_ref
    fasta_path  = options.fasta_ref
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"

    ## Setup workflow
    command_meta = {}
    sample_name = "dnase_seq"
    #file_urls = _get_input_urls(options.sample_id, options.config_file)

    #
    step = 0
    command_meta[step] = []
    outputF1 = "%s/%s-R1.fastq.gz" % (output_dir, sample_name)
    input_files = []
    output_files = [outputF1]
    #for url in file_urls:
    cmd = "python /opt/galaxy/tools/optimized_workflows/encode-bag-client.py %s %s/" % (options.sample_id, output_dir)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #
    step = 1
    command_meta[step] = []
    input_files = sorted(glob.glob("%s/*/*/data/*.fastq.gz" % output_dir))
    outputF1 = "%s/%s-R1.fastq.gz" % (output_dir, sample_name)
    output_files = [outputF1]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s/*/*/data/*.fastq.gz" % (outputF1, output_dir)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/snap/snap-alignment.py %s %s %s" % (reference_db, outputF, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sambamba/sambamba_sort.py --input=%s --order=coordinate --output=%s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "bamToBed -i %s > %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.fseq.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "mkdir tmp-fseq; cd tmp-fseq; fseq -f 0 -l 600 -of bed -s 1 -t 2.5 %s; cat *.bed >  %s; cd ..; rm -rf tmp-fseq" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.fseq.cut1.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c2,c3\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.fseq.cut1.compute1.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/stats/column_maker.py %s %s \"c2-c3\" no 3 \"str,int,int\"" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.fseq.cut1.compute1.filter1.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/stats/filtering.py %s %s \"c4__lt__-400\" 4 \"str,int,int,float\" 0" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 9
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.fseq.cut1.compute1.filter1.cut2.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c2,c3\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 10
    command_meta[step] = []
    inputF = command_meta[step-7][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 11
    command_meta[step] = []
    inputBed = command_meta[step-2][0]["output_file"]
    inputBam = command_meta[step-8][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.wellington.bed" % (output_dir, sample_name)
    input_files = [inputBed, inputBam, inputIndex]
    output_files = [outputF]
    tmpDir = "./wellingtonTmp"
    cmd = "mkdir %s; wellington_footprints.py %s %s %s; cp %s/p\ value\ cutoffs/%s.%s.WellingtonFootprints.-10.bed %s; rm -rf %s" % (tmpDir, inputBed, inputBam, tmpDir, tmpDir, os.path.basename(inputBam), os.path.basename(inputBed), outputF, tmpDir)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 12
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.wellington.bed.compute1.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/stats/column_maker.py %s %s \"c2-5\" yes 6 \"str,int,int,str,float,str\"" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 13
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.wellington.bed.compute2.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/stats/column_maker.py %s %s \"c3+5\" yes 7 \"str,int,int,str,float,str,int\"" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 14
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.wellington.bed.compute2.cut1.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c7,c8,c4,c5,c6\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 15
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.wellington.bed.compute2.cut1.converted.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/convert_characters.py %s U %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 16
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-snap-out.sorted.wellington.bed.compute2.cut1.converted.cut2.bed" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c5,c6,c7,c8,c9\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 17
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-extracted.fasta" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/extract/extract_genomic_dna.py %s %s -o fasta -d hg38 -1 1,2,3,0 -F %s" % (inputF, outputF, fasta_path)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 18
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-meme.txt" % (output_dir, sample_name)
    outputHtml = "%s/%s-meme.html" % (output_dir, sample_name)
    outputInterval = "%s/%s-meme.interval" % (output_dir, sample_name)
    outputGff = "%s/%s-meme.gff" % (output_dir, sample_name)
    outputXml = "%s/%s-meme.xml" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputHtml, outputInterval, outputGff, outputXml]
    tmpDirMeme = "./meme-tmp"
    cmd = "python /opt/galaxy/tools/meme/fimo_wrapper.py \'fimo --o %s --verbosity 1 %s %s\' %s %s %s %s %s %s; rm -rf %s" % (tmpDirMeme, options.meme_db, inputF, tmpDirMeme, outputHtml, outputInterval, outputF, outputXml, outputGff, tmpDirMeme)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 19
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-fimo.cut1.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c2\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 20
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_file"]
    outputF = "%s/%s-fimo.cut2.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c3,c4,c5,c6,c7,c8,c9\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 21
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_file"]
    outputF = "%s/%s-fimo.cut1.converted.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/convert_characters.py %s U %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 22
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-fimo.cut.converted.cut2.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c2,c3,c4\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 23
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    inputF2 = command_meta[step-3][0]["output_file"]
    outputF = "%s/%s-fimo.cut.converted.cut2.pasted.txt" % (output_dir, sample_name)
    input_files = [inputF, inputF2]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/pasteWrapper.pl %s %s T %s" % (inputF, inputF2, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 24
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-all-footprints.bed" % (output_dir, sample_name)
    outputFRest = options.output_bed
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/sorter.py --input=%s --out_file1=%s --column=1 --style=alpha --order=ASC 2 num ASC; cp -r %s %s" % (inputF, outputF, outputF, outputFRest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 25
    command_meta[step] = []
    inputF = command_meta[step-1][0]['output_file']
    outputF = options.output_log
    input_files = [inputF]
    output_files = [outputF]
    cmd = "minid --test --register --title 'DNASE Fingerprint for accession ID %s' %s 2> %s " % (options.sample_id, outputFRest, options.output_minid)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    return command_meta

def job_step_dependency(meta, current_step_id):
    step_meta = meta[current_step_id]
    input_files = []
    for job in step_meta:
        for ifile in job['input_files']:
            if ifile not in input_files:
                input_files.append(ifile)

    #print "INPUT_FILES: %s" % input_files

    keep_file = []
    for input_file in input_files:
        flag = None
        for step_id, step_meta in meta.iteritems():
            if step_id > current_step_id:
                for job in step_meta:
                    if input_file in job['input_files']:
                        keep_file.append(input_file)
                        flag = 1
                        break
                if flag:
                    #print "Can't delete file uses in step %s: %s" % (step_id, input_file)
                    break
    return keep_file

def delete_file(input_files, keep_files, step_id, final_outputs):
    for ifile in input_files:
        if ifile not in keep_files and ifile not in final_outputs:
            print "DELETING INPUT: %s" % ifile
            #os.remove(ifile)

def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--snap-ref', dest="snap_ref", help='The BWA reference genome to use or index' )
    parser.add_option( '--input-fasta', dest="fasta_ref", help='The GATK reference genome to use or index' )
    parser.add_option( '--meme-db', dest="meme_db", help='The (forward) fastq file to use for the mapping' )
    parser.add_option( '--config', dest="config_file", help='The reverse fastq file to use for mapping if paired-end data' )
    parser.add_option( '', '--sample', dest="sample_id", help='The bam input' )
    parser.add_option( '', '--output-bed', dest="output_bed", help="The file to save the output (BAM format)" )
    parser.add_option( '-l', '--output-log', dest="output_log", help='The log output (Txt format)' )
    parser.add_option( '', '--output-minid', dest="output_minid", help="The minid output" )
    (options, args) = parser.parse_args()

    final_outputs = [options.output_bed, options.output_log, options.output_minid]
    if len(args) > 0:
        parser.error('Wrong number of arguments')

#    # make temp directory 
    #output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-dnaExome")
    #tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimized-tmp-")
    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")
    input_dir = "%s/inputs" % output_dir
    os.makedirs(input_dir)

    #file_urls = _get_input_urls(options.sample_id, options.config_file, input_dir)

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir, input_dir)

    #print wf_meta
    for step_id, step_list in wf_meta.iteritems():
        print "START STEP %s: %s" % (step_id, time_stamp())
        step_input_files = []
        if len(step_list) == 1:
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                run_cmd ( job['cl'], tmp_dir, "running some job")
                for ifile in job['input_files']:
                    if ifile not in step_input_files:
                        step_input_files.append(ifile)
        else:   # run jobs in parallel
            ps = {}
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                p = run_cmd_parallel ( job['cl'], tmp_dir, "running some job")
                ps[p.pid] = p
            print "Waiting for %d processes..." % len(ps)
            while ps:
                pid, status = os.wait()
                if pid in ps:
                    del ps[pid]
                    print "Waiting for %d processes..." % len(ps)
        print "END STEP %s: %s" % (step_id, time_stamp())


        # job completed, check to see if there are inputs that should be deleted
        if step_id > 0:
            keep_files = job_step_dependency(wf_meta, step_id)
            delete_file(step_input_files, keep_files, step_id, final_outputs)


    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
