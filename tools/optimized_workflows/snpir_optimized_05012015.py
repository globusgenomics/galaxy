# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
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

def run_cmd_parallel( cmd, wd_tmpdir, descriptor):
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    p = subprocess.Popen(args=cmd, stderr=stderr, shell=True)
    return p

def get_sampleName(input_fastq):
    if input_fastq.endswith("fastq.gz"):
        name = os.path.basename(input_fastq)
        return (".").join(name.split(".")[:-2])
    else:
        name = os.path.basename(input_fastq)
        return (".").join(name.split(".")[:-1])

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

def create_wf_meta(options, output_dir, tmp_dir):

    # reads
    fastq = options.fastq
    if options.rfastq:
        rfastq = options.rfastq
    read_length = get_readLength(fastq)

    # get sample name
    sample_name = get_sampleName(fastq)

    # get ncores
    ncpu = get_ncores()
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()

    # get the reference datasets
    hg19_reference="/mnt/galaxyIndices/genomes/Hsapiens/hg19/seq/hg19_jxn.fa"
    RepeatMasker="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/repeatmasker.bed" 
    gene_annotation="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/hg19_gene_annotation.gtf" 
    rnaedit="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/rnaedit.bed"
    dbsnp="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/dbsnp_138.hg19.vcf"
    mills="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf"
    g1000="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/1000G_phase1.indels.hg19.vcf"
    GATK3_PATH = "/mnt/galaxyTools/tools/gatk3/GenomeAnalysisTK-3.4-46"
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"

    #bwa_index = None
    #if read_length <= 45:
    #    bwa_index="/scratch/averahealth/galaxy/data/genomes/Hsapiens/hg19/bwa_0.7.10_junctions/45bp/hg19_genome_junctions_45.fa"
    #elif read_length <= 75:
    #    bwa_index="/scratch/averahealth/galaxy/data/genomes/Hsapiens/hg19/bwa_0.7.10_junctions/75bp/hg19_genome_junctions_75.fa"
    #elif read_length <= 95:
    #    bwa_index="/scratch/averahealth/galaxy/data/genomes/Hsapiens/hg19/bwa_0.7.10_junctions/95bp/hg19_genome_junctions_95.fa"
    #elif read_length <= 140:
    #    bwa_index="/scratch/averahealth/galaxy/data/genomes/Hsapiens/hg19/bwa_0.7.10_junctions/140bp/hg19_genome_junctions_140.fa"

    bwa_index=options.ref
    ## Setup workflow
    command_meta = {}

    read_group = "--rgid=\"%s\" --rglb=\"%s\" --rgpl=\"%s\" --rgsm=\"%s\"" % (options.rgid, options.rglb, options.rgpl, options.rgsm)
    if options.rgcn:
        read_group += " --rgcn=\"%s\"" % options.rgcn
    if options.rgds:
        read_group += " --rgds=\"%s\"" % options.rgds
    if options.rgdt:
        read_group += " --rgdt=\"%s\"" % options.rgdt
    if options.rgfo:
        read_group += " --rgfo=\"%s\"" % options.rgfo
    if options.rgks:
        read_group += " --rgks=\"%s\"" % options.rgks
    if options.rgpg:
        read_group += " --rgpg=\"%s\"" % options.rgpg
    if options.rgpi:
        read_group += " --rgpi=\"%s\"" % options.rgpi
    if options.rgpu:
        read_group += " --rgpu=\"%s\"" % options.rgpu


    #
    step = 0
    command_meta[step] = []
    inputF = fastq
    input_files = [inputF]
    outputF = "%s/bwa-out1.sam" % output_dir
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sr_mapping/bwa_wrapper.py --threads=%s --fileSource=indexed --ref=%s --do_not_build_index --input1=%s --output=%s --genAlignType=single --params=full --maxEditDist=0 --fracMissingAligns=0.04 --maxGapOpens=1 --maxGapExtens=-1 --disallowLongDel=16 --disallowIndel=5 --seed=-1 --maxEditDistSeed=2 --mismatchPenalty=3 --gapOpenPenalty=11 --gapExtensPenalty=4 --suboptAlign=\"\" --noIterSearch=false --outputTopN=3 --outputTopNDisc=10 --maxInsertSize=500 --maxOccurPairing=100000 %s --suppressHeader=false" % (nthreads, bwa_index, inputF, outputF, read_group)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    inputF = rfastq
    input_files = [inputF]
    outputF = "%s/bwa-out2.sam" % output_dir
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sr_mapping/bwa_wrapper.py --threads=%s --fileSource=indexed --ref=%s --do_not_build_index --input1=%s --output=%s --genAlignType=single --params=full --maxEditDist=0 --fracMissingAligns=0.04 --maxGapOpens=1 --maxGapExtens=-1 --disallowLongDel=16 --disallowIndel=5 --seed=-1 --maxEditDistSeed=2 --mismatchPenalty=3 --gapOpenPenalty=11 --gapExtensPenalty=4 --suboptAlign=\"\" --noIterSearch=false --outputTopN=3 --outputTopNDisc=10 --maxInsertSize=500 --maxOccurPairing=100000 %s --suppressHeader=false" % (nthreads, bwa_index, inputF, outputF, read_group)
    command_meta[step].append({"cl":cmd, "input_files": input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = command_meta[step-1][1]["output_file"]
    outputF = "%s/grep1.sam" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/grep.py -i %s -o %s -pattern '^@' -v true" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 2
    command_meta[step] = []
    outputF = "%s/merged.sam" % output_dir
    input_files = [command_meta[0][0]["output_file"], command_meta[1][0]["output_file"]]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s %s" % (outputF, command_meta[0][0]["output_file"], command_meta[1][0]["output_file"])
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/merged.conv.sam" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "cat %s | java -Xmx2g convertCoordinates  > %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/header.sam" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/grep.py -i %s -o %s -pattern '^@' -v false" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[4][0]["output_file"]
    outputF = "%s/t1.sam" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/headWrapper.pl %s 25 %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    inputF = command_meta[4][0]["output_file"]
    outputF = "%s/t2.sam" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/tailWrapper.pl %s 2 %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    input1 = command_meta[5][0]["output_file"]
    input2 = command_meta[5][1]["output_file"]
    outputF = "%s/new_header.sam" % output_dir
    input_files = [input1, input2]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s %s" % (outputF, input1, input2)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    input1 = command_meta[3][0]["output_file"]
    input2 = command_meta[6][0]["output_file"]
    outputF = "%s/merged.conv.nh.sam" % output_dir
    input_files = [input1, input2]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --input %s -o %s --header-file %s --output-format sam -j \"%s/picard.jar ReplaceSamHeader\" --tmpdir %s" % (input1, outputF, input2, JAVA_JAR_PATH, tmp_dir)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputF = command_meta[7][0]["output_file"]
    outputF = "%s/merged.conv.nh.sort.bam" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/samtools/samtools_filter_bam.py -p '-bS' -p '-q 20' -p '-F 4' --input %s --output %s --sorted-bam" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 9
    command_meta[step] = []
    inputF = command_meta[8][0]["output_file"]
    outputF = "%s/merged.conv.nh.sort.rd.bam" % output_dir
    outputMetrics = "%s/SM1.dups" % output_dir
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --remdups true --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --optdupdist 100 -j \"%s/picard.jar MarkDuplicates\" -d %s -t %s -e bam" % (jvm_heap, inputF, tmp_dir, outputF, JAVA_JAR_PATH, tmp_dir, outputMetrics)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 10
    command_meta[step] = []
    inputF = command_meta[9][0]["output_file"]
    outputF = "%s/merged.conv.nh.sort.rd.bam.bai" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 11
    command_meta[step] = []
    inputF = command_meta[9][0]["output_file"]
    inputIndex = command_meta[10][0]["output_file"]
    outputF = "%s/output.intervals" % output_dir
    outputLog = "%s/output.intervals.log" % output_dir
    input_files = [inputF, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator -o %s --num_threads %s -R %s\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_1\" -p \"--filter_reads_with_N_cigar\"" % (outputLog, inputF, inputIndex, GATK3_PATH, outputF, ncpu, hg19_reference, mills, g1000)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 12
    command_meta[step] = []
    inputBam = command_meta[9][0]["output_file"]
    inputF = command_meta[11][0]["output_file"]
    inputIndex = command_meta[10][0]["output_file"]
    outputF = "%s/merged.conv.sort.rd.realigned.bam" % output_dir
    outputLog = "%s/merged.conv.sort.rd.realigned.log" % output_dir
    input_files = [inputF, inputIndex, inputBam]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T IndelRealigner -o %s -R %s -LOD 0.4 --filter_reads_with_N_cigar --consensusDeterminationModel KNOWNS_ONLY\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_1\"  -d \"-targetIntervals\" \"%s\" \"gatk_interval\" \"gatk_target_intervals\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, outputF, hg19_reference, mills, g1000, inputF )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 13
    command_meta[step] = []
    inputBam = command_meta[12][0]["output_file"]
    inputIndex = "%s/merged.conv.sort.rd.realigned.bam.bai" % output_dir
    outputF = "%s/recal_data.table" % output_dir
    outputLog = "%s/merged.conv.sort.rd.realigned.log" % output_dir
    input_files = [inputBam]
    output_files = [outputF, outputLog, inputIndex]
    cmd = "samtools index %s && python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\"  \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator -nct %s -R %s --out %s\" -d \"--knownSites:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_1\"" % (inputBam, outputLog, inputBam, inputIndex, GATK3_PATH, ncpu, hg19_reference, outputF, dbsnp, g1000, mills )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 14
    command_meta[step] = []
    inputBam = command_meta[12][0]["output_file"]
    inputIndex = "%s/merged.conv.sort.rd.realigned.bai" % output_dir
    inputTable = command_meta[13][0]["output_file"]
    outputF = "%s/post_recal_data.table" % output_dir
    outputLog = "%s/merged.conv.sort.rd.realigned.log" % output_dir
    input_files = [inputBam, inputTable, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator -nct %s --BQSR %s -R %s --out %s\" -d \"--knownSites:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_1\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, ncpu, inputTable, hg19_reference, outputF, dbsnp, g1000, mills )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 15
    command_meta[step] = []
    inputTable1 = command_meta[13][0]["output_file"]
    inputTable2 = command_meta[14][0]["output_file"]
    outputF = "%s/recalibration_plots.pdf" % output_dir
    outputLog = "%s/recalibration_plots.log" % output_dir
    input_files = [inputTable1, inputTable2]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 -p \"java -jar %s/GenomeAnalysisTK.jar -T AnalyzeCovariates -before %s -after %s -R %s -csv %s -plots %s\" " % (GATK3_PATH, inputTable1, inputTable2, hg19_reference, outputLog, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 16
    command_meta[step] = []
    inputBam = command_meta[12][0]["output_file"]
    inputTable = command_meta[13][0]["output_file"]
    inputIndex = "%s/merged.conv.sort.rd.realigned.bai" % output_dir
    #outputF = "%s/merged.conv.sort.rd.realigned.recal.bam" % output_dir
    outputF = options.output_bam
    outputLog = "%s/merged.conv.sort.rd.realigned.recal.log" % output_dir
    input_files = [inputTable, inputBam, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T PrintReads -o %s -nct %s -R %s --BQSR %s\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, outputF, ncpu, hg19_reference, inputTable)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 17
    command_meta[step] = []
    inputBam = command_meta[16][0]["output_file"]
    inputIndex = "%s/merged.conv.sort.rd.realigned.recal.bai" % output_dir
    outputF = "%s/raw_variants.vcf" % output_dir
    outputLog = "%s/raw_variants.log" % output_dir
    outputMetrics = "%s/raw_variants_metrics.txt" % output_dir
    input_files = [inputBam]
    output_files = [outputF, outputLog, outputMetrics, inputIndex]
    cmd = "samtools index %s && python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input_0\"  -d \"\" \"%s\" \"bam_index\" \"gatk_input_0\" -d \"--dbsnp:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T UnifiedGenotyper --num_threads %s --out %s --metrics_file %s -R %s --read_filter BadCigar --output_mode EMIT_VARIANTS_ONLY -stand_call_conf 0 -stand_emit_conf 0\"" % (inputBam, jvm_heap, outputLog, inputBam, inputIndex, dbsnp, GATK3_PATH, ncpu, outputF, outputMetrics, hg19_reference)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 18
    command_meta[step] = []
    inputF = command_meta[17][0]["output_file"]
    outputF = "%s/raw_variants.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "convertVCF.sh %s %s 20" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 19
    command_meta[step] = []
    inputF = command_meta[18][0]["output_file"]
    inputBam = command_meta[16][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.txt" % output_dir
    input_files = [inputF, inputBam]
    output_files = [outputF]
    cmd = "filter_mismatch_first6bp.pl -infile %s -outfile %s -bamfile %s" % (inputF, outputF, inputBam)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 20
    command_meta[step] = []
    inputF = command_meta[19][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.add_column.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/stats/column_maker.py %s %s \"c2-1\" yes 6 \"str,int,list,str,str,float\"" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 21
    command_meta[step] = []
    inputF = command_meta[20][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.cut_column.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c7,c2,c3,c4,c5,c6\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 22
    command_meta[step] = []
    inputF = command_meta[21][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.intersect.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "intersectBed -a %s  -b %s -v > %s" % (inputF, RepeatMasker, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 23
    command_meta[step] = []
    inputF = command_meta[22][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c3,c4,c5,c6,c7\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 24
    command_meta[step] = []
    inputF = command_meta[23][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.rmintron.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "filter_intron_near_splicejuncts.pl -infile %s -outfile %s -genefile %s -splicedist 4" % (inputF, outputF, gene_annotation)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 25
    command_meta[step] = []
    inputF = command_meta[24][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.rmintron.rmhom.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "filter_homopolymer_nucleotides.pl -infile %s -outfile %s -refgenome %s" % (inputF, outputF, hg19_reference)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 26
    command_meta[step] = []
    inputF = command_meta[25][0]["output_file"]
    inputBam = command_meta[16][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt" % output_dir
    input_files = [inputF, inputBam]
    output_files = [outputF]
    cmd = "samtools index -b %s && BLAT_candidates.pl -infile %s -outfile %s -bamfile %s -refgenome %s" % (inputBam, inputF, outputF, inputBam, hg19_reference)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 27
    command_meta[step] = []
    inputF = command_meta[26][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.add_column.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/stats/column_maker.py %s %s \"c2-1\" yes 6 \"str,int,list,str,str,float\"" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 28
    command_meta[step] = []
    inputF = command_meta[27][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.cut_column.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/filters/cutWrapper.pl %s \"c1,c7,c2,c3,c4,c5,c6\" T %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 29
    command_meta[step] = []
    inputF = command_meta[28][0]["output_file"]
    outputF = "%s/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "intersectBed -a %s -b %s -v > %s" % (inputF, rnaedit, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 30
    command_meta[step] = []
    inputF = command_meta[17][0]["output_file"]
    outputF = "%s/raw_variants.header.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/grep.py -i %s -o %s -pattern '^#' -v false" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    outputF = "%s/raw_variants.non_header.txt" % output_dir
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/grep.py -i %s -o %s -pattern '^#' -v true" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 31
    command_meta[step] = []
    inputF = command_meta[29][0]["output_file"]
    inputVCF = command_meta[30][1]["output_file"]
    outputF = "%s/red_variants.vcf" % output_dir
    input_files = [inputF, inputVCF]
    output_files = [outputF]
    cmd = "Rscript /opt/galaxy/tools/snpir/extractvcf.R %s %s %s" % (inputF, inputVCF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 32
    command_meta[step] = []
    inputF = command_meta[31][0]["output_file"]
    inputHeader = command_meta[30][0]["output_file"]
    #outputF = "%s/final_variants.vcf" % output_dir
    outputF = options.output_vcf
    input_files = [inputF, inputHeader]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s %s" % (outputF, inputHeader, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 33
    command_meta[step] = []
    inputF = command_meta[32][0]["output_file"]
    input_files = [inputF]
    file_name = "%s/%s.vcf" % (options.destination_dir, sample_name)
    cmd = "python /opt/galaxy/tools/globus/transfer.py %s -k -c -o transfer.log --goauth-token \"%s\" --source-ep=\"%s\" --source-path=\"%s\"  --destination-ep=\"%s\" --destination-path=\"%s\" --deadline=%s" % (options.username, options.goauth_token, options.source_ep, inputF, options.destination_ep, file_name, options.deadline)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 34
    command_meta[step] = []
    inputF = command_meta[16][0]["output_file"]
    input_files = [inputF]
    file_name = "%s/%s.bam" % (options.destination_dir, sample_name)
    cmd = "python /opt/galaxy/tools/globus/transfer.py %s -k -c -o transfer.log --goauth-token \"%s\" --source-ep=\"%s\" --source-path=\"%s\"  --destination-ep=\"%s\" --destination-path=\"%s\" --deadline=%s" % (options.username, options.goauth_token, options.source_ep, inputF, options.destination_ep, file_name, options.deadline)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 35
    command_meta[step] = []
    inputF = options.output_log
    input_files = [inputF]
    file_name = "%s/%s.log" % (options.destination_dir, sample_name)
    cmd = "python /opt/galaxy/tools/globus/transfer.py %s -k -c -o transfer.log --goauth-token \"%s\" --source-ep=\"%s\" --source-path=\"%s\"  --destination-ep=\"%s\" --destination-path=\"%s\" --deadline=%s" % (options.username, options.goauth_token, options.source_ep, inputF, options.destination_ep, file_name, options.deadline)
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
    parser.add_option( '--ref', dest="ref", help='The reference genome to use or index' )
    parser.add_option( '-f', '--fastq', help='The (forward) fastq file to use for the mapping' )
    parser.add_option( '-F', '--rfastq', help='The reverse fastq file to use for mapping if paired-end data' )
    parser.add_option( '', '--output-vcf', dest="output_vcf", help="The file to save the output (BAM format)" )
    parser.add_option( '-u', '--output-bam', dest="output_bam", help='The file to save the output (VCF format)' )
    parser.add_option( '-l', '--output-log', dest="output_log", help='The log output (Txt format)' )
    parser.add_option( '--rgid', help='Read group identifier' )
    parser.add_option( '--rgsm', help='Sample' )
    parser.add_option( '--rgpl', choices=[ 'CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'PACBIO' ], help='Platform/technology used to produce the reads' )
    parser.add_option( '--rglb', help='Library name' )
    parser.add_option( '--rgpu', help='Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD)' )
    parser.add_option( '--rgcn', help='Sequencing center that produced the read' )
    parser.add_option( '--rgds', help='Description' )
    parser.add_option( '--rgdt', help='Date that run was produced (ISO8601 format date or date/time, like YYYY-MM-DD)' )
    parser.add_option( '--rgfo', help='Flow order' )
    parser.add_option( '--rgks', help='The array of nucleotide bases that correspond to the key sequence of each read' )
    parser.add_option( '--rgpg', help='Programs used for processing the read group' )
    parser.add_option( '--rgpi', help='Predicted median insert size' )
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--goauth-token', dest="goauth_token", help='Globus token' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')

    final_outputs = [options.output_bam, options.output_log, options.output_vcf]

#    # make temp directory 
    #output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-snpir")
    #tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimized-tmp-")
    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir)

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
