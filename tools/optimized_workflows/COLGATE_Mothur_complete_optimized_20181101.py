# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the Mother pipeline in optimized mode in one Node
See below for options
"""

import workflows, time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
import ast
CHUNK_SIZE = 2**20 #1mb

def create_wf_meta(options, output_dir, tmp_dir):
    # reads
    sample_name = options.sample_name

    if options.transfer_info:
        transfer_info_in = ast.literal_eval(options.transfer_info)
        in_transfer_info = transfer_info_in
        in_transfer_info['to_path'] = "%s/%s-samples_dir" % (output_dir, sample_name)

        transfer_info_out = ast.literal_eval(options.transfer_info_out)
        out_transfer_info = transfer_info_out
        #out_transfer_info['from_path'] = "%s/%s-samples_dir" % (output_dir, sample_name)

        cred_file_dir = '/home/galaxy/.globusgenomics/tmp'
        cred_file_loc = os.path.join(cred_file_dir, in_transfer_info['cred_file'])
        local_cred_file = os.path.join(output_dir, os.path.basename(cred_file_loc))
        shutil.copyfile(cred_file_loc, local_cred_file)

        cred_file_loc_out = os.path.join(cred_file_dir, transfer_info_out['cred_file'])
        local_cred_file_out = os.path.join(output_dir, os.path.basename(cred_file_loc_out))
        shutil.copyfile(cred_file_loc_out, local_cred_file_out)

    lcm = 1
    #lcm = 4

    # get ncores
    ncpu = workflows.get_ncores()/lcm
    nthreads = ncpu/2

    # get memory
    jvm_heap = workflows.get_linuxRAM()/lcm

    # get the reference datasets
    #picard_reference_db = options.picard_ref

    ## Setup workflow
    command_meta = {}
    ##outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    outputLog = options.output_log

    #
    step = 0
    command_meta[step] = []
    #input_files = [options.input_dir]
    input_files = []
    outputF1 = "%s/%s-samples_dir" % (output_dir, sample_name)
    outputF2 = "%s/hellow.world.txt" % (output_dir)
    output_files = [outputF1, outputF2]
    #cmd = "echo S3 transfer" % (outputF1, " ".join(input_files), outputF1, outputFrest)
    if options.input_dir:
        cmd = "mkdir %s; cp -r %s/* %s/." % (outputF1, options.input_dir, outputF1)
    else:
        cmd = "echo hello > %s; cp %s %s; mkdir %s; python /opt/galaxy/tools/globus/s3_transfer.py --transfer-info \"%s\"" % (outputF2, local_cred_file, cred_file_loc, outputF1, in_transfer_info)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #
    step = 1
    ncpu_step = int(ncpu)
    command_meta[step] = []
    input1 = command_meta[step-1][0]["output_file"]
    input_files = [input1]
    outputF = "%s/combo_fastq.dat" % (output_dir)
    outputF2 = "%s/sample_mapping.dat" % (output_dir)
    output_files = [outputF, outputF2]
    cmd = "gunzip %s/*.gz; gzip %s/*; echo \"#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\" > %s ;for i in %s/*R1*; do samplename=`basename $i|cut -f1 -d '_'`;R1=`ls %s/$samplename'_'*R1*`; R2=`ls %s/$samplename'_'*R2*`; echo \"$samplename\t$R1\t$R2\" >> %s; echo \"$samplename\t\t\t$samplename\" >> %s;done" % (input1, input1, outputF2, input1, input1, input1, outputF, outputF2)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    ## s3 transfer test
    ##step = 1
    ##command_meta[step] = []
    #inputF = command_meta[step-1][0]["output_files"][1] # taxonomy - otu
    #input_files = [inputF]
    #outputF = "%s/s3_transfer.log" % (output_dir)
    #output_files = [outputF]
    #transfer_info_curr = transfer_info_out
    #transfer_info_curr['rename'] = 'HELLO_WORLD.txt'
    ##transfer_info_curr['sse'] = True
    #transfer_info_curr['from_path'] = inputF
    #transfer_info_curr['object_type'] = 'file'
    #cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
    #command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # make.contigs
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][0]
    outputF_trim_contigs_fasta = "%s/combo_fastq.trim.contigs.fasta" % (output_dir)
    outputF_trim_contigs_qual = "%s/combo_fastq.trim.contigs.qual" % (output_dir)
    outputF_scrap_contigs_fasta = "%s/combo_fastq.scrap.contigs.fasta" % (output_dir)
    outputF_scrap_contigs_qual = "%s/combo_fastq.scrap.contigs.qual" % (output_dir)
    outputF_contigs_report = "%s/combo_fastq.contigs.report" % (output_dir)
    outputF_contigs_groups = "%s/combo_fastq.contigs.groups" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputLog, outputF_trim_contigs_fasta, outputF_trim_contigs_qual, outputF_scrap_contigs_fasta, outputF_scrap_contigs_qual, outputF_contigs_report, outputF_contigs_groups]
    cmd = "echo 'make.contigs( file=%s, align=needleman, match=1, mismatch=-1, gapopen=-2, gapextend=-1, insert=20, processors=%s )' | sed 's/ //g' | mothur | tee -a %s;" % (inputF, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.seqs
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1]
    input_files = [inputF]
    outputF_summary = "%s/combo_fastq.trim.contigs.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary]
    cmd = "echo 'summary.seqs( fasta=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Screen.seqs
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1]
    inputF2 = command_meta[step-2][0]["output_files"][6]
    input_files = [inputF, inputF2]
    outputF_good_fasta = "%s/combo_fastq.trim.contigs.good.fasta" % (output_dir)
    outputF_bad_acc_nos = "%s/combo_fastq.trim.contigs.bad.accnos" % (output_dir)
    outputF_good_groups = "%s/combo_fastq.contigs.good.groups" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_good_fasta, outputF_bad_acc_nos, outputF_good_groups]
    cmd = "echo 'screen.seqs( fasta=%s, maxlength=475 ,maxambig=0, group=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Unique.seqs
    step = 5
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1]
    input_files = [inputF]
    outputF_unique_fasta = "%s/combo_fastq.trim.contigs.good.unique.fasta" % (output_dir)
    outputF_unique_names = "%s/combo_fastq.trim.contigs.good.names" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_unique_fasta, outputF_unique_names]
    cmd = "echo 'unique.seqs( fasta=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.seqs
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1]
    input_files = [inputF]
    outputF_summary = "%s/combo_fastq.trim.contigs.good.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary]
    cmd = "echo 'summary.seqs( fasta=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # count.seqs
    step = 7
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][2]
    inputF2 = command_meta[step-3][0]["output_files"][3]
    input_files = [inputF]
    outputF_counts = "%s/combo_fastq.trim.contigs.good.count_table" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_counts]
    cmd = "echo 'count.seqs( name=%s, group=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})


    # align.seqs
    step = 8
    command_meta[step] = []
    inputF = command_meta[step-3][0]["output_files"][1]
    input_files = [inputF]
    outputF_seqs = "%s/combo_fastq.trim.contigs.good.unique.align" % (output_dir)
    outputF_reports = "%s/combo_fastq.trim.contigs.good.unique.align.report" % (output_dir)
    outputF_flips = "%s/combo_fastq.trim.contigs.good.unique.flip.accnos" % (output_dir)
    outputLogAlign = "%s/%s.mother.out.align.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_seqs, outputF_reports, outputF_flips, outputLogAlign]
    cmd = "echo 'align.seqs( fasta=%s, reference=%s, align=needleman, ksize=8, processors=%s )' | sed 's/ //g' | mothur | tee -a %s; cat %s >> %s" % (inputF, options.align_reference, ncpu, outputLogAlign, outputLogAlign, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.seqs
    step = 9
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1]
    inputF2 = command_meta[step-2][0]["output_files"][1]
    input_files = [inputF]
    outputF_summary = "%s/combo_fastq.trim.contigs.good.unique.summary" % (output_dir)
    outputLogSummary = "%s/%s.mother.out.summary.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary, outputLogSummary]
    cmd = "echo 'summary.seqs( fasta=%s, count=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s; cat %s >> %s" % (inputF, inputF2, ncpu, outputLogSummary, outputLogSummary, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Screen.seqs
    step = 10
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1]  # align
    inputF2 = command_meta[step-1][0]["output_files"][1] # summary file
    inputF3 = command_meta[step-1][0]["output_files"][2] # summary log file
    inputF4 = command_meta[step-3][0]["output_files"][1] # counts file
    input_files = [inputF, inputF2, inputF3, inputF4]
    outputF_good_fasta = "%s/combo_fastq.trim.contigs.good.unique.good.align" % (output_dir)
    outputF_bad_acc_nos = "%s/combo_fastq.trim.contigs.good.unique.bad.accnos" % (output_dir)
    outputF_count_table = "%s/combo_fastq.trim.contigs.good.good.count_table" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_good_fasta, outputF_bad_acc_nos, outputF_count_table]
    cmd = "echo \"screen.seqs( fasta=%s,start=$(grep \"Median\" %s | tr -d \'\\n\' | awk {\'print int($2)\'}) ,end=$(grep \"Median\" %s | tr -d \'\\n\' | awk {\'print int($3)\'}) ,maxhomop=8 ,count=%s ,summary=%s, processors=%s )\" | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF3, inputF3, inputF4, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.seqs
    step = 11
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # fasta
    inputF2 = command_meta[step-1][0]["output_files"][3] # count table
    input_files = [inputF]
    outputF_summary = "%s/combo_fastq.trim.contigs.good.unique.good.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary]
    cmd = "echo 'summary.seqs( fasta=%s, count=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # count.seqs
    step = 12
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1]
    input_files = [inputF]
    outputF_filter = "%s/combo_fastq.filter" % (output_dir)
    outputF_fasta = "%s/combo_fastq.trim.contigs.good.unique.good.filter.fasta" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_filter, outputF_fasta]
    cmd = "echo 'filter.seqs( fasta=%s, vertical=true, trump=., soft=0, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Unique.seqs
    step = 13
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][2] # counts fasta
    inputF2 = command_meta[step-3][0]["output_files"][3] # counts table
    input_files = [inputF]
    outputF_unique_fasta = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.fasta" % (output_dir)
    outputF_unique_table = "%s/combo_fastq.trim.contigs.good.unique.good.filter.count_table" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_unique_fasta, outputF_unique_table]
    cmd = "echo 'unique.seqs( count=%s, fasta=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF2, inputF, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # precluster.seqs
    step = 14
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # fasta
    inputF2 = command_meta[step-1][0]["output_files"][2] # table
    input_files = [inputF, inputF2]
    outputF_unique_fasta = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.fasta" % (output_dir)
    outputF_unique_table = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.count_table" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_unique_fasta, outputF_unique_table]
    cmd = "echo 'pre.cluster( fasta=%s, count=%s, diffs=2, match=1, mismatch=-1, gapopen=-2, gapextend=-1, topdown=true, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # vchimera.seqs
    step = 15
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # fasta
    inputF2 = command_meta[step-1][0]["output_files"][2] # table
    input_files = [inputF, inputF2]
    outputF_chimera_search = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.chimeras" % (output_dir)
    outputF_chimera_table = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table" % (output_dir)
    outputF_chimera_accnos = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_chimera_search, outputF_chimera_table, outputF_chimera_accnos]
    cmd = "echo 'chimera.vsearch( fasta=%s, reference=self, abskew=1.9, count=%s, minh=0.3, mindiv=0.5, xn=8.0, dn=1.4, dereplicate=true, mindiffs=3, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})


    # remove.seqs
    step = 16
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1] # fasta
    inputF2 = command_meta[step-1][0]["output_files"][2] # table
    inputF3 = command_meta[step-1][0]["output_files"][3] # accnos
    input_files = [inputF, inputF2, inputF3]
    outputF_unique_fasta = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta" % (output_dir)
    outputF_unique_table = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_unique_fasta, outputF_unique_table]
    cmd = "echo 'remove.seqs( accnos=%s, fasta=%s,count=%s  )' | sed 's/ //g' | mothur | tee -a %s" % (inputF3, inputF, inputF2, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.seqs
    step = 17
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # fasta
    inputF2 = command_meta[step-1][0]["output_files"][2] # count table
    input_files = [inputF, inputF2]
    outputF_summary = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary]
    cmd = "echo 'summary.seqs( fasta=%s, count=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # classify.seqs
    step = 18
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1] # fasta
    inputF2 = command_meta[step-2][0]["output_files"][2] # count table
    input_files = [inputF, inputF2]
    outputF_taxonomy = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.taxonomy" % (output_dir)
    outputF_tax_summary = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.tax.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_taxonomy, outputF_tax_summary]
    cmd = "echo 'classify.seqs( fasta=%s, reference=%s, taxonomy=%s, method=wang, ksize=8, iters=100, cutoff=80, probs=true, count=%s, relabund=false, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, options.taxonomy_alignment_ref, options.taxonomy_ref, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # remove.lineage
    step = 19
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # tax seq
    inputF2 = command_meta[step-3][0]["output_files"][1] # fasta
    inputF3 = command_meta[step-3][0]["output_files"][2] # count
    input_files = [inputF, inputF2, inputF3]
    outputF_taxonomy = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.taxonomy" % (output_dir)
    outputF_fasta = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta" % (output_dir)
    outputF_table = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_taxonomy, outputF_fasta, outputF_table]
    cmd = "echo 'remove.lineage( taxonomy=%s ,taxon=\"Chloroplast-Mitochondria-unknown-Archaea-Eukaryota\" ,fasta=%s ,count=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, inputF3, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.tax
    step = 20
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # taxonomy
    inputF2 = command_meta[step-1][0]["output_files"][3] # count table
    input_files = [inputF, inputF2]
    outputF_summary = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tax.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary]
    cmd = "echo 'summary.tax( taxonomy=%s, count=%s, relabund=false )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # Summary.seqs
    step = 21
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][2] # fasta
    inputF2 = command_meta[step-2][0]["output_files"][3] # count table
    input_files = [inputF, inputF2]
    outputF_summary = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_summary]
    cmd = "echo 'summary.seqs( fasta=%s, count=%s, processors=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, ncpu, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # phylotype
    step = 22
    command_meta[step] = []
    inputF = command_meta[step-3][0]["output_files"][1] # taxonomy
    input_files = [inputF]
    outputF_rank = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.rabund" % (output_dir)
    outputF_species = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.sabund" % (output_dir)
    outputF_otu = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.list" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_rank, outputF_species, outputF_otu]
    cmd = "echo 'phylotype( taxonomy=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # count groups
    step = 23
    command_meta[step] = []
    inputF = command_meta[step-4][0]["output_files"][3] # count table
    input_files = [inputF]
    outputF_count = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_count]
    cmd = "echo 'count.groups( count=%s )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # classify otu
    step = 24
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][3] # otu list
    inputF2 = command_meta[step-5][0]["output_files"][3] # count table
    inputF3 = command_meta[step-5][0]["output_files"][1] # taxonomy
    input_files = [inputF, inputF2, inputF3]
    outputF_taxonomy = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.1.cons.taxonomy" % (output_dir)   ### to s3
    outputF_summary = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.1.cons.tax.summary" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_taxonomy, outputF_summary]
    cmd = "echo 'classify.otu( list=%s, taxonomy=%s, count=%s, label=1, basis=otu, probs=false, persample=false, cutoff=60 )' | sed 's/ //g' | mothur | tee -a %s; cp %s %s" % (inputF, inputF3, inputF2, outputLog, outputF_taxonomy, options.output_classify_otu)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # make.shared
    step = 25
    command_meta[step] = []
    inputF = command_meta[step-3][0]["output_files"][3] # otu list
    inputF2 = command_meta[step-6][0]["output_files"][3] # count table
    input_files = [inputF, inputF2]
    outputF_shared = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.shared" % (output_dir)   ### to s3
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_shared]
    cmd = "echo 'make.shared( count=%s, label=1, list=%s )' | sed 's/ //g' | mothur | tee -a %s; cp %s %s" % (inputF2, inputF, outputLog, outputF_shared, options.output_make_shared)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    if not options.input_dir:
        # s3 transfer
        #step = 25
        #command_meta[step] = []
        inputF = command_meta[step-1][0]["output_files"][1] # taxonomy - otu
        input_files = [inputF]
        outputF = "%s/s3_transfer.log" % (output_dir)
        output_files = [outputF]
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'Taxonomy.txt'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'file'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
        command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # make.biom
    step = 26
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # shared
    inputF2 = command_meta[step-2][0]["output_files"][1] # taxonomy
    input_files = [inputF, inputF2]
    outputF_biom = "%s/combo_fastq.trim.contigs.good.unique.good.filter.unique.precluster.pick.2.wang.pick.tx.1.biom" % (output_dir)   
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF_biom]
    cmd = "echo 'make.biom( shared=%s, constaxonomy=%s, matrixtype=sparse )' | sed 's/ //g' | mothur | tee -a %s" % (inputF, inputF2, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    if not options.input_dir:
        # s3 transfer
        #step = 26
        #command_meta[step] = []
        inputF = command_meta[step-1][0]["output_files"][1] # taxonomy - otu
        input_files = [inputF]
        outputF = "%s/s3_transfer.log" % (output_dir)
        output_files = [outputF]
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'Shared.txt'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'file'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
        command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # RCODE
    step = 27
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][1] # shared
    inputF2 = command_meta[step-3][0]["output_files"][1] # taxonomy
    input_files = [inputF, inputF2]
    outputF = "%s/otu_table.txt" % (output_dir)  #### to s3
    outputDirTmp = "%s/otu_tmp" % (output_dir)
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF]
    cmd = "echo python /opt/galaxy/tools/colgate/R_code_otu_table_wrapper.py -i %s -t %s -o %s --out-dir %s | tee -a %s; python /opt/galaxy/tools/colgate/R_code_otu_table_wrapper.py -i %s -t %s -o %s --out-dir %s | tee -a %s; cp %s %s" % (inputF, inputF2, outputF, output_dir, outputLog, inputF, inputF2, outputF, outputDirTmp, outputLog, outputF, options.output_otu_table)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})


    # KRONA
    step = 28
    command_meta[step] = []
    inputF = command_meta[step-3][0]["output_files"][1] # shared
    inputF2 = command_meta[step-4][0]["output_files"][1] # taxonomy
    input_files = [inputF, inputF2]
    outputF = "%s/text.krona.html" % (output_dir)   ### to s3
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)
    output_files = [outputLog, outputF]
    if not options.input_dir:
        cmd = "cp %s %s; cd %s; echo python /opt/galaxy/tools/krona/mycrobiota.py --command create_krona_plot_multisample --shared_file %s --level 1 --taxonomy %s -od %s | tee -a %s; python /opt/galaxy/tools/krona/mycrobiota.py --command create_krona_plot_multisample --shared_file %s --level 1 --taxonomy %s -od %s| tee -a %s; cd -; cp %s %s" % (local_cred_file_out, cred_file_loc_out, output_dir, inputF, inputF2, output_dir, outputLog, inputF, inputF2, output_dir, outputLog, outputF, options.output_krona_html)
    else:
        cmd = "cd %s; echo python /opt/galaxy/tools/krona/mycrobiota.py --command create_krona_plot_multisample --shared_file %s --level 1 --taxonomy %s -od %s | tee -a %s; python /opt/galaxy/tools/krona/mycrobiota.py --command create_krona_plot_multisample --shared_file %s --level 1 --taxonomy %s -od %s| tee -a %s; cd -; cp %s %s" % (output_dir, inputF, inputF2, output_dir, outputLog, inputF, inputF2, output_dir, outputLog, outputF, options.output_krona_html)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    if not options.input_dir:
        # s3 transfer
        #step = 28
        #command_meta[step] = []
        inputF = command_meta[step-1][0]["output_files"][1] # taxonomy - otu
        input_files = [inputF]
        outputF = "%s/s3_transfer.log" % (output_dir)
        output_files = [outputF]
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'OTU_table.txt'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'file'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
        command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    # core diversity
    step = 29
    command_meta[step] = []
    #inputF = command_meta[step-28][0]["output_files"][1] # sample list
    inputF = options.mapping_file
    inputF2 = command_meta[step-3][0]["output_files"][1] # biom
    inputF3 = command_meta[step-6][0]["output_files"][1] # count groups
    inputF4 = ""
    if options.categories is not None:
        inputF4 = "--categories %s" % options.categories
    input_files = [inputF, inputF2]
    outputDir = "%s/core_diversity_analyses" % (output_dir)  #### to s3
    outputF = "%s/index.html" % (outputDir) ###### to s3
    outputMatplotlib = "%s/matplotlibrc" % (output_dir) 
    #outputLog = "%s/%s.mother.out.log" % (output_dir, sample_name)  ### to s3
    output_files = [outputLog, outputF, outputDir]
    cmd = "cd %s; echo backend:agg > %s; value=`awk \'{print $2}\' %s | sort -n | head -1`; core_diversity_analyses.py --input_biom_fp %s -o %s --mapping_fp %s --sampling_depth $value --nonphylogenetic_diversity %s | tee -a %s; cd -; cp %s %s; cp -r %s %s" % (output_dir, outputMatplotlib, inputF3, inputF2, outputDir, inputF, inputF4, outputLog, outputF, options.output_core_diversity_html, outputDir, options.output_core_diversity_dir)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    if not options.input_dir:
        #step = 29
        #command_meta[step] = []
        inputF = command_meta[step-1][0]["output_files"][1] # taxonomy - otu
        input_files = [inputF]
        outputF = "%s/s3_transfer.log" % (output_dir)
        output_files = [outputF]
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'Krona_HTML.html'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'file'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
        command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})


    step = 30
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][1] # taxonomy - otu
    input_files = [inputF]
    outputF = "%s/s3_transfer.log" % (output_dir)
    output_files = [outputF]
    if not options.input_dir:
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'core_diversity.html'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'file'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
    else:
        cmd = "echo hello"
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    step = 31
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_files"][2] # core diversity dir
    input_files = [inputF]
    outputF = "%s/s3_transfer.log" % (output_dir)
    output_files = [outputF]
    if not options.input_dir:
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'core_diversity'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'dir'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
    else:
        cmd = "echo hello"
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})


    step = 32
    command_meta[step] = []
    inputF = outputLog # taxonomy - otu
    input_files = [inputF]
    outputF = "%s/s3_transfer.log" % (output_dir)
    output_files = [outputF]
    if not options.input_dir:
        transfer_info_curr = transfer_info_out
        transfer_info_curr['rename'] = 'Log.txt'
        #transfer_info_curr['sse'] = True
        transfer_info_curr['from_path'] = inputF
        transfer_info_curr['object_type'] = 'file'
        cmd = "cp %s %s; python /opt/galaxy/tools/globus/s3_transfer_out.py --transfer-info \"%s\"" % (local_cred_file_out, cred_file_loc_out, transfer_info_curr)
    else:
        cmd = "echo hello"
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #### testing up to this point. Look at "/scratch/galaxy/tmp/optimizedqPyXLR" for latest test run
    # hopefully that run produces non-empty accnos file.

    return command_meta


def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option("--transfer-info", dest="transfer_info", help="Detailed Transfer Information")
    parser.add_option("--transfer-info-out", dest="transfer_info_out", help="Detailed Transfer Information")
    parser.add_option("--input-dir", dest="input_dir", help="only for debugging")
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '--align-ref', dest="align_reference", help='The align reference to use or index' )
    parser.add_option( '--taxonomy-ref', dest="taxonomy_ref", help='The taxonomy reference to use or index' )
    parser.add_option( '--taxonomy-align-ref', dest="taxonomy_alignment_ref", help='The taxonomy alignment reference to use or index' )
    parser.add_option( '', '--mapping-file', dest="mapping_file", help='The path to fastq files to use ' )
    parser.add_option( '', '--categories', dest="categories", help="Categories for core_diversity tool" )
    parser.add_option( '--output-classify-otu', dest="output_classify_otu", help="the file to save")
    parser.add_option( '--output-make-shared', dest="output_make_shared", help="the file to save")
    parser.add_option( '--output-otu-table', dest="output_otu_table", help="the file to save")
    parser.add_option( '--output-krona-html', dest="output_krona_html", help='The file to save' )
    parser.add_option( '--output-core-diversity-html', dest="output_core_diversity_html", help='The file to save' )
    parser.add_option( '--output-core-diversity-directory', dest="output_core_diversity_dir", help='The dir output' )
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )
    (options, args) = parser.parse_args()

    final_outputs = []
    #final_outputs = [options.output_contatenate_fastq, options.output_tophat_bam, options.output_trinity_transcripts, options.output_picard_markDuplicates_html, options.output_picard_markDuplicates_directory, options.output_picard_alignmentSummary, options.output_picard_alignmentSummary_directory, options.output_picard_collectRnaSeqSummary, options.output_picard_collectRnaSeqSummary_directory, options.output_tophatStats, options.output_htseq_counts, options.output_htseq_log, options.output_fastqc_html, options.output_fastqc_text, options.output_log ]


    if len(args) > 0:
        parser.error('Wrong number of arguments')

#    # make temp directory 
    output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized")
    tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimizedtmp")
    # make temp directory 
#    output_dir = tempfile.mkdtemp(prefix="optimized")
#    tmp_dir = tempfile.mkdtemp(prefix="optimizedtmp")

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir)

    #print wf_meta
    os.chdir(output_dir)
    for step_id, step_list in wf_meta.iteritems():
        print "START STEP %s: %s" % (step_id, workflows.time_stamp())
        step_input_files = []
        if len(step_list) == 1:
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                workflows.run_cmd ( job['cl'], tmp_dir, "running some job")
                for ifile in job['input_files']:
                    if ifile not in step_input_files:
                        step_input_files.append(ifile)
        else:   # run jobs in parallel
            ps = {}
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                p = workflows.run_cmd_parallel ( job['cl'], tmp_dir, "running some job")
                ps[p.pid] = p
            print "Waiting for %d processes..." % len(ps)
            while ps:
                pid, status = os.wait()
                if pid in ps:
                    del ps[pid]
                    print "Waiting for %d processes..." % len(ps)
        print "END STEP %s: %s" % (step_id, workflows.time_stamp())


        # job completed, check to see if there are inputs that should be deleted
        if step_id > 0:
            keep_files = workflows.job_step_dependency(wf_meta, step_id)
            #workflows.delete_file(step_input_files, keep_files, step_id, final_outputs)


    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
