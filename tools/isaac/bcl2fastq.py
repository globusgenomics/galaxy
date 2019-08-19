#!/usr/bin/env python
import time
from datetime import datetime, timedelta
import optparse, os, shutil, subprocess, sys, tempfile, re
from bioblend import galaxy
import requests
from xml.dom.minidom import parse
import xml.dom.minidom
requests.packages.urllib3.disable_warnings()

def copy_samplesheet(original, dest):
    # instead of simply copying the samplesheet we need to parse the original file
    # because bcl2fastq does not recognize ^M characters as new lines
    original_ss = open(original, 'rU')
    dest_ss = open(dest, 'w')
    for line in fh_config:
        dest_ss.write("%s\n" % line.rstrip("\n"))
    dest_ss.close()
    original_ss.close()

def modifySamplesheet(samplesheet, config_file, tempdir, action):
    samples_dict = {}
    samples_list = []
    fh_config = open(config_file, 'rU')
    max_index_len = 0
    flag = 0
    for line in fh_config:
        if line.startswith("#") or line.startswith("Lane,"):
            flag = 1
        if line.startswith("#") or line.startswith("Lane,") or flag == 0:
            continue
        line = line.rstrip("\n")
        values = line.split(",")
        #samples_dict[values[1]] = {'type':values[4], 'vendor':values[3], 'lane':values[0], 'index':values[2]}
        #samples_list.append({'name':values[1], 'type':values[4], 'vendor':values[3], 'lane':values[0], 'index':values[2]})
        #length = values[2]
        samples_dict[values[1]] = {'type':values[5], 'vendor':values[6], 'lane':values[0], 'index':values[3]}
        samples_list.append({'name':values[1], 'type':values[5], 'vendor':values[6], 'lane':values[0], 'index':values[3]})
        length = values[3]
        if len(length) > max_index_len:
            max_index_len = len(length)
    
    # write the new samplesheet
    new_samplesheet = "%s/NewSampleSheet.csv" % tempdir
    fh_new_ss = open(new_samplesheet, 'w')
    fh_old_ss = open(samplesheet, 'rU')
    line_flag = False
    for line in fh_old_ss:
        if line_flag is False and not line.startswith("Lane"):
            fh_new_ss.write(line)
        elif line_flag is False and line.startswith("Lane"):
            fh_new_ss.write(line)
            line_flag = True
        elif line_flag is True:
            if action == "modify":
                line = line.rstrip("\n")
                values = line.split(",")
                values[-1] = "%s/%s" % (samples_dict[values[-1]]['type'], values[-1])
                fh_new_ss.write(",".join(values) + "\n")
            elif action == "create":
                for sample in samples_list:
                    #if len(samples[sample]['index']) != max_index_len:
                    #    diff = max_index_len - len(samples[sample]['index'])
                    #    samples[sample]['index'] = "%s%s" % (samples[sample]['index'], "N"*diff)
                    line = "%s,%s,%s,%s,%s/%s\n" % (sample['lane'], sample['name'], sample['name'], sample['index'],sample['type'], sample['name'])
                    fh_new_ss.write(line)
    fh_new_ss.close()
    fh_old_ss.close()

    return new_samplesheet,max_index_len

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    descr = "bcl2fastq_wrapper.py: Convert BCL to Fastq files"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-t', '--input-tar', dest="input_tar", help="Path to Illumina tar file" )
    parser.add_option( '', '--input-name', dest="input_name", help="name of input tar" )
    parser.add_option( '-i', '--input-dir', dest="input_dir", help="Path to Illumina directory data" )
    parser.add_option( '-s', '--samplesheet', dest="samplesheet", help="User provided samplesheet for demultiplexing" )
    parser.add_option( '', '--no-index', dest="non_index_flag", action="store_true", help="Do not run any indexes in bcl2fastq" )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly without any modification.' )
    parser.add_option( '', '--input-config', dest='config_file', help="File containing relationship of sample names to vendor and data type" )
    parser.add_option( '', '--reports-out', dest='html_reports', help='')
    parser.add_option( '', '--reports-out-dir', dest='reports_outdir', help='')
    parser.add_option( '', '--interop-out', dest='html_interop', help='')
    parser.add_option( '', '--interop-out-dir', dest='interop_outdir', help='')
    parser.add_option( '', '--stats-out', dest='html_stats', help='')
    parser.add_option( '', '--stats-out-dir', dest='stats_outdir', help='')
    parser.add_option( '', '--fastq-out', dest='html_fastq', help='')
    parser.add_option( '', '--fastq-out-dir', dest='fastq_outdir', help='')
    parser.add_option( '--userkey', dest='user_api_key', help="The user api key" )
    parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
    (options, args) = parser.parse_args()

    # create all output dirs
    if not os.path.exists(options.fastq_outdir):
        os.mkdir(options.fastq_outdir)
    if not os.path.exists(options.stats_outdir):
        os.mkdir(options.stats_outdir)
    if not os.path.exists(options.interop_outdir):
        os.mkdir(options.interop_outdir)
    if not os.path.exists(options.reports_outdir):
        os.mkdir(options.reports_outdir)

    # Untar if necessary
    illumina_dir = None
    tmpdir = None
    runinfo_path = None
    if options.input_dir:
        # assume directory is owned by Galaxy user
        # if Samplesheet.csv is provided it will need to be integrated into the directory
        illumina_dir = options.input_dir
    else:
        # untar
        #tmpdir = tempfile.mkdtemp( prefix='tmp-illumina-dir')
        tmpdir = "%s/tmp-illumina-dir" % options.fastq_outdir
        illumina_dir = tmpdir
        if not os.path.exists(illumina_dir):
            os.mkdir(illumina_dir)
        #tar_name = os.path.basename(options.input_tar)
        #project_name = os.path.splitext(tar_name)[0]
        cmd = "tar xvf %s -C %s" % ( options.input_tar, tmpdir )

        #set up stdout and stderr output options
        stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
        stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name

        print "...Untarring input Illumina Tar file: %s" % cmd
        buffsize = 1048576
        proc = subprocess.Popen( args=cmd, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True )

        return_code = proc.wait()
        if return_code:
            tmp_stderr = open( stderr_name, 'rb' )
            stderr = ''
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            if return_code != 0:
                raise Exception, stderr

        if options.non_index_flag:
            # modify the runParameters.xml file
            # look for line:    <IndexRead1>8</IndexRead1>
            # and change it to: <IndexRead1>0</IndexRead1>
            runparameters_path = "%s/runParameters.xml" % illumina_dir 
            before = open(runparameters_path, "r")
            post_file = runparameters_path + ".post"
            f = open(post_file, "w")
            for line in before:
                if "IndexRead1" in line:
                    line = "    <IndexRead1>0</IndexRead1>\n"
                f.write(line)
            f.close()
            before.close()
            shutil.move(post_file, runparameters_path)

            # modify the RunInfo.xml file
            # look for line:    <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
            # and change it to: <Read Number="2" NumCycles="0" IsIndexedRead="N" />
            runinfo_path = "%s/RunInfo.xml" % illumina_dir
            samplesheet_path = "%s/SampleSheet.csv" % illumina_dir
            before = open(runinfo_path, "r")
            post_file = runinfo_path + ".post"
            f = open(post_file, "w")
            for line in before:
                if '<Read Number="2"' in line:
                    line = '      <Read Number="2" NumCycles="0" IsIndexedRead="N" />\n'
                f.write(line)
            f.close()
            before.close()
            shutil.move(post_file, runinfo_path)
            
            # delete the SampleSheet.csv file
            os.remove("%s/SampleSheet.csv" % illumina_dir )
	else:
            # Copy Samplesheet.csv into the directory structure
            #illumina_dir = "%s/%s" % (tmpdir, project_name)
            dest = "%s/SampleSheet.csv" % illumina_dir
            newsamplesheet = None
            index_length = 8
            if options.config_file is not None and "createSamplesheet_template" not in options.samplesheet:
                (newsamplesheet, index_length) = modifySamplesheet(options.samplesheet, options.config_file, tmpdir, "modify")
            elif options.config_file is not None and "createSamplesheet_template" in options.samplesheet:
                (newsamplesheet,index_length) = modifySamplesheet(options.samplesheet, options.config_file, tmpdir, "create")
            else: 
                newsamplesheet = options.samplesheet
            if newsamplesheet is not None:
                copy_samplesheet(newsamplesheet, dest)
                #shutil.copyfile(newsamplesheet, dest)
            runinfo_path = "%s/RunInfo.xml" % illumina_dir
            samplesheet_path = "%s/SampleSheet.csv" % illumina_dir
            
            # modify RunInfo file if necessary
            if index_length != 8:
                before = open(runinfo_path, "r")
                post_file = runinfo_path + ".post"
                f = open(post_file, "w")
                for line in before:
                    if '<Read Number="2"' in line:
                        line = '      <Read Number="2" NumCycles="6" IsIndexedRead="Y" />\n'
                    f.write(line)
                f.close()
                before.close()
                shutil.move(post_file, runinfo_path)


    # Get RunInfo Flowcell id
    flowcell_id = None
    if runinfo_path is not None:
        DOMTree = xml.dom.minidom.parse(runinfo_path)    
        context = DOMTree.documentElement
        flowcell_id = context.getElementsByTagName("Flowcell")[0].childNodes[0].data

        # look if the samples in the samplesheet require creation of sub-directories:
        # i.e. sample name column is: exome/sample_name1
        ss_fh = open(samplesheet_path, "r")
        seen_flag = 0
        sub_dirs = {}
        for line in ss_fh:
            if seen_flag == 1:
                values = line.split(",")
                sample_sub_dirs = values[-1].split("/")
                if len(sample_sub_dirs) > 1:
                    sub_dirs[sample_sub_dirs[0]] = 1
                    continue
            if line.startswith("Lane"):
                seen_flag = 1

        # if flowcell_id is not None then create the necessary directories in the stats output directory
        if len(sub_dirs) > 0:
            for sub_dir in sub_dirs:
                os.makedirs("%s/html/%s/%s" % (options.reports_outdir, flowcell_id, sub_dir) )

    # Construct BCL2fastq command
    cmd = "bcl2fastq --minimum-trimmed-read-length 50 --runfolder-dir %s --output-dir %s --stats-dir %s --reports-dir %s --interop-dir %s %s" % (illumina_dir, options.fastq_outdir, options.stats_outdir, options.reports_outdir, options.interop_outdir, ' '.join( options.pass_through_options ))
    print "...Initiating bcl2fastq: %s" % cmd

    #set up stdout and stderr output options
    stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
    stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name

    buffsize = 1048576
    proc = subprocess.Popen( args=cmd, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True )
    return_code = proc.wait()                 
    if return_code:
        tmp_stderr = open( stderr_name, 'rb' )
        stderr = ''                           
        try:
            while True:
                stderr += tmp_stderr.read( buffsize ) 
                if not stderr or len( stderr ) % buffsize != 0: 
                    break
        except OverflowError:                 
            pass                                
        tmp_stderr.close()                      
        if return_code != 0:                    
            raise Exception, stderr
    #shutil.rmtree(tmpdir, ignore_errors=True)

    # create individual directories for each lane in the output directory
#-rw-r--r-- 1 galaxy galaxy 31765680039 Mar 15 21:47 Undetermined_S0_L008_R1_001.fastq.gz
#-rw-r--r-- 1 galaxy galaxy  1368217453 Mar 15 21:47 Undetermined_S0_L008_R2_001.fastq.gz
#-rw-r--r-- 1 galaxy galaxy 33953821320 Mar 15 21:47 Undetermined_S0_L008_R3_001.fastq.gz

    if options.non_index_flag:
        prefix_fastqs = "Undetermined_S0_L00"
        postfix_fastqs = "_001.fastq.gz"
        for i in range(1, 9):
           lane_prefix = "%s%s" % (prefix_fastqs, i) # "Undetermined_S0_L008"
           lane_dir = "%s/%s" % (options.fastq_outdir, lane_prefix)  # "/output_path/Undetermined_S0_L008"
           lane_read1 = "%s_R1%s" % (lane_dir, postfix_fastqs)   # "/output_path/Undetermined_S0_L008_R1_001.fastq.gz"
           lane_read2 = "%s_R2%s" % (lane_dir, postfix_fastqs)   # "/output_path/Undetermined_S0_L008_R2_001.fastq.gz"
           lane_read3 = "%s_R3%s" % (lane_dir, postfix_fastqs)   # "/output_path/Undetermined_S0_L008_R3_001.fastq.gz"

           if os.path.exists(lane_read1) and os.path.exists(lane_read2) and os.path.exists(lane_read3):
               # move read1 to lane directory, delete read2 (index file), rename read3 to read2 into lane directory
               if not os.path.exists(lane_dir):
                   os.mkdir(lane_dir)
               shutil.move(lane_read1, lane_dir)
               os.remove(lane_read2)
               shutil.move(lane_read3, "%s/%s_R2%s" % (lane_dir, lane_prefix, postfix_fastqs) )

    # print output files
    output_files = {options.html_fastq : options.fastq_outdir, options.html_stats : options.stats_outdir, options.html_reports : options.reports_outdir, options.html_interop : options.interop_outdir}
    for filename, directory in output_files.iteritems():
        f = open(filename,'w')
        for dirName, subdirList, fileList in os.walk(directory):
            f.write('\nFound directory: %s\n' % dirName)
            for fname in fileList:
                f.write('\t%s\n' % fname)
            # Remove the first entry in the list of sub-directories
            # if there are any sub-directories present
            if len(subdirList) > 0:
                del subdirList[0]
        f.close()

    #add to shared library if necessary
    if options.user_api_key:
        lib_name = options.input_name
        url = options.galaxyurl
        if "http:" in url:
            url = url.replace("http", "https")

        gi =  galaxy.GalaxyInstance(url=url, key=options.user_api_key)
        library = gi.libraries.create_library(lib_name, description="Results from bcl2fastq run")
        gi.libraries.upload_from_galaxy_filesystem(library["id"],options.fastq_outdir,file_type='fastqsanger',dbkey='hg19',link_data_only="True")

        # add the Stats output directory
        folder = gi.libraries.create_folder(library["id"], "Stats")
        gi.libraries.upload_form_galaxy_filesystem(library["id"],options.stats_outdir,folder_id=folder['id'],file_type='xml',dbkey='hg19',link_data_only="True")

        # add the interop output directory
        folder = gi.libraries.create_folder(library["id"], "InterOp")
        gi.libraries.upload_form_galaxy_filesystem(library["id"],options.interop_outdir,folder_id=folder['id'],file_type='xml',dbkey='hg19',link_data_only="True")

    # create symlinks to the fastq files so that isaac-align can be run
    directory =  options.fastq_outdir
    for dirName, subdirList, fileList in os.walk(directory):
        for fname in fileList:
           matchObj = re.match(r'.*_L00(.*)_R(.*)_.*', fname, re.M|re.I)
           if matchObj:
               lane = matchObj.group(1)
               read = matchObj.group(2)
               link_name = "%s/lane%s_read%s.fastq.gz" % (dirName, lane, read)
               os.symlink("%s/%s" % (dirName,fname), link_name)

if __name__ == "__main__":
    __main__()
