#!/usr/bin/python

"""
Runs wfdb commands on physionet db files
wfdb_wrapper.py --program sqrs --record /scratch/go/galaxy/cvrg/100.hea

usage: wfdb_wrapper.py [options]
   -p, --program wfdb_command_name
   -r, --record /filepath/to/record.hea
   -d, --datafile /filepath/to/record.dat
   -a, --annotator /filepath/to/record.[qrs,wqrs,atr]
       --inAnntype [qrs,wqrs,atr]
   -b, --begintime begin_time
   -e, --stoptime stop_time
       --resolution [true|false]
       --jpoints [true|false]
       --frequency frequency
       --resample [true|false]
       --insertQann M
       --outAnnType [annotator type]
       --allIntervals $selectOptions.allIntervals
       --consecutiveIntervals $selectOptions.consecutiveIntervals
       --iformat $options.iformat
       --itypeEnd $options.itypeEnd
       --itypeBegin $options.itypeBegin
       --timesFormatEnd $options.timesFormatEnd
       --timesFormatBegin $options.timesFormatBegin
       --annotationsEnd $options.annotationsEnd
       --annotationsBegin $options.annotationsBegin
     --intervals $options.intervals
     --precision $options.selectTimePrecision
     --timeFormat $options.selectTimeFormat
     --signalList $options.signalList
     --signalStart $options.signalStart
     --printHeader $options.printHeader
     --outputFormat $options.outputFormat
   -t, --threshold threshold
   -s, --signal signal
   -o, --outputfile outputfilename

list of commands:
sqrs
wqrs
nguess

"""

import optparse, os, shutil, sys, tempfile, re, subprocess, time

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def get_recordname( record_file ):
    recordFH = open(record_file, 'r')
    firstLine = recordFH.readline()
    recordFH.close()
    metadata = firstLine.rsplit(" ")
    record_name = metadata[0]
    print record_name
    return record_name


def __main__():
    # timestamp
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))    

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--program', dest='program', help='wfdb program to use' )
    parser.add_option( '-r', '--record', dest='recordfile', help='/filepath/to/record.hea' )
    parser.add_option( '-d', '--datafile', dest='datafile', help='/filepath/to/record.dat' )
    parser.add_option( '--annotator', dest='annotatorfile', help='/filepath/to/record.[qrs|wqrs|atr|nguess]' )
    parser.add_option( '-q', '--qrsfile', dest='qrsfile', help='/filepath/to/record.qrs' )
    parser.add_option( '-a', '--atrfile', dest='atrfile', help='/filepath/to/record.atr' )
    parser.add_option( '--inAnnType', dest='inAnnType', help='input annotator type [wqrs, qrs, etc]' )
    parser.add_option( '--outAnnType', dest='outAnnType', help='output annotator type [wqrs, qrs, etc]' )
    parser.add_option( '-b', '--begintime', dest='begintime', help='First time period to make calculations')
    parser.add_option( '-e', '--stoptime', dest='stoptime', help='Last time period to make calculations')
    parser.add_option( '--resolution', dest='resolution', help='resolution definition level')
    parser.add_option( '--jpoints', dest='jpoints', help='find and annotate J-points (QRS ends) as well as QRS onsets')
    parser.add_option( '--frequency', dest='frequency', help='specify power line (mains) frequency ')
    parser.add_option( '--resample', dest='resample', help='resample input')
    parser.add_option( '--insertQann', dest='insertQann', help='insert a Q if RR > M ')
    parser.add_option( '-s', '--signal', dest='signal', help='Analyze specified signal')
    parser.add_option( '-t', '--threshold', dest='threshold', help='Threshold')
    parser.add_option( '-o', '--outputfile', dest='outputfile', help='Output file')
    parser.add_option( '--allIntervals', dest='allIntervals', help='intervals [true|fals]')
    parser.add_option( '--consecutiveIntervals', dest='consecutiveIntervals', help='print consecutive intervals [true|false]')
    parser.add_option( '--iformat', dest='iformat', help='input format')
    parser.add_option( '--itypeEnd', dest='itypeEnd', help='input type end')
    parser.add_option( '--itypeBegin', dest='itypeBegin', help='input type begin')
    parser.add_option( '--timesFormatEnd', dest='timesFormatEnd', help='format time at end')
    parser.add_option( '--timesFormatBegin', dest='timesFormatBegin', help='format time at beginning')
    parser.add_option( '--annotationsEnd', dest='annotationsEnd', help='annotations at end')
    parser.add_option( '--annotationsBegin', dest='annotationsBegin', help='annotations at beginning')
    parser.add_option( '--intervals', dest='intervals', help='rdsam param' )
    parser.add_option( '--precision', dest='selectTimePrecision', help='rdsam param' )
    parser.add_option( '--timeFormat', dest='selectTimeFormat', help='rdsam param' )
    parser.add_option( '--signalList', dest='signalList', help='rdsam param' )
    parser.add_option( '--signalStart', dest='signalStart', help='rdsam param' )
    parser.add_option( '--printHeader', dest='printHeader', help='rdsam param' )
    parser.add_option( '--outputFormat', dest='outputFormat', help='rdsam param' )
    (options, args) = parser.parse_args()

    # set environment variables
    #os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":/mnt/galaxyTools/tools/wfdb/10.5.20/lib64"
    os.environ["LD_LIBRARY_PATH"] = "/mnt/galaxyTools/tools/wfdb/10.5.20/lib64"

    # get record name
    # open the record file and read the first line.
    if options.recordfile:
        record_name = get_recordname(options.recordfile)
    else:
        print "no record file provided"
        sys.exit()

    # create a temporary working directory
    try:
        tmp_wfdb_dir = tempfile.mkdtemp(dir="/tmp")

        #cp any files to this temporary directory
        if options.recordfile:
            #tmp_record_filename = os.path.basename(options.recordfile)
            tmp_record_filename = "%s.hea" % (record_name)
            dest = "%s/%s" % (tmp_wfdb_dir, tmp_record_filename)
            shutil.copyfile(options.recordfile, dest)
            
        if options.qrsfile:
            #tmp_qrs_filename = os.path.basename(options.qrsfile)
            tmp_qrs_filename = "%s.qrs" % (record_name)
            dest = "%s/%s" % (tmp_wfdb_dir, tmp_qrs_filename)
            shutil.copyfile(options.qrs, dest)

        if options.atrfile:
            #tmp_atr_filename = os.path.basename(options.atrfile)
            tmp_atr_filename = "%s.atr" % (record_name)
            dest = "%s/%s" % (tmp_wfdb_dir, tmp_atr_filename)
            shutil.copyfile(options.atr, dest)

        if options.datafile:
            #tmp_data_filename = os.path.basename(options.datafile)
            tmp_data_filename = "%s.dat" % (record_name)
            dest = "%s/%s" % (tmp_wfdb_dir, tmp_data_filename)
            shutil.copyfile(options.datafile, dest)

        if options.annotatorfile:
            #tmp_ann_filename = os.path.basename(options.annotatorfile)
            tmp_ann_filename = "%s.%s" % (record_name, options.inAnnType)
            dest = "%s/%s" % (tmp_wfdb_dir, tmp_ann_filename)
            shutil.copyfile(options.annotatorfile, dest)



        # set the temporary output filename for each program
        if options.program == 'sqrs':
            record_name = re.sub('.hea', '', tmp_record_filename )
            qrs_tmp_outfile = "%s/%s.qrs" % (tmp_wfdb_dir, record_name)
            #print qrs_tmp_outfile
        if options.program == 'wqrs':
            record_name = re.sub('.hea', '', tmp_record_filename )
            wqrs_tmp_outfile = "%s/%s.wqrs" % (tmp_wfdb_dir, record_name)
            #print qrs_tmp_outfile
        if options.program == 'nguess':
            record_name = re.sub('.hea', '', tmp_record_filename )
            if options.outAnnType:
                outType = options.outAnnType
            else:
                outType = "nguess"
            nguess_tmp_outfile = "%s/%s.%s" % (tmp_wfdb_dir, record_name, outType)
            #print nguess_tmp_outfile


        # build command line
        wfdb_dir = "/mnt/galaxyTools/tools/wfdb/10.5.20/bin"
        cmd = "%s/%s" % (wfdb_dir, options.program)
        
        if options.recordfile:
            if tmp_record_filename:
                record_name = re.sub('.hea', '', tmp_record_filename )
                #tmp_record = "%s/%s" % (tmp_wfdb_dir, record_name)
                tmp_record = "%s" % (record_name)
                cmd += " -r %s" % (tmp_record)
#        if options.datafile:
#            if tmp_data_filename:
#                if options.program != "sqrs" and options.program != "wqrs" and options.program != "nguess":
#                    cmd += " -something"
#        if options.qrsfile:
#            if tmp_qrs_filename:
#                if options.program != "sqrs" and options.program != "wqrs" and options.program != "nguess":
#                    cmd += " -something"
#        if options.atrfile:
#            if tmp_atr_filename:
#                if options.program != "sqrs" and options.program != "wqrs" and options.program != "nguess":
#                    cmd += " -something"
        if options.annotatorfile:
            cmd += " -a %s" % (options.inAnnType)
        if options.outAnnType:
            cmd += " -o %s" % (options.outAnnType)
        if options.begintime:
            cmd += " -f %s" % (options.begintime)
        if options.stoptime:
            cmd += " -t %s" % (options.stoptime)
        if options.resolution == "true":
            cmd += " -H"
        if options.jpoints == "true":
            cmd += " -j"
        if options.resample == "true":
            cmd += " -R"
        if options.frequency:
            cmd += " -p %s" % (options.frequency)
        if options.signal:
            cmd += " -s %s" % (options.signal)
        if options.threshold:
            cmd += " -m %s" % (options.threshold)
        if options.insertQann:
            cmd += " -m %s" % (options.insertQann)
        if options.allIntervals:
            cmd += " -A"
        if options.consecutiveIntervals:
            cmd += " -c"
        if options.iformat:
            cmd += " -d %s" % (options.iformat)
        if options.itypeEnd:
            cmd += " -p %s" % (options.itypeEnd)
        if options.itypeBegin:
            cmd += " -P %s" % (options.itypeBegin)
        if options.timesFormatEnd:
            cmd += " -v %s" % (options.timesFormatEnd)
        if options.timesFormatBegin:
            cmd += " -V %s" % (options.timesFormatBegin)
        if options.annotationsEnd:
            cmd += " -w"
        if options.annotationsBegin:
            cmd += " -W"
        if options.intervals:
            cmd += " -l %s" % (options.intervals)
        if options.selectTimePrecision:
            precision_param = ""
            if options.selectTimePrecision == "milliseconds":
                precision_param += "-p"
            elif options.selectTimePrecision == "highestPrecision":
                precision_param += "-P"
            if options.selectTimePrecision != "default":
                if options.selectTimeFormat:
                    if options.selectTimeFormat == "day":
                        precision_param += "d"
                    elif options.selectTimeFormat== "elapsed":
                        precision_param += "e"
                    elif options.selectTimeFormat== "hours":
                        precision_param += "h"
                    elif options.selectTimeFormat== "minutes":
                        precision_param += "m"
                    elif options.selectTimeFormat== "seconds":
                        precision_param += "s"
                    elif options.selectTimeFormat== "interval":
                        precision_param += "S"
            cmd += " %s" % (precision_param)
        if options.signalList:
            cmd += " -s %s" % (options.signalList)
        if options.signalStart:
            cmd += " -S %s" % (options.signalStart)
        if options.printHeader:
            cmd += " -v"
        if options.outputFormat:
            if options.outputFormat == 'csv':
                cmd += " -c"
            elif options.outputFormat == 'xml':
                cmd += " -X"


        if options.program == 'ann2rr' or options.program == 'rdsamp':
            cmd += " > %s" % (options.outputfile)

        print cmd
        #time.sleep(120)
        ## execute job in tempdir, move and rename output to appropriate dir, delete tempdir
        ## execute command in tempdir
        curr_wdir = os.getcwd()
        os.chdir(tmp_wfdb_dir)
        tmp_out = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp_out, 'w' )
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open( tmp_err, 'w' )
        proc = subprocess.Popen( args=cmd, shell=True, cwd="./", stdout=tmp_stdout, stderr=tmp_stderr )

        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp_err, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

        # if job is sqrs, copy the new qrs to the output file
        if options.program == 'sqrs':
            shutil.copyfile(qrs_tmp_outfile, options.outputfile)
        # if job is wqrs, copy the new wqrs to the output file
        if options.program == 'wqrs':
            shutil.copyfile(wqrs_tmp_outfile, options.outputfile)
        # if job is nguess, copy the new nguess to the output file
        if options.program == 'nguess':
            shutil.copyfile(nguess_tmp_outfile, options.outputfile)


    ## delete tempdir
    finally:
        try:
            shutil.rmtree(tmp_wfdb_dir) # delete directory
#            try:
#                with open(stderr): pass
#                os.unlink(stderr)
#            except IOError:
#                print 'Cleanup done: stderr'
#            try:
#                with open(stdout): pass
#                os.unlink(stdout)
#            except IOError:
#                print 'Cleanup done: stdout'
            #print tmp_wfdb_dir
        except OSError, e:
            if e.errno != 2: # code 2 - no such file or directory
                raise
    os.chdir(curr_wdir)
    print "End time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

if __name__=="__main__": __main__()
