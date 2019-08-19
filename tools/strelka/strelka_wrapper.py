#!/usr/bin/env python
#Gregoire Seguin-Henry (Engineer IT)
#Amine Sbitti (Data Scientist)
#Ludovic Marie-Sainte (Project Manager)
#For Geviteam 2014

"""
A wrapper script for running the GenomeAnalysisTK.jar commands.
"""
from __future__ import print_function
import sys, argparse, os, tempfile, subprocess, shutil
from binascii import unhexlify
from string import Template
from galaxy import eggs

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def _create_config(args, config_path):
    conf_file = open(config_path, "w")
    conf_file.write("[user]\n")
    args2 = vars(args)
    for option in args2:
        if not option in ["tumorBam", "normalBam", "refFile", "configFile", "scriptPath", "a", "b", "c", "d", "e", "extraStrelkaArguments"] and args2[option]!=None:
	    conf_file.write("%s=%s\n" % (option, args2[option]))
    if args.extraStrelkaArguments == "yes":
    	conf_file.write("extraStrelkaArguments=")
   	if args.a:
            conf_file.write("--ignore-conflicting-read-names ")
    	if args.b != None:
            conf_file.write("-used-allele-count-min-qscore %s " % (args.b))
    	if args.c != None:
            conf_file.write("--candidate-indel-input-vcf %s " % (args.c))
    	if args.d != None:
            conf_file.write("--force-output-vcf %s " % (args.d))
    	if args.e != None:
            conf_file.write("-min-small-candidate-indel-read-frac %s " % (args.e))
    	conf_file.write("\n")
    conf_file.close()

def my_Popen(cmd, prefix_for_stderr_name, tmp_dir, msg_error):
    stderr_name = tempfile.NamedTemporaryFile( prefix = prefix_for_stderr_name ).name
    proc = subprocess.Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
    return_code = proc.wait()                          
    if return_code:
	for line in open( stderr_name ):
           print(line, file=sys.stderr)
	os.unlink( stderr_name ) #clean up
 	cleanup_before_exit( tmp_dir )
 	raise Exception( msg_error )
    else:
        os.unlink( stderr_name )

def index_bam_files( bam_filenames, tmp_dir ):
    for bam_filename in bam_filenames:
        bam_index_filename = "%s.bai" % bam_filename
        print("bam_filename is: " + bam_filename + " bam_index_filename is: " + bam_index_filename + " test is: %s" % os.path.exists(bam_index_filename))
        if not os.path.exists( bam_index_filename ):
            #need to index this bam file
            command = 'samtools index %s %s' % ( bam_filename, bam_index_filename )
            my_Popen( command, "bam_index_stderr", tmp_dir, "Error during indexation of fasta file :" + bam_filename)

def index_fasta_files( fasta_filenames, tmp_dir ):
    for fasta_filename in fasta_filenames:
        fasta_index_filename = "%s.fai" % fasta_filename
        print("fasta_filename is: " + fasta_filename + " fasta_index_filename is: " + fasta_index_filename + " test is: %s" % os.path.exists(fasta_index_filename))
        if not os.path.exists( fasta_index_filename ):
            #need to index this bam file
            command = 'samtools faidx %s %s' % ( fasta_filename, fasta_index_filename )
            my_Popen( command, "fasta_index_stderr", tmp_dir, "Error during indexation of fasta file :" + fasta_filename)

def __main__():
    #Manage options
    print(os.environ['PATH'])
    parser = argparse.ArgumentParser()                                             
    parser.add_argument( '--tumorBam', help='path to tumor bam file', required = False )
    parser.add_argument( '--normalBam', help='', required = False )   
    parser.add_argument( '--refFile', help='', required = False )
    parser.add_argument( '--configFile', help='', required = False )
    parser.add_argument( '--depthFilterMultiple', help='', required = False )
    parser.add_argument( '--snvMaxFilteredBasecallFrac', help='', required = False )
    parser.add_argument( '--snvMaxSpanningDeletionFrac', help='', required = False )
    parser.add_argument( '--indelMaxRefRepeat', help='', required = False )
    parser.add_argument( '--indelMaxWindowFilteredBasecallFrac', help='', required = False )
    parser.add_argument( '--indelMaxIntHpolLength', help='', required = False )
    parser.add_argument( '--ssnvPrior', help='', required = False )
    parser.add_argument( '--sindelPrior', help='', required = False )
    parser.add_argument( '--ssnvNoise', help='', required = False )
    parser.add_argument( '--sindelNoise', help='', required = False )
    parser.add_argument( '--ssnvNoiseStrandBiasFrac', help='', required = False )
    parser.add_argument( '--minTier1Mapq', help='', required = False )
    parser.add_argument( '--minTier2Mapq', help='', required = False )
    parser.add_argument( '--ssnvQuality_LowerBound', help='', required = False )
    parser.add_argument( '--sindelQuality_LowerBound', help='', required = False )
    parser.add_argument( '--isWriteRealignedBam', help='', required = False )
    parser.add_argument( '--binSize', help='path to tumor bam file', required = False )
    parser.add_argument( '--extraStrelkaArguments', help='', required = False )
    parser.add_argument( '--isSkipDepthFilters', help='', required = False )
    parser.add_argument( '--maxInputDepth', help='', required = False )
    parser.add_argument( '--scriptPath', help='', required = False )
    parser.add_argument( '-a', action="store_true", help='', required = False )
    parser.add_argument( '-b', help='', required = False )
    parser.add_argument( '-c', help='', required = False )
    parser.add_argument( '-d', help='', required = False )
    parser.add_argument( '-e', help='', required = False )
    args = parser.parse_args()

    root_dir= args.scriptPath
    expected_dir="for_tests"
    job_dir=os.getcwd()
    analysis_dir=job_dir + "/StrelkaAnalysis"
    config_script=root_dir + "/configureStrelkaWorkflow.pl"
    tmp_dir = tempfile.mkdtemp( prefix='tmp-strelkaAnalysis-' )
    config_ini = "%s/config.ini" % (tmp_dir)

    print("root_dir: " + root_dir + "\njob_dir :" + job_dir + "\nanalysis_dir :" + analysis_dir + "\nconfig_script :" + config_script + "\ntmp_dir :" + tmp_dir + "\nconfig_ini :" +  config_ini)


    #verifying eveything's ok
    if not os.path.isfile(config_script):
    	sys.exit("ERROR: The strelka workflow must be built prior to running. See installation instructions in '$root_dir/README'")
    print("configuring...", file=sys.stdout)
    if os.path.exists(analysis_dir):
	sys.exit("'" + analysis_dir + "' already exist, if you are executing this tool from galaxy it should not happen")
    

    # creating index if needed
    bam_filenames = [ args.tumorBam, args.normalBam ]
    index_bam_files( bam_filenames, tmp_dir )
    fasta_files = [ args.refFile ]
    index_fasta_files( fasta_files, tmp_dir )
    
    #creating config file if needed
    if args.configFile == "Custom":
    	_create_config(args, config_ini)
    elif args.configFile in ["strelka_config_bwa_default.ini", "strelka_config_isaac_default.ini", "strelka_config_eland_default.ini"]:
        cmdbash="cp %s %s" % (root_dir + "/lib/" + args.configFile, config_ini)
        my_Popen(cmdbash, "copy_default_file_err", tmp_dir, "Error during the copy of default config file, maybe it was removed")
    else:
    	if not os.path.exists(args.configFile):
	     print( "The path to your configuration File seems to be wrong, use another one or custom option", file=sys.stderr)
    	cmdbash="cp %s %s" % (args.configFile, config_ini)
        my_Popen(cmdbash, "copy_default_file_err", tmp_dir, "Error during the copy of the selected config file")




    #configuration of workflow
    cmd="%s --tumor=%s --normal=%s --ref=%s --config=%s --output-dir=%s" % (config_script, args.tumorBam, args.normalBam, args.refFile, config_ini, analysis_dir)
    print( "**** Starting configuration.")
    print( "**** Configuration cmd: '" + cmd + "'")
    my_Popen( cmd, "cinfugation_stderr", tmp_dir, "Error during configuration !")
    print("completed configuration")
    
    #run the workflow !
    cmd="make -C " + analysis_dir
    print("**** starting workflow.")
    print("**** workflow cmd: '" + cmd + "'")
    my_Popen( cmd, "workflow_stderr", tmp_dir, "Error during workflow execution !")   
    print("**** completed workflow execution")
    
    cmdbash="cp %s %s" % (config_ini, analysis_dir + "/config.ini")
    my_Popen(cmdbash, "copy_final_conf_file_err", tmp_dir, "Error during the copy of conf file after job is done, quite strange...")  


if __name__=='__main__':
    __main__()
