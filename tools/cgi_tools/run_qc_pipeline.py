'''
# #===========================================================================##
# # HCP_pipeline.py
# # this is the HCP pipeline script that construct using the Exome pipeline backbone
# #
# # Version 1.0
# # 1. checking FASTQ name
# # 2. FASTQC
# # 3. BWA alignment
# # 4. PICARD sort, markduplicated, TargetPcrMetrics, collectAlignmentMetrics
# # 5. GATK HaplotypeCaller
# # 6. GATK SelectVariant to seperate SNP and INDEL
# #===========================================================================##
# #===========================================================================##
# # All of the general library calls and PATH
# #===========================================================================##
'''

'''
USAGE:
python pipeline.py -w <work_dir> -r <runName> -s <sampleId> -t <tmp_dir> <option:--trimming/--no-trimming> 
'''



#!/usr/bin/python
'''
import packages
'''
import argparse
import errno
import logging as log
import re
import shutil
import sys
import tempfile

#from slackclient import SlackClient

from target_qc_metrics import *

## working directory
#os.chdir(r'/Volumes/lab data/Rutherford_data/HCP_data/161130_NB501330_0059_AHGWVLAFXX/Data/Intensities/BaseCalls')

## add parameter
parser = argparse.ArgumentParser()
parser.add_argument('-w',dest='work_dir',help=('full path of working directory'),required=True)
#parser.add_argument('-r',dest='runName',help=('full path of running directory wit data'),required=True)
parser.add_argument('-s',dest='sampleId',help=('option to execute individually by sampleId'),required=True)
parser.add_argument('-i',dest='input_bam',help=('full path of sample bam file'),required=True)
parser.add_argument('-t',dest='target_file',help=('Target list file'), required=True)
parser.add_argument('-b',dest='target_bed',help=('Target bed file'), required=True)
#parser.add_argument('-t',dest='tmp_dir', default='/scratch')

parser.add_argument('--reference', dest='reference', help=('Fasta reference file'), required=True)

args = parser.parse_args()
sampleId = args.sampleId
bam_file = args.input_bam
#work_dir = args.work_dir
output_dir = args.work_dir
reference = args.reference
target_file = args.target_file
target_bed = args.target_bed


# make temp directory
#work_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-dnaExome")
#tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimized-tmp-")
# make temp directory
work_dir = tempfile.mkdtemp(prefix="optimized-")
tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")


## TODO change output_dir to NGS_data ...
#output_dir = os.path.join(work_dir,runName,sampleId,'pipeline_output')

'''
Define reference and software samtools, PICARD
'''

PICARD_JAR = '/mnt/galaxyTools/tools/picard/2.7.1/picard.jar'
picard              = 'java -jar '+ PICARD_JAR
picard_more_memory  = 'java -jar -Xmx20g ' + PICARD_JAR

## create sam file list
sam_file_list = []

try:
    os.makedirs(output_dir)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise


## Create a list of files to merge
## Picard SortSam can convert SAM to BAM and sort
dedup_bam_list = []
input_file_list = [bam_file]
for file in input_file_list:
    #samplename = os.path.basename(file).split('.')[0]
    #BAM = os.path.join(output_dir, sampleId + '_sorted.bam')
    BAM = file
    DEDUP_BAM = os.path.join(output_dir, sampleId + '_sorted_dedup.bam')
    DEDUP_METRICS = os.path.join(output_dir,sampleId+'_sorted_dedup_metrics.txt')
    TARGET_METRICS    = os.path.join(output_dir,sampleId+'_TargetPcrMetrics.txt')

    ALIGNMENT_METRICS = os.path.join(output_dir,sampleId+'_aligment_metrics.txt')
    option1 = ' R=%s' % reference
    option2 = ' I='
    option3 = ' O='
    option4 = ' CREATE_INDEX=TRUE'
    option5 = ' M='

    TARGET_option1 = ' CollectTargetedPcrMetrics VALIDATION_STRINGENCY=SILENT'
    TARGET_option2 = ' REFERENCE_SEQUENCE=%s' % reference
    TARGET_option3 = ' AMPLICON_INTERVALS=%s' % target_file
    TARGET_option4 = ' TARGET_INTERVALS=%s' % target_file

    #merge_input_string = merge_input_string + ' I=' + file
    tmp_dir = os.path.join(output_dir,'tmp')

    ## markDuplication
    dedup_cmd = picard_more_memory +' MarkDuplicates REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT'+option4+option2+BAM+option3+DEDUP_BAM+option5+DEDUP_METRICS
    print ('======= markDuplicates Running =======')
    print (dedup_cmd)
    subprocess.call(dedup_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_deDup_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## Generate Alignment summary metrics
    align_summary_cmd = picard_more_memory+' CollectAlignmentSummaryMetrics' + option1 + option2 + DEDUP_BAM + option3 + ALIGNMENT_METRICS
    print ('======= picard alignment_summary_metrics =======')
    print (align_summary_cmd)
    subprocess.call(align_summary_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_alignment_summary_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## Generate HSMetrics
    TARGET_Metrics_cmd = picard_more_memory+TARGET_option1+TARGET_option2+TARGET_option3+TARGET_option4+option2+DEDUP_BAM+option3+TARGET_METRICS
    print ('======= picard TargetPcr_metrics =======')
    print (TARGET_Metrics_cmd)
    subprocess.call(TARGET_Metrics_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_TargetPcrMetrics_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## remove PICARD tmp directory
    #shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)
    dedup_bam_list.append(DEDUP_BAM)
    print 'dedup BAM: ',DEDUP_BAM

## Generating QC metrics report\
## Step 1. Convert BAM to BED, make coverage depth file
## Step 2. Collect all QC metrics and put together to csv file
print ('======= Calculating QC_metrics Report =======')
collect_Bedtools_coverage(output_dir,sampleId,target_bed)
collect_Final_QCMetrics(output_dir,sampleId)

print ('======= HCP pipeline finished!!!=======')
