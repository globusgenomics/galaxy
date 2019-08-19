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
parser.add_argument('-i',dest='input_dir',help=('full path of sample fastq directory'),required=True)
parser.add_argument('-t',dest='target_file',help=('Target list file'), required=True)
parser.add_argument('-b',dest='target_bed',help=('Target bed file'), required=True)
#parser.add_argument('-t',dest='tmp_dir', default='/scratch')

## optional trimming
#parser.add_argument('--trimming', dest='trimming', action='store_true')
#parser.add_argument('--no-trimming', dest='trimming', action='store_false')
#parser.set_defaults(trimming=False)

## germline or somatic boolean option
#parser.add_argument('--germline',dest='germline',action='store_true')
#parser.add_argument('--somatic',dest='somatic',required='store_false')
#parser.set_defaults(germline=True)

parser.add_argument('--reference', dest='reference', help=('Fasta reference file'), required=True)
parser.add_argument('--bwa-reference', dest='bwa_reference', help=('BWA reference file'), required=True)

args = parser.parse_args()
#runName = args.runName
sampleId = args.sampleId
fastq_dir = args.input_dir
#work_dir = args.work_dir
output_dir = args.work_dir
reference = args.reference
bwa_reference = args.bwa_reference
target_file = args.target_file
target_bed = args.target_bed
#trimming = args.trimming


# make temp directory
#work_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-dnaExome")
#tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimized-tmp-")
# make temp directory
work_dir = tempfile.mkdtemp(prefix="optimized-")
tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

#### Slack authentication
#SLACK_TOKEN="xoxp-166523070336-166523070800-258980986992-dbbf6a0d6d8b92d316ea0b37543a7a2c"
#slack_client = SlackClient(SLACK_TOKEN)

def send_message(channel_id, message):
   slack_client.api_call(
         "chat.postMessage",
         channel=channel_id,
         text=message,
         username='pipeline',
         icon_emoji=':bulb:',
   )

#fastq_dir  = os.path.join(work_dir,runName,sampleId)


## TODO change output_dir to NGS_data ...
#output_dir = os.path.join(work_dir,runName,sampleId,'pipeline_output')

print 'the FASTQ directory is at: ',fastq_dir
## reference
#reference_test = r'/Volumes/lab data/reference/hg19/hg19.fasta'
#reference = r'/home/lee/reference/hg19/hg19.fasta'
#adapter = 'ILLUMINACLIP:/home/lee/NGS/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10'

# this is the QXT
#adapter = 'ILLUMINACLIP:/home/lee/NGS/Trimmomatic-0.36/adapters/QXTIllumina-PE.fa:2:30:10'

'''
Define reference and software BWA, samtools, PICARD, GATK, Trimmomatic
'''


PICARD_JAR = '/mnt/galaxyTools/tools/picard/2.7.1/picard.jar'
picard              = 'java -jar '+ PICARD_JAR
picard_more_memory  = 'java -jar -Xmx20g ' + PICARD_JAR

#TRIMMOMATIC_JAR = '/home/lee/NGS/Trimmomatic-0.36/trimmomatic-0.36.jar'
#trimmomatic     = 'java -jar '+TRIMMOMATIC_JAR

#SURECALLTRIMMER_JAR = '/home/lee/NGS/AGeNT/SurecallTrimmer_v4.0.1.jar'
#surecalltrimmer     = 'java -jar '+SURECALLTRIMMER_JAR

GATK_JAR = '/mnt/galaxyTools/tools/gatk4/gatk-4.0.1.2/gatk-package-4.0.1.2-local.jar'
gatk     = 'java -jar -Xmx2g '+GATK_JAR
gatk_more_memory = 'java -jar -Xmx20g '+GATK_JAR

SNPEFF_JAR = '/mnt/galaxyTools/tools/snpeff/snpEff_4.1/snpEff.jar'
SNPEFF_CONFIG = '/mnt/galaxyTools/tools/snpeff/snpEff_4.1/snpEff.config'
snpeff = 'java -jar -Xmx10g '+SNPEFF_JAR


## FASTQ - checking FASTQ file format & Perform FASTQC
if not os.path.exists(fastq_dir):
    log.error("FASTQ_Dir does not exist: %s" % fastq_dir)
    sys.exit(1)

read1_lanes = []
read2_lanes = []
#r1_pattern = re.compile(r'.*L\d{3}_R1_\d{3}_001.fastq\.gz')
r1_pattern = re.compile(r'.*_R1_001.fastq.gz')
#r2_pattern = re.compile(r'.*L\d{3}_R2_\d{3}_001.fastq\.gz')
r2_pattern = re.compile(r'.*_R2_001.fastq.gz')


for fastq in glob.glob(os.path.join(fastq_dir,'*.fastq.gz')):
    print fastq
    #print FASTQC
    if r1_pattern.match(fastq):
        read1_lanes.append(fastq)
    elif r2_pattern.match(fastq):
        read2_lanes.append(fastq)
    else:
        log.warn("FASTQ File %s did not match pattern", fastq)

print("Found %d lanes for sample %s", len(read1_lanes), sampleId)
print 'FASTQ files: ',read1_lanes, read2_lanes



## create sam file list
sam_file_list = []

try:
    os.makedirs(output_dir)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise


## Run FASTQC
#fastqc_cmd = 'fastqc ' + read1_lanes + ' '+ read2_lanes +' -o ' + output_dir
fastqc_cmd = 'fastqc ' + ''.join(map(str, read1_lanes)) + ' '+ ''.join(map(str, read2_lanes)) +' -o ' + output_dir

print "====== Running FASTQC ======="
subprocess.call(fastqc_cmd, shell=True,
                stdout=open(os.path.join(output_dir, sampleId + '_fastqc_pipeline_log.txt'), 'wt'),
                stderr=subprocess.STDOUT)


## Align FASTQ
for index,file in enumerate(read1_lanes):
    # Alignment inputs
    read1 = read1_lanes[index]
    read2 = read2_lanes[index]

    # Alignment outputs
    lane = index + 1
    #sam        = os.path.join(output_dir,'lane' + str(lane) + '.sam')

    #samplename = os.path.basename(read1).split("_R1_001")[0]
    sam        = os.path.join(output_dir,sampleId + '.sam')

    print sampleId


    ## BWA align
    ## Step 1. Run bwa mem, FASTQ --> SAM
    ## Step 2. samtools to convert SAM --> BAM

    RG_Tag = '"@RG\tID:' + sampleId + '\tPU:' + str(lane) + '\tSM:' + sampleId + '"'
    #RG_Tag = '"@RG\tID:' + runName + '\tPU:' + '\tSM:' + samplename + '"'
    bwa_cmd = 'bwa mem' + ' -t 20' + ' -v 3 '+ ' -R '+ RG_Tag + ' '+'"' +bwa_reference+ '"' + ' ' + '"' +read1 + '"' +' ' + '"'+read2 + '"' + ' > ' + sam
    #samtools_convert = 'samtools view -S -b'
    #samtools_sort    = 'samtools sort -o '+ ' ' +'>'+ ' ' + samplename+'_sorted.bam'
    print 'bwa comand: ',bwa_cmd
    ## use subprocess to execute the bash commands
    #outfile = open(samplename + '_sorted.bam','w')
    print "====== Running BWA Mem into memory ======="
    subprocess.call(bwa_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_bwa_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)
    #bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    #stdout, stderr = bwa_process.communicate()
    sam_file_list.append(sam)
    print 'sam files: ',sam

## Merge SAM Files and Sort

## Create a list of files to merge
## Picard SortSam can convert SAM to BAM and sort
dedup_bam_list = []
for file in sam_file_list:
    #samplename = os.path.basename(file).split('.')[0]
    BAM = os.path.join(output_dir, sampleId + '_sorted.bam')
    BAI = os.path.join(output_dir, sampleId + '_sorted.bai')
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
    sort_cmd = picard_more_memory +' SortSam VALIDATION_STRINGENCY=SILENT Sort_Order=coordinate'+ option2 + file + option3 + BAM + ' TMP_DIR=' + tmp_dir + option4
    print sort_cmd

    print ('======= sorted BAM Running =======')
    subprocess.call(sort_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_sortSam_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## markDuplication
    dedup_cmd = picard_more_memory +' MarkDuplicates REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT'+option4+option2+BAM+option3+DEDUP_BAM+option5+DEDUP_METRICS
    print ('======= markDuplicates Running =======')
    print (dedup_cmd)
    subprocess.call(dedup_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_deDup_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## Generate Alignment summary metrics
    align_summary_cmd = picard_more_memory+' CollectAlignmentSummaryMetrics' + option1 + option2 + DEDUP_BAM + option3 + ALIGNMENT_METRICS
    print ('======= picard alignment_summary_metrics =======')
    subprocess.call(align_summary_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_alignment_summary_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## Generate HSMetrics
    TARGET_Metrics_cmd = picard_more_memory+TARGET_option1+TARGET_option2+TARGET_option3+TARGET_option4+option2+DEDUP_BAM+option3+TARGET_METRICS
    print ('======= picard TargetPcr_metrics =======')
    subprocess.call(TARGET_Metrics_cmd,shell=True,stdout=open(os.path.join(output_dir,sampleId+'_TargetPcrMetrics_pipeline_log.txt'),'wt'),stderr=subprocess.STDOUT)

    ## remove sorted_BAM file and _sorted_BAI index
    #os.remove(os.path.abspath(BAM))
    #os.remove(os.path.abspath(BAI))

    ## remove PICARD tmp directory
    #shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)
    dedup_bam_list.append(DEDUP_BAM)
    print 'dedup BAM: ',DEDUP_BAM

VCF = os.path.join(output_dir, sampleId + '.vcf')
## RUN GATK on dedup BAM against target intervals from Agilent SureSelect QXT Targets
for bam in glob.glob(os.path.join(output_dir,'*_sorted_dedup.bam')):
    print "Input BAM file is: ",bam
    #samplename = os.path.basename(bam).split('_sorted.')[0]

    GVCF_SPARK = os.path.join(output_dir, sampleId + '_spark.vcf')
    option1 = ' -R %s' % reference
    #option1_spark = ' -R /home/lee/reference/hg19/hg19.2bit'
    option2 = ' -I '
    option3 = ' -O '
    #option4 = ' --emit-ref-confidence GVCF '
    option5 = ' -L %s' % target_file
    option6 = ' -stand_call_conf 30'
    option7 = ' -stand_emit_conf 10'

    GATK_cmd =gatk_more_memory + ' HaplotypeCaller' + option1 + option2 + bam + option3 + VCF + option5
    #GATK_SPARK='/home/lee/NGS/GATK/gatk-4.0.1.1/gatk HaplotypeCallerSpark '+option1_spark+option2+bam+option3+GVCF_SPARK+option4+option5

    print ('======= GATK HaplotypeCaller Running =======',GATK_cmd)
    subprocess.call(GATK_cmd, shell=True,cwd=os.path.join(output_dir),
                    stdout=open(os.path.join(output_dir, sampleId + '_GATK_pipeline_log.txt'), 'wt'),
                    stderr=subprocess.STDOUT)

SNP_VCF = os.path.join(output_dir, sampleId + '_snp.vcf')
INDEL_VCF = os.path.join(output_dir, sampleId + '_indel.vcf')
## Run GATK SNP and INDEL VariantSelect
for vcf in glob.glob(os.path.join(output_dir,'*.vcf')):
    print "Input VCF file is: ",vcf

    option1 = ' -R %s' % reference
    option2 = ' -V '
    option_snp   = ' --select-type-to-include SNP'
    option_indel = ' --select-type-to-include INDEL'
    option3 = ' -O '
    GATK_select_SNP = gatk_more_memory + ' SelectVariants' + option1 + option2 + vcf + option_snp  + option3 + SNP_VCF
    GATK_select_INDEL=gatk_more_memory + ' SelectVariants' + option1 + option2 + vcf + option_indel+ option3 + INDEL_VCF
    print ('======= GATK selectSNP Running =======',GATK_select_SNP)
    subprocess.call(GATK_select_SNP, shell=True,
                    stdout=open(os.path.join(output_dir, sampleId + '_GATK_pipeline_log.txt'), 'a+'),
                    stderr=subprocess.STDOUT)
    print ('======= GATK selectINDEL Running =======',GATK_select_INDEL)
    subprocess.call(GATK_select_INDEL, shell=True,
                    stdout=open(os.path.join(output_dir, sampleId + '_GATK_pipeline_log.txt'), 'a+'),
                    stderr=subprocess.STDOUT)


## remove sam file
#for file in sam_file_list:
#    ## remover sam file to save space
#    os.remove(os.path.abspath(file))


## Generating QC metrics report\
## Step 1. Convert BAM to BED, make coverage depth file
## Step 2. Collect all QC metrics and put together to csv file
print ('======= Calculating QC_metrics Report =======')
collect_Bedtools_coverage(output_dir,sampleId,target_bed)
collect_Final_QCMetrics(output_dir,sampleId)

#if os.path.exists(SNP_VCF):
#    ## send slack message
#    message = 'Congratulations! Analysis for Run: '+runName+' Sample: '+sampleId+' is Finished! :beers: :thumbsup: :facepunch:'
#    #send_message('HCP-pipeline', message)
#    print 'Slack Message has been sent'

print ('======= HCP pipeline finished!!!=======')
