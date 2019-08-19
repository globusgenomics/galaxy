#!/usr/bin/python
import tempfile, sys, os, glob, datetime, time, shutil, optparse
from bioblend.galaxy import GalaxyInstance
import requests
requests.packages.urllib3.disable_warnings()
sql_path = '/opt/galaxy/.venv/lib/python2.7/site-packages'
sys.path.append(sql_path)
import sqlalchemy as sa
from sqlalchemy import Table

def get_filepath(UUID):
    GALAXY_DATABASE_CONN = "postgresql://galaxy:globus_genomics_pass@rds.ops.globusgenomics.org:5432/galaxy_pilot1"
    galaxy_db_conn = sa.create_engine(GALAXY_DATABASE_CONN).connect()
    galaxy_db_meta = sa.MetaData(bind=galaxy_db_conn)

    dataset_table = sa.Table("dataset", galaxy_db_meta, autoload=True)

    dataset_id = galaxy_db_conn.execute(sa.sql.select([dataset_table.c.id]).where(dataset_table.c.uuid == UUID)).fetchone()[0]
    dir_num =  get_dir_num(dataset_id)
    file_path = "/scratch/galaxy/files/{0}/dataset_{1}.dat".format(get_dir_num(dataset_id), dataset_id)
    #print file_path
    return file_path

def get_dir_num(dataset_id):
  if dataset_id < 1000:
    return "000"
  tmp = str(dataset_id)[:-3]
  if len(tmp) == 1:
    return "00{0}".format(tmp)
  else:
    return "0{0}".format(tmp)


def file_as_bytes(file):
    with file:
        return file.read()

def is_complete (historyID, gi):
    status = gi.histories.get_status(historyID)
    if status['percent_complete'] == "100":
        return True
    elif status['state'] == 'ok':
        return True
    else:
        return False

#####  When the pipeline completes, the Reports,
##### VCFs, consensus, Clims and access files are sent by email to the user and copied to their
##### final repository for storage.
##### usage : python Parallel_TB_runner.py final_progress.txt pascal.lapierre__at__health.ny.gov /path/to/files


##### Version 3.0 August 15th, 2018
##### Third major update of the pipeline
##### Improve script annotation and mechanics my major help from Alex R.


##### on September 7th, edited out line $pm -> set_waitpid_blocking_sleep(0); in an attempt to solve the missing sample issue

## Script version
version = "2.0_beta"

## Variable to set start time of the runclock  runtime calculation
t0 = time.time()

## Variable to set end time of the runclock for runtime calculation
t1 = None

## time difference between $t1 - $t0
td = None

## variable to set the full date of the analyses
## get localtime
date = datetime.datetime.now().strftime('%Y_%m_%d')

## parse input vars
parser = optparse.OptionParser()
parser.add_option( '-k', '--api-key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
parser.add_option( '-u', '--url', dest='url_name', help='Galaxy URL' )
parser.add_option( '--input', dest='table_file', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--out', dest='log', help='Log file' )
parser.add_option( '--out-dir', dest='output_dir', help='Output directory to store BAM, FASTQ, VCF and more' )
parser.add_option( '--email', dest='email', help='email to send alert to' )
parser.add_option( '--access-log', dest='access_log', help='path to access log' )
parser.add_option( '--clims-log', dest='clims_log', help='path to clims log' )
(options, args) = parser.parse_args()

## variable defining the directory where the running scripts are located
outdir = options.output_dir # File(s) location for the analysis
if not os.path.exists(outdir):
    os.makedirs(outdir)
    os.makedirs("%s/VCF" % outdir)
    os.makedirs("%s/Consensus" % outdir)
    os.makedirs("%s/CLIMS" % outdir)
    os.makedirs("%s/Reports" % outdir)

## Email address(s) where to send the reports. Can contain multiple addresses
## modify the provided email address into right format for sendmail command
email = None
if options.email is not None:
    email = options.email.replace("__at__", "@").replace(" ", ",")  # Email address provided by the user

## logfile placeholder
log_file = options.log # contain the runlog name and location as provided by Galaxy

## batch file log
batch_log_file = options.table_file

## api key
api_key = options.api_key

#url
url = options.url_name

## ARGV[0] contain the runlog name and location as provided by Galaxy
fh_log = open(log_file, "w")
fh_log.write("Script ./Parallel_TB_runner_v2.0_beta.pl Version %s\n\nCurrent working directory %s\nTemporary running directory %s\n" % (version, outdir, outdir))

############## reporting ##########################
## creation of a CLIMS importable report
## and of an Access compatible table

fh_clims = open("%s/CLIMS_Report_%s.out" % (outdir,date), "w")
fh_clims.write("sample_name\tAnalysis date\tpercent reads mapped\tpercent coverage\tavg. depth\tpercent pos. passing threshold\tnumber of snps found\tpipeline_status\tin_silico_genotyping_results\tbioproject_accession\ttitle\tlibrary_id\tdesign_description\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tfiletype\tfilename1\tfilename2\tfilename3\tfilename4\ttoxins_profiling\tresistance_profiling\n")
fh_clims.close()

## open list of all the temp directory that were created
## this file contain the sample names and directory locations
if "http:" in url:
    url = url.replace("http", "https")
if "dev1" in url:
    url = url.replace("dev1", "dev")

# monitor the histories and loop until jobs are complete
gi = GalaxyInstance(url=url, key=api_key)
fh_batch = open(batch_log_file, "r")
monitor_meta = {}
for line in fh_batch:
    line = line.rstrip("\n")
    if line.startswith("SUBMITTED"):
        status,sampleName,wfName,wfID,historyName,historyID = line.split("\t")
        monitor_meta[historyID] = { 'sampleName':sampleName, 'wfName':wfName, 'wfID':wfID, 'historyName':historyName }

completed_meta = {}
outputs = []
while len(monitor_meta) != len(completed_meta):
    for hKey in monitor_meta:
        #print "Analyzing HKEY %s" % hKey
        if hKey in completed_meta.keys():
            continue
        elif is_complete(hKey, gi) == True:
            #fh_out.write("%s\n" % hKey)
            completed_meta[hKey] = {'status' : 'complete', 'sampleName': monitor_meta[hKey]['sampleName'], 'datasets' : {}}

            # get metadata about the completed jobs
            for dataset in gi.histories.show_history(hKey, contents=True):
                meta = gi.histories.show_dataset(hKey, dataset['id'])
                #print meta
                if meta['visible'] is True and meta['deleted'] is False and meta['name'] == "progress log output":
                    job_name = meta['name'].replace(" ", "_")
                    job_name = job_name.replace("/", "_")
                    job_name = job_name.replace(",", "_")
                    UUID = meta['uuid'].replace("-", "")
                    #print UUID
                    file_path = get_filepath(UUID)
                    extra_files_path = "%s/%s_files" % (os.path.dirname(file_path), ('.').join(os.path.basename(file_path).split('.')[:-1]))
                    completed_meta[hKey]['datasets'][dataset['id']] = { 'data_type' : meta['data_type'], 'file_name': file_path, 'job_name' : job_name, 'extra_files_path' : extra_files_path}
                    outputs.append([monitor_meta[hKey]['sampleName'].split("_")[0], file_path, extra_files_path])
    if len(monitor_meta) != len(completed_meta):
        time.sleep(300)


##  Block that copy the reports, consensus, vcf to their final destination $dest
## and add the individual CLIMS line to the final CLIMS report

clims_output = "%s/CLIMS_Report_%s.out" % (outdir, date)
reports_locations = []
for sample_name, file_path, extra_file_path in outputs:
    vcf_file = "%s/VCF/%s.vcf.gz" % (extra_file_path, sample_name)
    ## if VCF file present, mean pipeline completed analysis without fails and these files
    ## were created
    if os.path.exists(vcf_file):
        shutil.copyfile(vcf_file, "%s/VCF/%s" % (outdir,os.path.basename(vcf_file)))
        shutil.copyfile("%s/VCF/%s_INDELS.vcf.gz" % (extra_file_path, sample_name), "%s/VCF/%s_INDELS.vcf.gz" % (outdir, sample_name))
        shutil.copyfile("%s/Consensus/%s.fa" % (extra_file_path, sample_name), "%s/Consensus/%s.fa" % (outdir,sample_name))

    ## Append each individual clims report to main CLIMS file
    clims_report_file = glob.glob("%s/CLIMS/CLIMS_report_*.out" % extra_file_path)[0]
    os.system("cat %s >> %s/CLIMS_Report_%s.out" % (clims_report_file, outdir, date));
    shutil.copy("%s/CLIMS_Report_%s.out" % (outdir, date), "%s/CLIMS/" % outdir)

    ## Store report location for each samples into array for later access and email
    reports_file = glob.glob("%s/Reports/%s_*.txt" % (extra_file_path, sample_name))[0]
    reports_locations.append(reports_file)
    shutil.copy(reports_file, "%s/Reports" % (outdir))

join_reports_location = " ".join(reports_locations)
shutil.copy(reports_file, options.clims_log)

### create access report ##################################################
### this blog scrubs the reports generated by the pipeline and keep all the
### identified mutations in the loci list of interest

output_access = "%s/access_reports.txt" % outdir
fh_outacc = open(output_access, "w")
fh_outacc.write("IDR number\tNucleotide position\tCodon #\tNucleotide change\tAmino acid change\tPv gene notation\tGene name\tDrug\thetero/homozygous\n")

for reports_tmp_directory in reports_locations:
    ## only process samples which has IDR or AFB in their sample name
    ## since reports generated with earlier version are not compatible with this script
    print reports_tmp_directory
    if 'IDR' in reports_tmp_directory or 'AFB' in reports_tmp_directory:
        base = os.path.basename(reports_tmp_directory)
        sample_name= base.split('.')[0].split("_")[0].split("-")[0]
        print sample_name

        fh_report = open(reports_tmp_directory, "r")
        for inline in fh_report:
            if 'All Mutations in screened loci' in inline:
                inline = fh_report.next()
                while '-------' not in inline:
                    if '->' in inline and 'locus missing' not in inline:
                        fh_outacc.write("%s\t%s" % (sample_name, inline))
                    inline = fh_report.next()
                break
        fh_report.close()
fh_outacc.close()
shutil.copyfile(output_access, options.access_log)

################################################################

#system ("sendemail -t '$email' -m \"Please do not reply to this email.  For any questions or issues, contact the pipeline administrator at pascal.lapierre\@health.ny.gov\" -u \"MTB pipeline reports $date\" -a $rundir/CLIMS_Report\_$date\.out $rundir/access_reports.txt $join_reports_location -f bioinfo\@health.ny.gov");

### Cleanup temp run directory
#system ("rm -r $rundir");

t1 = time.time()
td = t1 - t0
fh_log.write("\n\nPipeline runtime: %s seconds \n" % td)
fh_log.close()



