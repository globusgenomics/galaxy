#!/usr/bin/python
import sys
import tarfile
import requests
import os
from bdbag import bdbag_api
import urllib
import time
import json
from shutil import copyfile, rmtree
import subprocess
from subprocess import STDOUT,PIPE
import boto3
from botocore.client import Config
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
sys.path.append("/mnt/galaxyTools/tools/pymodules/python2.7/lib/python/globus_sdk-1.5.0-py2.7.egg")
import globus_sdk
from datetime import datetime
import hashlib
import uuid 
from bdbag import bdbag_ro as ro
import getpass
import tempfile
from identifier_client.identifier_api import IdentifierClient
from minid_client import minid_client_api
import fnmatch
from pprint import pprint
import ConfigParser

from optparse import OptionParser



TMP_DIR = "/ephemeral/0/tmp/"

def check_globus_token(token, config_file, globus_dir, username, token_type):
   token_type_to_server = {
       'transfer': 'transfer.api.globus.org',
       'auth': 'auth.globus.org',
       'identifier': 'identifiers.globus.org'
   }
   if token_type == 'identifier':
       Config = ConfigParser.ConfigParser()
       Config.read(config_file)
       client_id = Config.get('commons', 'client_id')
       client_secret = Config.get('commons', 'client_secret')
       client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
       token_file = os.path.join(globus_dir, 'tokens-identifiers', username)
   else:
       Config = ConfigParser.ConfigParser()
       Config.read(config_file)
       client_id = Config.get('globus', 'client_id')
       client_secret = Config.get('globus', 'client_secret')
       client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
       if token_type == 'transfer':
           token_file = os.path.join(globus_dir, 'tokens', username)
       elif token_type == 'auth':
           token_file = os.path.join(globus_dir, 'tokens-auth', username)
       else:
           sys.exit('token_type unsupported')
   if not client.oauth2_validate_token(token)['active']:
       refresh_token_file = token_file + '_refresh'
       if os.path.isfile(refresh_token_file):
           with open(refresh_token_file, 'r') as f:
               refresh_token = f.readline().strip()
           token_response = client.oauth2_refresh_token(refresh_token)
           token = token_response.by_resource_server[token_type_to_server[token_type]]['access_token']
           with open(token_file, 'w') as write_file:
               write_file.write(token)
   return token

def create_settings(setting_file, input_path, output_dir, work_dir, tmp_dir, vtype, db_list):
    option = None
    if vtype == "snps":
        option = "s"
    elif vtype == "indels":
        option = "i"
    elif vtype == "both":
        option = "b"

    fh = open(setting_file, "w")
    fh.write("input file name:    %s\n" % input_path)
    fh.write("output file name:   %s/out.vcf.annotated" % output_dir)
    fh.write("resources dir:    %s\n" % "/WGSA/resources/")
    fh.write("annovar dir:    %s\n" % "/WGSA/annovar20160201/annovar/")
    fh.write("snpeff dir:    %s\n" % "/WGSA/snpeff/snpEff/")
    fh.write("vep dir:    %s\n" % "/WGSA/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/")
    fh.write(".vep dir:    %s\n" % "/WGSA/.vep/")
    fh.write("tmp dir:    %s\n" % tmp_dir)
    fh.write("work dir:    %s\n" % work_dir)
    fh.write("retain intermediate file:           n\n")
    for db in db_list:
        fh.write("%s:                    %s\n" % (db, option))

    #fh.write("ANNOVAR/Ensembl:                    %s\n" % option)
    #fh.write("ANNOVAR/RefSeq:                     %s\n" % option)
    #fh.write("ANNOVAR/UCSC:                       %s\n" % option)
    #fh.write("SnpEff/Ensembl:                     %s\n" % option)
    #fh.write("SnpEff/RefSeq:                      %s\n" % option)
    #fh.write("VEP/Ensembl:                        %s\n" % option)
    #fh.write("VEP/RefSeq:                         %s\n" % option)
    #fh.write("gene-model annotation on the fly:   n\n")
    #fh.write("dbSNP:                              %s\n" % option)
    #fh.write("snoRNA miRNA:                       %s\n" % option)
    #fh.write("UTR3 miRNA target:                  %s\n" % option)
    #fh.write("scSNV deleteriousness prediction:   n\n")
    #fh.write("SPIDEX:                             n\n")
    #fh.write("Aloft:                              %s\n" % option)
    #fh.write("GWAS catalog:                       n\n")
    #fh.write("GRASP:                              n\n")
    #fh.write("Clinvar:                            %s\n" % option)
    #fh.write("COSMIC:                             n\n")
    #fh.write("GTEx:                               %s\n" % option)
    #fh.write("Geuvadis:                           %s\n" % option)
    #fh.write("Duke mappability:                   n\n")
    #fh.write("Duke mappability (+-149bp):         n\n")
    #fh.write("Genome mappability score:           n\n")
    #fh.write("1000G mask:                         %s\n" % option)
    #fh.write("RepeatMasker:                       %s\n" % option)
    #fh.write("EPO ancestral:                      n\n")
    #fh.write("AltaiNeandertal genotypes:          n\n")
    #fh.write("Denisova genotypes:                 n\n")
    #fh.write("VindijiaNeandertal genotypes:       n\n")
    #fh.write("PhyloP_primate:                     n\n")
    #fh.write("PhyloP_placental:                   n\n")
    #fh.write("PhyloP_vertebrate:                  n\n")
    #fh.write("PhastCons_primates:                 n\n")
    #fh.write("PhastCons_placental:                n\n")
    #fh.write("PhastCons_vertebrate:               n\n")
    #fh.write("GERP++:                             n\n")
    #fh.write("SiPhy:                              n\n")
    #fh.write("bStatistic:                         n\n")
    #fh.write("fitCons:                            %s\n" % option)
    #fh.write("LINSIGHT:                           %s\n" % option)
    #fh.write("GenoCanyon:                         %s\n" % option)
    #fh.write("1000G phase 3 allele frequencies:   %s\n" % option)
    #fh.write("UK10K allele frequencies:           %s\n" % option)
    #fh.write("ESP6500 allele frequencies:         %s\n" % option)
    #fh.write("ExAC frequencies:                   %s\n" % option)
    #fh.write("ExAC nonTCGA subset frequencies:    %s\n" % option)
    #fh.write("ExAC nonpsych subset frequencies:   %s\n" % option)
    #fh.write("gnomAD exomes frequencies:          %s\n" % option)
    #fh.write("gnomAD genomes frequencies:         %s\n" % option)
    #fh.write("RegulomeDB:                         %s\n" % option)
    #fh.write("funseq like non-coding:             n\n")
    #fh.write("funseq2 non-coding:                 n\n")
    #fh.write("CADD:                               n\n")
    #fh.write("CADDindel:                          n\n")
    #fh.write("DANN:                               %s\n" % option)
    #fh.write("fathmm-MKL:                         %s\n" % option)
    #fh.write("fathmm-XF coding:                   %s\n" % option)
    #fh.write("fathmm-XF non-coding:               %s\n" % option)
    #fh.write("Eigen and EigenPC:                  n\n")
    #fh.write("ORegAnno:                           n\n")
    #fh.write("Topologically_Associating_Domains:  n\n")
    #fh.write("ENCODE_TFBS:                        %s\n" % option)
    #fh.write("ENCODE_Dnase:                       %s\n" % option)
    #fh.write("EnhancerFinder:                     %s\n" % option)
    #fh.write("SuperEnhancer:                      %s\n" % option)
    #fh.write("Genehancer:                         %s\n" % option)
    #fh.write("FANTOM5_enhancer_permissive:        %s\n" % option)
    #fh.write("FANTOM5_enhancer_robust:            %s\n" % option)
    #fh.write("FANTOM5_enhancer_target:            %s\n" % option)
    #fh.write("FANTOM5_enhancer_expression:        %s\n" % option)
    #fh.write("FANTOM5_CAGE_peak_permissive:       %s\n" % option)
    #fh.write("FANTOM5_CAGE_peak_robust:           %s\n" % option)
    #fh.write("EnsemblRB_Overviews:                %s\n" % option)
    #fh.write("EnsemblRB_TFBS:                     %s\n" % option)
    #fh.write("dbNSFP3_variant:                     n\n")
    #fh.write("EnsemblRB_Cell_Type_Activity:       %s\n" % option)
    #fh.write("EnsemblRB_Cell_Type_Segmentations:  %s\n" % option)
    #fh.write("ENCODE_Cell_Type_Segmentations:     %s\n" % option)
    #fh.write("Roadmap-15-state_model:             %s\n" % option)
    #fh.write("Roadmap-25-state_model:             %s\n" % option)
    #fh.write("Roadmap_peak_calls:                 %s\n" % option)
    #fh.write("GenoSkyline-Plus:                   n\n")
    if "Roadmap" in db_list:
        fh.write("Roadmap_sample_ids:                 E001-E129\n")
    fh.close()

def file_as_bytes(file):
    with file:
        return file.read()


args=None
parser = OptionParser()

parser.add_option("--input-minid", dest="input_minid", help="Input minid")
parser.add_option("--tmp-dir", dest="tmpdir", help="tmpdir")
parser.add_option("--wgsa-path", dest="wgsa_path", help="wgsa_path_class")
parser.add_option("--variant-type", dest="vtype", help="Variant type")
parser.add_option("--output", dest="output", help="Output")
parser.add_option("--output-dir", dest="output_dir", help="Output_dir")
parser.add_option( '--history', dest="history_id", help="the history api id" )
parser.add_option("--stdout", dest="output_stdout", help="Output log stdout")
parser.add_option("--stderr", dest="output_stderr", help="Output log stderr")
parser.add_option( '-t', '--token-auth', dest="goauth_token", help='Globus auth token' )
parser.add_option( '-s', '--token-service', dest="goauth_service", help='Globus auth token' )
parser.add_option( '--db', dest="db_type", help="annotation db selected", action="append", type="string" )
parser.add_option( '--xmx', dest='memory_req', help="")

options, args = parser.parse_args(args)

#cwl_file = options.cwl_file
#cwl_file = "/home/galaxy/sbg_dockstore_tools/topmed-workflows/alignment/topmed-alignment.cwl"
cwl_file = "/home/galaxy/topmed-workflows/aligner/sbg-alignment-cwl/topmed-alignment.cwl"
#inputs = json.loads(json_acceptable_string)
output = options.output
inputs_tmp = "%s/input" % options.tmpdir
outputs_tmp = "%s/output" % options.tmpdir
tmp_tmp = "%s/tmp" % options.tmpdir
work_tmp = "%s/work" % options.tmpdir

#if not os.path.exists(output_dir):
#    os.makedirs(output_dir)

# work on the cwl-worflow


# work on the json inputs file

def decide_path_type(path):
    if path.startswith("ark:"):
        return "minid"
    elif path.startswith("s3:"):
        return "s3"
    else:
        msg = "Path type not supported: {0}".format(path)
        sys.exit(msg)

def download_file_from_minid(minid, local_dir):
    QUERY_BASE = None
    if len(minid) > 17:
        QUERY_BASE = "https://identifiers.globus.org"
    else:
        QUERY_BASE = "http://minid.bd2k.org/minid/landingpage"
    minid = minid
    BASE_DOWNLOAD_PATH = os.path.join(local_dir, str(int(time.time()*10000)))

    if not os.path.exists(BASE_DOWNLOAD_PATH):
        os.makedirs(BASE_DOWNLOAD_PATH)

    query = "%s/%s" % (QUERY_BASE, minid)
    #print "Executing query: %s" % query

    r = requests.get(query,  headers = {"Accept" : "application/json"})
    #print r.json()

    location = None
    if len(minid) > 17:
        location = r.json()["location"][0]
    else:
        location = r.json()["locations"][0]['link']

    filename = location.split("/")[-1]
    path = "%s/%s" % (BASE_DOWNLOAD_PATH, filename)

    #print "Downloading result: %s" % location

    testfile = urllib.URLopener()
    testfile.retrieve(location, path)

    # BDBag tooling doesnt let us unzip into a single dir.. so
    # we use a /base/uuid/uuid
    extract_path = ".".join(path.split(".")[0:-1])
    #output_path = "%s/%s" %(extract_path, filename.split(".")[0])
    output_path = "%s/%s" %(extract_path, os.path.basename(extract_path))
    uncomp_path = "%s/%s" %(extract_path, "real_data")

    print "Extracting bag and resolving fetch: %s" % output_path
    bdbag_api.extract_bag(path, extract_path)
    os.mkdir("%s/data" % output_path)
    bdbag_api.resolve_fetch(output_path, True)

    tmp_path = os.path.join(output_path, "data")
    file_list = os.listdir(tmp_path)

    def find_cram_file(f_l):
        for item in f_l:
            if item.endswith('.vcf.gz') or item.endswith('.vcf') or item.endswith('vcf.tar.gz'):
                return item
        return None

    file_name = find_cram_file(file_list)
    if file_name == None:
        file_name = file_list[0]

    if file_name.endswith('vcf.tar.gz'):
        # uncompress/untar and concatenate all vcf.gz files into
        gtex_id = os.path.basename(file_name).split(".")[0]
        tar = tarfile.open(os.path.join(tmp_path, file_name))
        tar.extractall(uncomp_path)
        tar.close()

        # get list of all vcf.gz files
        matches = []
        for root, dirnames, filenames in os.walk(uncomp_path):
            for filename in fnmatch.filter(filenames, '*.vcf.gz'):
                matches.append(os.path.join(root, filename))

        # concatenate all *.vcf.gz files
        my_vcf = "%s/%s.out.vcf.gz" % (output_path, gtex_id)
        cmd = "vcf-concat %s | gzip -c > %s" % (" ".join(matches), my_vcf)
        os.system(cmd)
        file_name = my_vcf

    if "/" in file_name:
        file_path = file_name
    else:
        file_path = os.path.join(tmp_path, file_name)
    print "VCF file is at %s" % file_path
    return file_path

def download_file_from_s3(s3_path, local_dir):
    tmp1 = s3_path.split('//')
    tmp2 = tmp1[1].split('/', 1)
    bucket = tmp2[0]
    file_path = tmp2[1]
    file_name = tmp1[1].split('/')[-1]

    local_dir = os.path.join(local_dir, str(int(time.time()*10000)))

    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    local_file_path = os.path.join(local_dir, file_name)

    #s3 = boto3.resource('s3', config=Config(signature_version='s3v4'))
    #s3.meta.client.download_file(bucket, file_path, local_file_path)

    # Downloading from a public S3 bucket
    #s3 = boto3.resource('s3')
    #bucket = s3.Bucket(bucket)
    #bucket.download_file(file_path, local_file_path)

    # Downloading from NC
    session = boto3.session.Session(profile_name='ncs3conn')
    client = session.client('sts')
    rolearn = 'arn:aws:iam::600168050588:role/developer_access_gtex'
    #rolearn = 'arn:aws:iam::600168050588:role/developer_access_topmed'
    assumeRoleObject = response = client.assume_role(RoleArn=rolearn, RoleSessionName ='NIH-Test', DurationSeconds=3600 )
    credentials = assumeRoleObject['Credentials']
    s3 = boto3.client('s3',aws_access_key_id = credentials['AccessKeyId'],
            aws_secret_access_key = credentials['SecretAccessKey'],
            aws_session_token = credentials['SessionToken'])
    s3.download_file(Bucket=bucket, Key=file_path, Filename=local_file_path)

    return local_file_path
    
def translate_path(path, input_tmp):
    path_type = decide_path_type(path)
    # minid
    if path_type == "minid":
        file_path = download_file_from_minid(path, input_tmp)
        return file_path
    elif path_type == "s3":
        file_path = download_file_from_s3(path, input_tmp)
        return file_path


local_path = translate_path(options.input_minid, inputs_tmp)
sample_name = os.path.basename(local_path).split(".")[0]
#########

## create settings file for wgsa
setting_file = "%s/settings" % inputs_tmp
create_settings(setting_file, local_path, outputs_tmp, work_tmp, tmp_tmp, options.vtype, options.db_type)

## create bash setting file
setting_bash_file = "%s.sh" % setting_file
setting_bash_tmp_file = "%s.tmp.sh" % setting_file
settings_bash_cmd = ['java', '-Xmx%s' % options.memory_req, '-Xms1g','-cp', options.wgsa_path, "WGSA075", setting_file, "-v", "hg38", "-i", "vcf"]
print "WGSA Setting Command: %s" % " ".join(settings_bash_cmd)
proc = subprocess.Popen(settings_bash_cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
stdout,stderr = proc.communicate("Understand")

## modify bash setting file
copyfile(setting_bash_file, setting_bash_tmp_file)
fhr = open(setting_bash_tmp_file, "r")
fhw = open(setting_bash_file, "w")
for line in fhr:
    if "./" in line:
        line = line.replace("./", "/WGSA/")
    fhw.write(line)
fhw.close()
fhr.close()
os.remove(setting_bash_tmp_file)

if os.path.exists(setting_bash_file):
    # submit the wgsa annotation
    cmd = ["bash", setting_bash_file, ">", options.output_stdout, "2>", options.output_stderr]
    print "WGSA Command: %s" % " ".join(cmd)
    os.system(" ".join(cmd))
    #stderr_reformat = tempfile.NamedTemporaryFile( prefix = "reformat_stderr" ).name
    #process = subprocess.Popen(cmd, shell=True, stderr=open( stderr_reformat, 'wb' ), cwd=os.getcwd())
    #rval = process.wait()

else:
    sys.exit("bash script was not properly created: %s, %s" % (stdout, stderr))

output_list = os.listdir(outputs_tmp)
print output_list
if output_list == []:
    sys.exit("No output file.")
else:
    for i in output_list:
        if " " in i:
            os.rename("%s/%s" % (outputs_tmp, i), "%s/%s" % (outputs_tmp, i.replace(" ", ".")))
    output_list = os.listdir(outputs_tmp)
print output_list

#
### s3 transfer of result files

S3_BUCKET = 'gg-commons'
results_url = "https://results.fair-research.org"
results_path = "/results"
bdbags_url = "https://bags.fair-research.org"
bdbags_path = "/bags"

identifier = str(int(time.time()*10000))
s3_dir_to_use = "{5}_{6}/{0}_{1}_{2}_{3}_{4}".format(sample_name, options.vtype, "_".join(options.db_type).replace("/", "_").replace(" ", "_"), os.path.basename(outputs_tmp).strip()+"_1", identifier, sample_name, identifier)
files_to_bdbag = []
remote_manifest = []

for outf in output_list:
    source_path = os.path.join(outputs_tmp, outf)
    to_loc = "{0}/{1}".format(s3_dir_to_use, outf)
    to_path = "%s/%s" % (results_path, to_loc)
    real_url_path = "%s/%s" % (results_url, to_loc)
    file_size = os.path.getsize(source_path)
    md5sum = hashlib.md5(file_as_bytes(open(source_path, 'rb'))).hexdigest()
    filename = os.path.basename(source_path)

    #s3 transfer

    #s3 = boto3.resource('s3')
    session = boto3.session.Session(profile_name='gs3conn')
    s3_client = session.client("s3", config=Config(signature_version="s3v4"))
    #bucket = s3.Bucket(S3_BUCKET)
    url_path = "https://s3.amazonaws.com/%s/results/%s" % (S3_BUCKET,  to_loc)
    s3_client.upload_file(source_path, S3_BUCKET, "results/%s" % (to_loc))

    remote_manifest.append({'url' : real_url_path, "length": file_size, "filename" : filename, "md5" : md5sum})

#print remote_manifest
tmp_file_path = os.path.join(tmp_tmp, "dbbag_history_info_{0}".format(identifier))
#print tmp_file_path
with open(tmp_file_path, 'w') as fp:
    json.dump(remote_manifest, fp)
bdbag_name = "%s.%s" % (s3_dir_to_use, "outputs.bdbag")
bdbag_path = os.path.join(tmp_tmp, bdbag_name)
if not os.path.isdir(bdbag_path):
    os.makedirs(bdbag_path, 0755)
bdbag = bdbag_api.make_bag(bdbag_path, remote_file_manifest=tmp_file_path)

bdbag_fetch = "%s/fetch.txt" % bdbag_path

ro_metadata = dict()
ro_author_name = "Globus Genomics"
ro_author_orcid = "https://orcid.org/0000-0003-0723-665X"
ro_manifest = ro.init_ro_manifest(author_name=ro_author_name, author_orcid=ro_author_orcid)
ro_annotation_about = list()
ro_annotation_filename = "outputs.metadata.json"
ro_annotation_content_path = "annotations/%s" % ro_annotation_filename

fh = open(bdbag_fetch, "r")
for line in fh:
    line = line.rstrip("\n")
    values = line.split("\t")
    url = values[0]
    filename = values[2]
    if filename.startswith("data/"):
        filename = filename.replace("data/", "")
    ro_annotation_about.append(ro.ensure_payload_path_prefix(filename))
    ro.add_file_metadata(ro_manifest,
                         source_url=url, # url parameter from the RFM
                         bundled_as=ro.make_bundled_as(
                                    folder=os.path.dirname(filename), # filename parameter from the RFM
                                    filename=os.path.basename(filename))) # filename parameter from the RFM
#ro.add_annotation(ro_manifest, ro_annotation_about, uri="urn:uuid:%s" % str(uuid.uuid4()), content=ro_annotation_content_path)
ro_metadata["manifest.json"] = ro_manifest

#ro_metadata[ro_annotation_content_path] = history_info_to_write # ALEX: this is the metadata from Galaxy in dict form

# save the BDBAG with the RO data and archive
ro.serialize_bag_ro_metadata(ro_metadata, bdbag_path)
bdbag.save()
# archive the bag
archived_path = bdbag_api.archive_bag(bdbag_path, "zip")

# upload output bdbag to s3
remote_path = "bags/%s" % ( os.path.basename(archived_path))
#bucket.upload_file(archived_path, remote_path, ExtraArgs={'ACL':'public-read'})
#bucket.upload_file(archived_path, remote_path)
#s3_client.upload_file(source_path, S3_BUCKET, remote_path)
s3_client.upload_file(archived_path, S3_BUCKET, remote_path)

s3_url = "%s/%s" % (bdbags_url, os.path.basename(archived_path))

# get user email and name from auth token
#globus_auth_token = '%s \'' % options.goauth_token
#service_token = options.goauth_token
service_token = check_globus_token(options.goauth_token, "/home/galaxy/.globusgenomics/globus_creds", "/home/galaxy/.globusgenomics", "sulakhe", "identifier")

minid_title = "TOPMED workflow output results for %s" % os.path.basename(source_path)

ic = IdentifierClient('Identifier', base_url='https://identifiers.globus.org/', app_name='WES Service', authorizer=globus_sdk.AccessTokenAuthorizer(service_token))
#checksum = minid_client_api.compute_checksum(archived_path)
md5sum = hashlib.md5(file_as_bytes(open(archived_path, 'rb'))).hexdigest()
# HHxPIZaVDh9u - namespace for test minids
# kHAAfCby2zdn - namespace for production once ready
kwargs = {
        'namespace': 'kHAAfCby2zdn',
        'visible_to': json.dumps(['public']),
        'location': json.dumps([s3_url]),
        'checksums': json.dumps([{
              'function': 'md5',
              'value': md5sum
          }]),
        'metadata': json.dumps({
                'Title': minid_title
        })
}
response_with_minid = ic.create_identifier(**kwargs)
minid = response_with_minid['identifier']

out_dict = { 'minid' : minid, 'BDbag_files' : remote_manifest, 'BDbag_path' : s3_url, 'title': minid_title}
fh_out = open(output, "w")
fh_out.write( json.dumps(out_dict))
fh_out.close()
#copyfile(source_path, out_final)

#os.system("ls -l %s > %s" % (output_dir, output))


# delete all content in the tmp dir
#for the_file in os.listdir(TMP_DIR):
#    file_path = os.path.join(TMP_DIR, the_file)
#    try:
#        if os.path.isfile(file_path):
#            os.remove(file_path)
#        elif os.path.isdir(file_path):
#            rmtree(file_path)
#    except Exception as e:
#        print(e)

"""
{
  "output_name": "output1",
  "output": {
    "path": "/home/ubuntu/test/output1",
    "class": "File"
  },
    "bwa_index": {
        "class": "File",
        "path": "ark:/57799/b91414"
    },
    "dbsnp": {
        "class": "File",
        "path": "ark:/57799/b9wb0k"
    },
    "input_file": {
        "class": "File",
        "path": "ark:/57799/b9rm50"
    },
    "reference_genome": {
        "class": "File",
        "path": "ark:/57799/b9mt4f"
    }
}

hs38DH.fa.tar    hs38DH_fa_tar_bdbag.zip    ark:/57799/b91414
Homo_sapiens_assembly38.dbsnp138.vcf.gz    Homo_sapiens_assembly38_dbsnp138_bdbag.zip    ark:/57799/b9wb0k
NWD176325.0005.recab.cram    NWD176325_0005_recab_cram_bdbag.zip    ark:/57799/b9rm50
hs38DH.fa    hs38DH_fa_bdbag.zip    ark:/57799/b9mt4f

{"bwa_index": {"class": "File", "path": "ark:/57799/b91414" }, "dbsnp": {"class": "File", "path": "ark:/57799/b9wb0k"}, "input_file": {"class": "File", "path": "ark:/57799/b9rm50"}, "reference_genome": {"class": "File", "path": "ark:/57799/b9mt4f"}}

"""
