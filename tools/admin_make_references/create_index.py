import os, sys, glob, argparse, shutil, time, re
from bioblend.galaxy import GalaxyInstance
import requests
import tempfile
from slackclient import SlackClient
import globus_sdk
from datetime import datetime

##
# python /scratch/galaxy/test/create_index.py --input-directory /scratch/go/arodri7/tmp/ --genome-directory /mnt/galaxyIndices2/genomes/ --config-directory /mnt/galaxyIndices2/galaxy/tool-data.new2/ --api-key somekey --url https://eupathdbworkshop.globusgenomics.org --token "sometoken"
###


#### Slack authentication
#S_L_A_C_K_T_O_K_E_N = "<YOUR_TOKEN>"
slack_client = SlackClient(SLACK_TOKEN)

def send_message(channel_id, message):
   slack_client.api_call(
         "chat.postMessage",
         channel=channel_id,
         text=message,
         username='pipeline',
         icon_emoji=':bulb:',
         link_names=True
   )

def genome_exists(name, genome_dir):
    genome_id = "_".join(name.split("_")[1:-1])
    delete = []
    existing = glob.glob("%s/*_%s_*" % (genome_dir, genome_id))
    if len(existing) > 0:
        for i in existing:
            delete.append(os.path.basename(i))
    #print genome_dir
    #print name
    #print genome_id
    #print existing
    return delete

def create_len_files(genome_directory, config_directory, existing):
    genomes = sorted(glob.glob("%s/*_Genome/seq/*.fai" % args.genome_directory))
    for genome in genomes:
        name = ".".join(os.path.basename(genome).split(".")[0:-2])
        if name not in existing:
            shutil.copy(genome, "%s/%s.len" % (config_directory, name))

def snpeff_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s" % (config_directory, tool), "w")
    fh.write("?\tunspecified (?)\n")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            fh.write("%s\t%s\n" % (name, gname))
    fh.close()

def snpeffdb_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s" % (config_directory, tool), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            fh.write("%s\t%s\t%s/%s\n" % (name, gname, "/mnt/galaxyIndices2/genomes/snpeff/data/", name))
    fh.close()

def gff_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s.loc" % (config_directory, tool), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            gff_name = "_".join(name.split("_")[0:-1])
            gff = "%s/%s/annotation/%s.gff" % (genome_directory, name, gff_name)
            fh.write("%s\t%s\t%s\t%s\n" % (name, name, gname, gff))
    fh.close()

def sam_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s_indices.loc" % (config_directory, tool), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            fasta = "%s/%s/seq/%s.fasta" % (genome_directory, name, name)
            fh.write("index\t%s\t%s\n" % (name, fasta))
    fh.close()


def gmap_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s_indices.loc" % (config_directory, tool), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            ref_dir = "%s/%s/%s/%s" % (genome_directory, name, "gsnap", name)
            fh.write("%s\t%s\t%s\t15\tsplicesties,introns,snps\tsnps\t%s\n" % (name, name, gname, ref_dir))
    fh.close()

def star_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s_indices.loc" % (config_directory, tool), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            ref_dir = "%s/%s/%s" % (genome_directory, name, tool)
            fasta = "%s/%s/seq/%s.fasta" % (genome_directory, name, name)
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (name, name, gname, ref_dir, fasta))
    fh.close()

def tool_ref_loc(genome_directory, config_directory, tool, suffix, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s_%s.loc" % (config_directory, tool, suffix), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            ref_dir = "%s/%s/%s/%s" % (genome_directory, name, tool, name)
            fh.write("%s\t%s\t%s\t%s\n" % (name, name, gname, ref_dir))
    fh.close()

def picard_loc(genome_directory, config_directory, tool, existing):
    genomes = sorted(glob.glob("%s/*_Genome" % args.genome_directory))
    fh = open("%s/%s.loc" % (config_directory, tool), "w")
    for genome in genomes:
        name = os.path.basename(genome)
        values = name.split("_")
        gname = "%s_%s_%s" % ("_".join(values[1:-1]), values[0], values[-1])
        if name not in existing:
            fasta = "%s/%s/seq/%s.fasta" % (genome_directory, name, name)
            fh.write("%s\t%s\t%s\t%s\n" % (name, name, gname, fasta))
    fh.close()


# prepare directory for genome
def prepare_genome_dir(tup):
    name, fasta, gff, genome_dir = tup
    picard_jar = "/mnt/galaxyTools/tools/picard/2.7.1/picard.jar"
    if not os.path.exists(genome_dir):
        seq_dir = "%s/seq" % genome_dir
        ann_dir = "%s/annotation" % genome_dir
        build_dirs = [genome_dir, seq_dir, ann_dir]
        for build in build_dirs:
            os.makedirs(build)
        
        # copy fasta and annotation files into directories
        shutil.copy(fasta, seq_dir)
        shutil.copy(gff, ann_dir)

        # create samtools faidx file, dict file
        os.system('samtools faidx %s/%s' % (seq_dir, os.path.basename(fasta)))
        os.system('java -jar  %s CreateSequenceDictionary R=%s/%s O=%s/%s.dict' % (picard_jar, seq_dir, os.path.basename(fasta), seq_dir, "".join(os.path.basename(fasta).split(".")[0:-1])))
        return([name, "%s/%s" % (seq_dir, os.path.basename(fasta)), "%s/%s" % (ann_dir, os.path.basename(gff)), genome_dir])

def is_complete (historyID, gi):
    status = gi.histories.get_status(historyID)
    if status['percent_complete'] == "100":
        return True
    elif status['state'] == 'ok' or status['state'] == 'error':
        return True
    else:
        return False

def get_current_genomes (path):
    #current_genomes = os.listdir(path)
    current_genomes = []
    for i in glob.glob("%s/*-*" % path):
        current_genomes.append(os.path.basename(i))
    return current_genomes

def get_remote_genome_list (ep_name, token, remote_path):
    transfer_client = globus_sdk.TransferClient(authorizer=globus_sdk.AccessTokenAuthorizer(token))
    ep = transfer_client.endpoint_search(ep_name)[0]
    ep_id = ep['id']
    contents = transfer_client.operation_ls(ep_id, path=remote_path)['DATA']
    genome_list = {}
    for i in contents:
        if i['name'].endswith('fasta'):
            p = re.compile('(.*).fasta')
            m = p.match(i['name'])
            if m.group() and m.group(1):
                genome_name = m.group(1)
                fasta = i['name']
                if genome_name in genome_list:
                    genome_list[genome_name]['fasta'] = fasta 
                else:
                    genome_list[genome_name] = { 'fasta' : fasta, 'path' : remote_path }
        elif i['name'].endswith('gff'):
            p2 = re.compile('(.*).gff')
            g = p2.match(i['name'])
            if g.group() and g.group(1):
                genome_name = "%s_Genome" % g.group(1)
                gff = i['name']
                if genome_name in genome_list:
                    genome_list[genome_name]['gff'] = gff
                else:
                    genome_list[genome_name] = { 'gff' : gff, 'path' : remote_path }
    return(genome_list)

def returnNotMatches(a, b):
    a = set(a)
    b = set(b)
    return [list(b - a), list(a - b)]

def transfer_genomes(missing, input_directory, token, ep_name):
    transfer_client = globus_sdk.TransferClient(authorizer=globus_sdk.AccessTokenAuthorizer(token))
    s_ep = transfer_client.endpoint_search(ep_name)[0]
    d_ep = transfer_client.endpoint_search("galaxy#eupathdbstaging_0")[0]

    tdata = globus_sdk.TransferData(transfer_client, s_ep['id'], d_ep['id'])

    for i in missing:
        print i
        tdata.add_item("%s/%s" % (i['path'], i['fasta']), "%s/%s" % (input_directory, os.path.basename(i['fasta'])))    
        tdata.add_item("%s/%s" % (i['path'], i['gff']), "%s/%s" % (input_directory, os.path.basename(i['gff'])))

    transfer_result = transfer_client.submit_transfer(tdata)

    while not transfer_client.task_wait(transfer_result['task_id'], timeout=60):
        status = transfer_client.get_task(transfer_result['task_id'])
        #print status
        nice_status_short_description = None
        if status['nice_status_short_description'] is None:
            nice_status_short_description = "active - no errors"
        else:
            nice_status_short_description = status['nice_status_short_description']

        print("%s\t%s\t%s\t%s\t%s" % (datetime.now().strftime("%Y-%m-%d %H:%M"), transfer_result['task_id'], status['status'], status['bytes_transferred'], nice_status_short_description))

        # if there is a permissions denied error, kill the transfer job and report error
        if status['nice_status_short_description'] == "permission denied":
            transfer_client.cancel_task(transfer_result['task_id'])
            continue

# python /scratch/galaxy/test/create_index.py --input-directory /scratch/go/arodri7/tmp/ --genome-directory /mnt/galaxyIndices2/genomes/ --config-directory /mnt/galaxyIndices2/galaxy/tool-data.new2/ --api-key somekey --url https://eupathdbworkshop.globusgenomics.org --token "sometoken"

parser = argparse.ArgumentParser(description='Run genome indexing process')
parser.add_argument('-i', '--input-directory', help='input_directory where raw FASTA and GFF files are located', required=True)
parser.add_argument('-g', '--genome-directory', help='path where genome indexes will be placed', required=True)
parser.add_argument('-c', '--config-directory', help='path to location files', required=True )
parser.add_argument('-k', '--api-key', help='api key for gg instance', required=True )
parser.add_argument('-u', '--url', help='instance url', required=True ) 
parser.add_argument('-t', '--token', help='globus oauth token', required=True )
args = parser.parse_args()

fh_token = open("/home/galaxy/.globusgenomics/tokens/arodri7", "r")
args.token = fh_token.readline().rstrip("\n")

## get list of FASTA, GFF pairs from remote endpoint ("EuPathDB Data for Globus Galaxy")
ep_name = "EuPathDB Data for Globus Galaxy"
remote_path = "/var/www/Common/apiSiteFilesMirror/globusGenomesShare/"
local_genomes = get_current_genomes(args.genome_directory)
remote_genomes = get_remote_genome_list(ep_name, args.token, remote_path)
missing_genomes, delete_genomes = returnNotMatches(local_genomes, remote_genomes.keys())

print "MISSING Genomes: %s" % missing_genomes
print "NEED TO DELETE Genomes: %s" % delete_genomes

# transfer files
missing = []
for i in missing_genomes:
    if "gff" in remote_genomes[i] and "fasta" in remote_genomes[i]:
        missing.append(remote_genomes[i])
    else:
        print "Source is MISSING some data for %s" % remote_genomes[i]

transfer_genomes(missing, args.input_directory, args.token, ep_name)

## get list of FASTA, GFF pairs
raw_tuples = []
fasta_files = glob.glob("%s/*_Genome.fasta" % args.input_directory)
gff_files = glob.glob("%s/*.gff" % args.input_directory)
for genome_name in missing_genomes:
    if "gff" in remote_genomes[genome_name] and "fasta" in remote_genomes[genome_name]:
        print genome_name
        fasta = "%s/%s" % (args.input_directory, remote_genomes[genome_name]['fasta'])
        print fasta
        gff = "%s/%s" % (args.input_directory, remote_genomes[genome_name]['gff'])
        print gff
        genome_path = "%s/%s" % (args.genome_directory, genome_name)
        if os.path.exists(fasta):
            raw_tuples.append([genome_name, fasta, gff, genome_path])
        else:
            print("Could not file FASTA for genome: %s - %s" % (genome_name, fasta))

## Create genome workspace, fasta index file and fasta dictionary file
new_tuples = []
new_genomes = []
#delete_genomes = []
for tup in raw_tuples:
    name, fasta, gff, genome_dir = tup
    #delete_genomes.extend(genome_exists(name, args.genome_directory))
    new_tuples.append(prepare_genome_dir(tup))
    new_genomes.append(name)
#print new_tuples
#print delete_genomes


#####
# create individual loc files, len files
# do no add for genomes that are to be deleted
# loc file to create:
# picard_index.loc
picard_loc(args.genome_directory, args.config_directory, "picard_index", delete_genomes)

# all_fasta.loc
picard_loc(args.genome_directory, args.config_directory, "all_fasta", delete_genomes)

# bowtie_indices.loc
tool_ref_loc(args.genome_directory, args.config_directory, "bowtie", "indices", delete_genomes)

# bowtie2_indices.loc
tool_ref_loc(args.genome_directory, args.config_directory, "bowtie2", "indices", delete_genomes)

# bwa_index.loc
tool_ref_loc(args.genome_directory, args.config_directory, "bwa", "index", delete_genomes)

# star_indices.loc
star_loc(args.genome_directory, args.config_directory, "star", delete_genomes)

# gmap_indices.loc
gmap_loc(args.genome_directory, args.config_directory, "gmap", delete_genomes)

# sam_fa_indices.loc
sam_loc(args.genome_directory, args.config_directory, "sam_fa", delete_genomes)

# all_gff3.loc
gff_loc(args.genome_directory, args.config_directory, "all_gff3", delete_genomes)

# snpeff_databases.loc
snpeff_loc(args.genome_directory, args.config_directory, "snpeff_databases.loc", delete_genomes)

# hisat2_index.loc
tool_ref_loc(args.genome_directory, args.config_directory, "hisat2", "index", delete_genomes)

# builds.txt
os.makedirs("%s/shared" % args.config_directory)
os.makedirs("%s/shared/ucsc" % args.config_directory)
snpeff_loc(args.genome_directory, args.config_directory, "shared/ucsc/builds.txt", delete_genomes)

# create len file
os.makedirs("%s/len_files" % args.config_directory)
create_len_files(args.genome_directory, "%s/len_files" % args.config_directory, delete_genomes)

# download ga workflow file from github repo
ga_file = "https://raw.githubusercontent.com/arodri7/galaxy_workflows/master/Galaxy-Workflow-eupath_indexes-5.ga"

fd, tmp_ga_file = tempfile.mkstemp()
r = requests.get(ga_file, allow_redirects=True)
with os.fdopen(fd, 'w') as tmp:
    # do stuff with temp file
    tmp.write(r.content)


# submit each genome through workflow
gi = GalaxyInstance(url=args.url, key=args.api_key)
workflow = gi.workflows.import_workflow_from_local_path(tmp_ga_file)
os.remove(tmp_ga_file)
#print workflow
workflow = gi.workflows.show_workflow(workflow['id'])
invocations = []
histories = []
for tup in new_tuples:
    a = {}
    a['ref_id'], a['fasta_input'], a['gff_input'], a['genome_dir'] = tup
    a['ref_name'] = a['ref_id']
    a['prefix'] = a['ref_id']
    a['build_dir'] = None

    ds_map = {}
    parameters = {}
    history_name = "Batch_index-%s-%s" % (a['ref_id'], time.strftime("%a_%b_%d_%Y_%-I:%M:%S_%p",time.localtime(time.time())))
    history = gi.histories.create_history(name=history_name)
    histories.append(history['id'])
    for key,val in workflow['steps'].iteritems():
        parameters[str(key)] = {}
        if val['tool_id'] == 'delete_history':
            parameters[str(key)] = {'historyid' : history['id'] }
        else:
            tool = val['tool_id'].split("_")[0]
            if tool == "gmap":
                tool = "gsnap"
        a['build_dir'] = "%s/%s" % (a['genome_dir'], tool)
        a['ref_path'] = "%s/%s/%s" % (a['genome_dir'], tool, a['ref_id'])
        for step in val['tool_inputs']:
            if isinstance(val['tool_inputs'][step], dict) and '__class__' in val['tool_inputs'][step] and 'RuntimeValue' == val['tool_inputs'][step]['__class__'] and step != 'input':
                #print "%s: %s" % (step, val['tool_inputs'][step])
                parameters[str(key)][step] = a[step]
    #print parameters

    # submit batch file
    invocations.append(gi.workflows.invoke_workflow(workflow['id'], inputs=ds_map, params=parameters, history_id=history['id'], import_inputs_to_history=False))

#print invocations

# wait for batch job to complete
completed_meta = {}
while len(completed_meta) != len(histories):
    for history_id in histories:
        if history_id in completed_meta.keys():
            continue
        elif is_complete(history_id, gi) == True:
            completed_meta[history_id] = {'status' : 'complete'}
    if len(histories) != len(completed_meta):
        time.sleep(300)

# send message via slack channel
## send slack message
message = "@sulakhe @Alex @ravi @markxiao\nCongratulations! Indexes have been built on URL: %s\nThe following genomes are new:\n %s.\n\nThe following genomes need to be deleted from %s:\n %s\n\n:beers: :thumbsup: :facepunch:\n\nNow we need to send over to Rsync, update loc files and restart!" % (args.url, "\n".join(new_genomes), args.genome_directory, "\n".join(delete_genomes) )
send_message('globus-genomics', message)
print 'Slack Message has been sent'
