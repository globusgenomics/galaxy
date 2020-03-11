import logging
import os
import time
import urllib
from galaxy.tools.actions import DefaultToolAction
import globus_sdk
import boto3
import botocore

log = logging.getLogger( __name__ )


class GlobusTransferInAction( DefaultToolAction ):

    def generate_transfer_info(self, tool_id, incoming, galaxy_config, username):
        transfer_job_info = []
        if tool_id == 'globus_get_data_flowcell_text':
            for i in range(1,9):
                from_path_key = 'from_path' + str(i)
                if from_path_key in incoming and incoming[from_path_key] != '':
                    job = {}
                    job['from_path'] = incoming[from_path_key].strip().rstrip('/').rstrip('\\').replace(':\\', '/').replace('\\', '/')
                    job['to_path'] = os.path.join(galaxy_config['ftp_upload_dir'], username, str(int(time.time()*10000)), job['from_path'].lstrip('/~').replace(' ', '_'))
                    job['output_name'] = 'out_file' + str(i)
                    transfer_job_info.append(job)

        elif tool_id == 'globus_get_data_text':
            job = {}
            job['from_path'] = incoming['from_path'].strip().rstrip('/').rstrip('\\').replace(':\\', '/').replace('\\', '/')
            job['to_path'] = os.path.join(galaxy_config['ftp_upload_dir'], username, str(int(time.time()*10000)), job['from_path'].lstrip('/~').replace(' ', '_'))
            job['output_name'] = 'out_file1'
            transfer_job_info.append(job)

        elif tool_id == 'globus_get_bioproject_data_from_ebi':
            timestamp = str(int(time.time()*10000))
            accessions = incoming["accession"].split(",")
            url_path = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=fastq_ftp' % incoming["accession"]
            response = urllib.urlopen(url_path)
            html = response.read()
            reads = html.split("\n")[1:-1]
            for read_set in reads:
              for read in read_set.split(";"):
                job = {}
                path_values = read.split("/")
                path_values.pop(0)
                from_path = "/" + "/".join(path_values)
                job["from_path"] = from_path.replace("/vol1/", "/gridftp/ena/")
                job["to_path"] = os.path.join(galaxy_config['ftp_upload_dir'], username, timestamp, os.path.basename(job['from_path']).lstrip('/~')).replace(' ', '_')
                job['output_name'] = 'out_file1'
                job['output_ext'] = "fastqsanger"
                transfer_job_info.append(job)

        elif tool_id == 'globus_get_data_from_ebi' or tool_id == 'globus_get_data_from_ebi_to_collections':
            timestamp = str(int(time.time()*10000))
            accessions = incoming["accession"].split(",")
            for acc in accessions:
                url_path = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=fastq_ftp' % acc
                response = urllib.urlopen(url_path)
                html = response.read()
                if incoming['datatype_cond']['single_paired'] == "paired":
                    job0 = {}
                    job1 = {}
                    ebi_paths = html.split("\n")[1].split(";")
                    paired_paths = []
                    for path in ebi_paths:
                        path_values = path.split("/")
                        path_values.pop(0)
                        paired_paths.append("/" + "/".join(path_values))
                    job0["from_path"] = paired_paths[0].replace("/vol1/", "/gridftp/ena/")
                    job1["from_path"] = paired_paths[1].replace("/vol1/", "/gridftp/ena/")
                    job0["to_path"] = os.path.join(galaxy_config['ftp_upload_dir'], username, timestamp, os.path.basename(job0['from_path']).lstrip('/~')).replace(' ', '_')
                    job1["to_path"] = os.path.join(galaxy_config['ftp_upload_dir'], username, timestamp, os.path.basename(job1['from_path']).lstrip('/~')).replace(' ', '_')
                    job0['output_name'] = 'out_file1'
                    job1['output_name'] = 'out_file2'
                    if incoming['datatype_cond']['datatype'] == "fastqsanger":
                        job0['output_ext'] = "fastqsanger"
                    transfer_job_info.append(job0)
                    transfer_job_info.append(job1)
                else:
                    job = {}
                    path_values = html.split("\n")[1].split("/")
                    path_values.pop(0)
                    job["from_path"] = ("/" + "/".join(path_values)).replace("/vol1/", "/gridftp/ena/") 
                    job["to_path"] = os.path.join(galaxy_config['ftp_upload_dir'], username, timestamp, os.path.basename(job['from_path']).lstrip('/~')).replace(' ', '_').replace("/vol1/", "/gridftp/ena/")
                    job['output_name'] = 'out_file1'
                    if incoming['datatype_cond']['datatype'] == "fastqsanger":
                        job['output_ext'] = "fastqsanger"
                    transfer_job_info.append(job)
        else:
            raise Exception('Missing tool id')

        return transfer_job_info
        


    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):

        tmp_file_dir = '/home/galaxy/.globusgenomics/tmp'

        identifier = str(int(time.time()*10000))

        galaxy_config = trans.app.config.config_dict

        if 'username' in incoming and incoming['username'] not in ['**', '', None]:
            username = incoming['username']
        else:
            #username = trans.user.username
            user_email = trans.user.email
            username = user_email[0:user_email.find('@')]
   
        if 'goauth' in incoming and incoming['goauth'] not in ['**', '', None]:
            goauth_token = incoming['goauth']
        elif hasattr(trans, 'environ') and 'X-GLOBUS-TOKEN' in trans.environ:
            goauth_token = trans.environ['X-GLOBUS-TOKEN']
        elif os.path.isfile(os.path.join(galaxy_config['globus_dir'], 'tokens', username)):
            with open(os.path.join(galaxy_config['globus_dir'], 'tokens', username), 'r') as f:
                goauth_token = f.readline().strip()
        else:
            goauth_token = None        

        transfer_job_info = self.generate_transfer_info(tool.id, incoming, galaxy_config, username)

        if 'path_hidden' in incoming and incoming['path_hidden']:
            job_to_path = transfer_job_info[0]["to_path"]
            incoming['path_hidden'] = os.path.dirname(job_to_path)


        transfer_info = {'username': username,
                         'globus_cred_file': galaxy_config['globus_cred_file'].strip(),
                         'goauth_token': goauth_token, 
                         'from_endpoint': incoming['from_endpoint'].strip(),
                         'to_endpoint': galaxy_config['globus_endpoint'].strip(),
                         'jobs': transfer_job_info}

        # record transfer_info
        if not os.path.isdir(tmp_file_dir):
            os.mkdir(tmp_file_dir, 0700)

        transfer_info_file_name = 'globus_transfer_info_{0}_{1}'.format(username, identifier)
        transfer_info_file_loc = os.path.join(tmp_file_dir, transfer_info_file_name) 

        with open(transfer_info_file_loc, 'w') as f:
            f.write(str(transfer_info))
        
        incoming['transfer_info'] = transfer_info_file_name
        
        if 'transfer_direction' not in incoming:
            incoming['transfer_direction'] = 'in'

        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)

        transfer_client = globus_sdk.TransferClient(authorizer=globus_sdk.AccessTokenAuthorizer(transfer_info['goauth_token']))
        try:
            for ep in transfer_client.endpoint_search(transfer_info['from_endpoint']):
                if ep["display_name"] == transfer_info['from_endpoint'] or ep['canonical_name'] == transfer_info['from_endpoint']:
                    ep_id = ep['id']
                    r = transfer_client.endpoint_autoactivate(ep_id, if_expires_in=3600)
        except:
            pass

        if tool.id == "globus_get_bioproject_data_from_ebi":
            hda = output['out_file1']
            transfer_type = 'dir'
            dataset = hda.dataset
            if not os.path.exists(os.path.dirname(transfer_info['jobs'][0]["to_path"])):
                os.makedirs(os.path.dirname(transfer_info['jobs'][0]["to_path"]))
                symlink_dir = str(dataset.get_file_name())[0:-4] + "_files"
                os.symlink(os.path.dirname(transfer_info['jobs'][0]["to_path"]), symlink_dir)
            if 'symlink' in incoming and incoming['symlink'] == "**":
                incoming['symlink'] = symlink_dir
            ext = "txt"
            hda.extension = ext
            hda.name = "Get BioProject data %s from EBI server" % incoming["accession"]
            info = "Globus transfer summary:\n From: {0}\n To: {1}\n".format(transfer_info['from_endpoint'], transfer_info['to_endpoint'])
            info = info + 'Directory Object ({0})'.format(os.path.dirname(transfer_info['jobs'][0]["to_path"]))
            hda.info = info
        else:
            for i in transfer_info['jobs']:
                # Check transfer type
                transfer_type = 'file'
                try:
                    transfer_client.operation_ls(ep_id, path=i['from_path'])
                    transfer_type = 'dir'
                except:
                    pass
    
                if tool.id != 'globus_get_data_from_ebi_to_collections' and tool.id != "globus_get_bioproject_data_from_ebi":
                    hda = output[i['output_name']]
                    dataset = hda.dataset
    
                    if transfer_type == 'dir':
                        if not os.path.exists(i["to_path"]):
                            os.makedirs(i["to_path"])
                            symlink_dir = str(dataset.get_file_name())[0:-4] + "_files"
                            msg = "!!!!!!!!!!!!!!!!!os.symlink(i[\"to_path\"], symlink_dir): {0}, {1}".format(i["to_path"],symlink_dir)
                            log.info(msg)
                            os.symlink(i["to_path"], symlink_dir)
                        ext = "txt"
                    else:
                        if 'output_ext' in i:
                            ext = i['output_ext']
                        elif "fastq.gz" in dataset.file_name:
                            ext = "fastqsanger"   
                        else:
                            ext = os.path.splitext(i["to_path"])[-1]
                            if ext == ".gz":
                                ext = "fastqsanger"
                            else:
                                ext = ext[1:] # remove the leading dot.
                        dataset.set_file_name(i["to_path"])

                    hda.extension = ext
                    hda.name = os.path.basename(i["to_path"])
                    info = "Globus transfer summary:\n From: {0}\n To: {1}\n".format(transfer_info['from_endpoint'], transfer_info['to_endpoint'])

                    if transfer_type == 'dir':
                        info = info + 'Directory Object ({0})'.format(i["to_path"])
                    hda.info = info

        trans.sa_session.flush()

        return job, output


class GlobusTransferOutAction( DefaultToolAction ):

    def get_file_path_from_dataset(self, dataset):
        time.sleep(30)
        dataset_path = dataset.get_file_name()
        msg = "!!!!!!!!!!!!!!!!X!!dataset_path: {0}".format(dataset_path)
        log.info(msg)

        # Check whether it is a Directory object thansferred through Globus Genomics
        dir_link = dataset_path[0:-4] + "_files"
        if os.path.islink(dir_link):
            dir_realpath = os.path.realpath(dir_link)
            dataset_path = dir_realpath
        elif os.path.isdir(dir_link):
            dataset_path = dir_link

        dataset_path = dataset_path.rstrip('/').rstrip('\\')
        return dataset_path

    def generate_transfer_info(self, tool_id, incoming):
        transfer_job_info = []
        
        def handle_send_bam(dataset, to_path):
            if str(dataset.ext) == 'bam':
                if dataset.metadata.bam_index not in ['None', ''] and hasattr(dataset.metadata.bam_index, 'file_name'):
                    if str(dataset.metadata.bam_index.file_name) not in ['None', '']: 
                        print dataset.metadata.bam_index.file_name
                        job = {}
                        job['from_path'] = str(dataset.metadata.bam_index.file_name)
                        job['to_path'] = to_path + '.bai'
                        transfer_job_info.append(job)

        if tool_id == 'globus_send_data_multiple':
            for i in incoming['src_dataset']:
                job = {}
                job['from_path'] = self.get_file_path_from_dataset(i['from_path'].dataset)
                job['to_path'] = i['to_path'].strip().rstrip('/').rstrip('\\').replace(' ', '_').replace(':\\', '/').replace('\\', '/')
                msg = "!!!!!!!!!!!!!!!!X!!globus_send_data_multiple: {0}, {1}".format(job['from_path'], job['to_path'])
                log.info(msg)
                transfer_job_info.append(job)
                handle_send_bam(i['from_path'], job['to_path'])

        elif tool_id == 'globus_send_data_multiple2':
            stamp = time.strftime('_%Y_%m_%d', time.localtime(time.time()))
            for i in incoming['src_dataset']:
                job = {}
                job['from_path'] = self.get_file_path_from_dataset(i['from_path'].dataset)
                if i["include_datestamp"] == "none":
                    stamp = ""
                to_path = "%s/%s%s.%s" % (incoming['to_directory'], i['to_name'], stamp, i['to_extension'])
                job['to_path'] = to_path.strip().rstrip('/').rstrip('\\').replace(' ', '_').replace(':\\', '/').replace('\\', '/')
                msg = "!!!!!!!!!!!!!!!!X!!globus_send_data_multiple2: {0}, {1}".format(job['from_path'], job['to_path'])
                log.info(msg)
                transfer_job_info.append(job)
                handle_send_bam(i['from_path'], job['to_path'])

        elif tool_id == 'globus_send_data':
            job = {}
            job['from_path'] = self.get_file_path_from_dataset(incoming['from_dataset'].dataset)
            job['to_path'] = incoming['to_path'].strip().rstrip('/').rstrip('\\').replace(' ', '_').replace(':\\', '/').replace('\\', '/')
            msg = "!!!!!!!!!!!!!!!!!globus_send_data: {0}, {1}".format(job['from_path'], job['to_path'])
            log.info(msg)
            transfer_job_info.append(job)
            #handle_send_bam(incoming['from_dataset'], job['to_path'])

        elif tool_id == 'globus_send_data2':
            job = {}
            job['from_path'] = self.get_file_path_from_dataset(incoming['from_dataset'].dataset)
            stamp = ""
            if incoming["include_datestamp"] == "include":
                stamp = time.strftime('_%Y_%m_%d', time.localtime(time.time()))

            to_path = "%s/%s%s.%s" % (incoming['to_directory'], incoming['to_name'], stamp, incoming['to_extension'])
            job['to_path'] = to_path.strip().rstrip('/').rstrip('\\').replace(' ', '_').replace(':\\', '/').replace('\\', '/')
            msg = "!!!!!!!!!!!!!!!!!globus_send_data2: {0}, {1}".format(job['from_path'], job['to_path'])
            log.info(msg)
            transfer_job_info.append(job)
            #handle_send_bam(incoming['from_dataset'], job['to_path'])  

        else:
            raise Exception('Missing tool id')

        return transfer_job_info


    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):
        tmp_file_dir = '/home/galaxy/.globusgenomics/tmp'

        identifier = str(int(time.time()*10000))

        galaxy_config = trans.app.config.config_dict

        if 'username' in incoming and incoming['username'] not in ['**', '', None]:
            username = incoming['username']
        else:
            #username = trans.user.username
            user_email = trans.user.email
            username = user_email[0:user_email.find('@')]

        if 'goauth' in incoming and incoming['goauth'] not in ['**', '', None]:
            goauth_token = incoming['goauth']
        elif hasattr(trans, 'environ') and 'X-GLOBUS-TOKEN' in trans.environ:
            goauth_token = trans.environ['X-GLOBUS-TOKEN']
        elif os.path.isfile(os.path.join(galaxy_config['globus_dir'], 'tokens', username)):
            with open(os.path.join(galaxy_config['globus_dir'], 'tokens', username), 'r') as f:
                goauth_token = f.readline().strip()
        else:
            goauth_token = None   

        transfer_job_info = self.generate_transfer_info(tool.id, incoming)

        transfer_info = {'username': username, 
                         'globus_cred_file': galaxy_config['globus_cred_file'].strip(),
                         'goauth_token': goauth_token, 
                         'from_endpoint': galaxy_config['globus_endpoint'].strip(),
                         'to_endpoint': incoming['to_endpoint'].strip(),
                         'jobs': transfer_job_info}

        # record transfer_info
        if not os.path.isdir(tmp_file_dir):
            os.mkdir(tmp_file_dir, 0700)

        transfer_info_file_name = 'globus_transfer_info_{0}_{1}'.format(username, identifier)
        transfer_info_file_loc = os.path.join(tmp_file_dir, transfer_info_file_name) 

        with open(transfer_info_file_loc, 'w') as f:
            f.write(str(transfer_info))
        
        incoming['transfer_info'] = transfer_info_file_name
        
        if 'transfer_direction' not in incoming:
            incoming['transfer_direction'] = 'out'

        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)
        
        hda = output["out_file1"]
        hda.extension = "txt"
        hda.name = "Send to {0}".format(transfer_info['to_endpoint'])
        hda.info = "Globus transfer summary:\n From: {0}\n To: {1}\n".format(transfer_info['from_endpoint'], transfer_info['to_endpoint'])

        trans.sa_session.flush()
        
        return job, output


class S3Transfer(DefaultToolAction):
    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):

        cred_file_dir = '/home/galaxy/.globusgenomics/tmp'

        if 'username' in incoming and incoming['username'] not in ['**', '', None]:
            username = incoming['username']
        else:
            username = trans.user.username

        identifier = str(int(time.time()*10000))

        incoming['bucket'] = incoming['bucket'].strip().strip('/').strip('\\')
        incoming['from_path'] = incoming['from_path'].strip().strip('/').rstrip('\\')
        incoming['to_path'] = os.path.join(trans.app.config.ftp_upload_dir, username, identifier, incoming['bucket'], incoming['from_path'].lstrip('/~')).rstrip('/').replace(' ', '_')

        if incoming['from_path'] == '':
            object_type = 'bucket'
        else:
            object_type = 'file'
            try:
                s3 = boto3.resource('s3', aws_access_key_id=incoming['aws_access_key_id'],
                                    aws_secret_access_key=incoming['aws_secret_access_key'])
                bucket = s3.Bucket(incoming['bucket'])
                prefix = incoming['from_path'] + '/'
                #bucket.Object(prefix).load()
                for obj in bucket.objects.filter(Prefix=prefix):
                    object_type = 'dir'
                    break
            except botocore.exceptions.ClientError as e:
                print '!!!!!!!!!!!boto'
                print e
                print vars(e)
                pass

        incoming['object_type'] = object_type

        # record creds
        if not os.path.isdir(cred_file_dir):
            os.mkdir(cred_file_dir, 0700)

        cred_file_name = 's3_transfer_key_{0}_{1}'.format(username, identifier)
        cred_file_loc = os.path.join(cred_file_dir, cred_file_name) 

        with open(cred_file_loc, 'w') as f:
            if 'aws_access_key_id' in incoming:
                f.write('{0}\n'.format(incoming['aws_access_key_id']))
                incoming['aws_access_key_id'] = '**'
            else:
                f.write('\n')
            if 'aws_secret_access_key' in incoming:
                f.write('{0}'.format(incoming['aws_secret_access_key']))
                incoming['aws_secret_access_key'] = '**'
            else:
                f.write('')

        incoming['cred_file'] = cred_file_name

        transfer_info = {
            'bucket': incoming['bucket'],
            'from_path': incoming['from_path'],
            'to_path': incoming['to_path'],
            'cred_file': incoming['cred_file'],
            'object_type': incoming['object_type'],
            'include_subdir': incoming['include_subdir'],
            'tags': incoming['tags']
        }
        incoming['transfer_info'] = transfer_info

        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)
        
        hda = output["out_file1"]
        dataset = hda.dataset

        if incoming['object_type'] in ['dir', 'bucket']:
            if not os.path.exists(incoming['to_path']):
                os.makedirs(incoming['to_path'])
                symlink_dir = str(dataset.get_file_name())[0:-4] + "_files"
                os.symlink(incoming['to_path'], symlink_dir)
            ext = "txt"
        else:
            if "fastq.gz" in dataset.file_name:
                ext = "fastqsanger"   
            else:
                ext = os.path.splitext(incoming["to_path"])[-1]
                if ext == ".gz":
                    ext = "fastqsanger"
                else:
                    ext = ext[1:] # remove the leading dot.
            dataset.set_file_name(incoming["to_path"])

        hda.extension = ext
        hda_name = os.path.basename(incoming["to_path"])
        if hda_name == identifier:
            hda.name = ''
        else:
            hda.name = hda_name
        info = "S3 transfer\n from: {0}\n path: {1}\n".format(incoming['bucket'], incoming['from_path'])
        if incoming['object_type'] == 'dir':
            info = info + "Directory Object\n"
        elif incoming['object_type'] == 'bucket':
            info = info + "Bucket Object\n"
        hda.info = info

        trans.sa_session.flush()

        return job, output


class S3TransferOut(DefaultToolAction):

    def get_file_path_from_dataset(self, dataset):
        dataset_path = dataset.get_file_name()
        object_type = 'file'
        # Check whether it is a Directory object thansferred through Globus Genomics
        dir_link = dataset_path[0:-4] + "_files"
        if os.path.islink(dir_link):
            dir_realpath = os.path.realpath(dir_link)
            dataset_path = dir_realpath
            object_type = 'dir'
        elif os.path.isdir(dir_link):
            dataset_path = dir_link
            object_type = 'dir'

        dataset_path = dataset_path.rstrip('/').rstrip('\\')
        return (dataset_path, object_type)

    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):

        cred_file_dir = '/home/galaxy/.globusgenomics/tmp'

        if 'username' in incoming and incoming['username'] not in ['**', '', None]:
            username = incoming['username']
        else:
            username = trans.user.username

        identifier = str(int(time.time()*10000))

        # record creds
        if not os.path.isdir(cred_file_dir):
            os.mkdir(cred_file_dir, 0700)

        cred_file_name = 's3_transfer_key_{0}_{1}'.format(username, identifier)
        cred_file_loc = os.path.join(cred_file_dir, cred_file_name) 

        with open(cred_file_loc, 'w') as f:
            if 'aws_access_key_id' in incoming:
                f.write('{0}\n'.format(incoming['aws_access_key_id']))
                incoming['aws_access_key_id'] = '**'
            else:
                f.write('\n')
            if 'aws_secret_access_key' in incoming:
                f.write('{0}'.format(incoming['aws_secret_access_key']))
                incoming['aws_secret_access_key'] = '**'
            else:
                f.write('')

        incoming['cred_file'] = cred_file_name

        incoming['bucket'] = incoming['bucket'].strip().strip('/').strip('\\')
        incoming['to_path'] = incoming['to_path'].strip().strip('/').rstrip('\\')
        dataset_info = self.get_file_path_from_dataset(incoming['from_dataset'].dataset)
        incoming['from_path'] = dataset_info[0]
        incoming['object_type'] = dataset_info[1]

        transfer_info = {
            'bucket': incoming['bucket'],
            'from_path': incoming['from_path'],
            'to_path': incoming['to_path'],
            'cred_file': incoming['cred_file'],
            'object_type': incoming['object_type'],
            'rename': incoming['rename'],
            'sse': incoming['sse'],
            'tags': incoming['tags']
        }
        incoming['transfer_info'] = transfer_info

        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)
        
        hda = output["out_file1"]
        hda.extension = "txt"
        hda.name = "Send to {0}".format(incoming['bucket'])
        info = "S3 transfer\n to: {0}\n path: {1}\n".format(incoming['bucket'], incoming['to_path'])
        if incoming['object_type'] == 'dir':
            info = info + "Directory Object\n"
        hda.info = info

        trans.sa_session.flush()
        
        return job, output


class S3TransferOutMultiple(DefaultToolAction):

    def get_file_path_from_dataset(self, dataset):
        dataset_path = dataset.get_file_name()
        object_type = 'file'
        # Check whether it is a Directory object thansferred through Globus Genomics
        dir_link = dataset_path[0:-4] + "_files"
        if os.path.islink(dir_link):
            dir_realpath = os.path.realpath(dir_link)
            dataset_path = dir_realpath
            object_type = 'dir'
        elif os.path.isdir(dir_link):
            dataset_path = dir_link
            object_type = 'dir'

        dataset_path = dataset_path.rstrip('/').rstrip('\\')
        return (dataset_path, object_type)

    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):

        cred_file_dir = '/home/galaxy/.globusgenomics/tmp'

        if 'username' in incoming and incoming['username'] not in ['**', '', None]:
            username = incoming['username']
        else:
            username = trans.user.username

        identifier = str(int(time.time()*10000))

        # record creds
        if not os.path.isdir(cred_file_dir):
            os.mkdir(cred_file_dir, 0700)

        cred_file_name = 's3_transfer_key_{0}_{1}'.format(username, identifier)
        cred_file_loc = os.path.join(cred_file_dir, cred_file_name) 

        with open(cred_file_loc, 'w') as f:
            if 'aws_access_key_id' in incoming:
                f.write('{0}\n'.format(incoming['aws_access_key_id']))
                incoming['aws_access_key_id'] = '**'
            else:
                f.write('\n')
            if 'aws_secret_access_key' in incoming:
                f.write('{0}'.format(incoming['aws_secret_access_key']))
                incoming['aws_secret_access_key'] = '**'
            else:
                f.write('')

        incoming['cred_file'] = cred_file_name

        incoming['bucket'] = incoming['bucket'].strip().strip('/').strip('\\')

        transfer_job_info = []
        for i in incoming['src_dataset']:
            job = {}
            dataset_info = self.get_file_path_from_dataset(i['from_dataset'].dataset)
            job['from_path'] = dataset_info[0]
            job['object_type'] = dataset_info[1]
            job['to_path'] = i['to_path'].strip().strip('/').rstrip('\\')
            job['rename'] = i['rename']
            transfer_job_info.append(job)

        transfer_info = {
            'bucket': incoming['bucket'],
            'cred_file': incoming['cred_file'],
            'sse': incoming['sse'],
            'tags': incoming['tags'],
            'jobs': transfer_job_info
        }
        incoming['transfer_info'] = transfer_info

        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)
        
        hda = output["out_file1"]
        hda.extension = "txt"
        hda.name = "Send to {0}".format(incoming['bucket'])
        info = "S3 transfer to: {0}\n".format(incoming['bucket'])
        hda.info = info

        trans.sa_session.flush()
        
        return job, output

class S3TransferOptimized(DefaultToolAction):
    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):

        cred_file_dir = '/home/galaxy/.globusgenomics/tmp'

        if 'username' in incoming and incoming['username'] not in ['**', '', None]:
            username = incoming['username']
        else:
            username = trans.user.username

        identifier = str(int(time.time()*10000))

        incoming['bucket'] = incoming['bucket'].strip().strip('/').strip('\\')
        incoming['from_path'] = incoming['from_path'].strip().strip('/').rstrip('\\')
        incoming['to_path'] = os.path.join(trans.app.config.ftp_upload_dir, username, identifier, incoming['bucket'], incoming['from_path'].lstrip('/~')).rstrip('/').replace(' ', '_')

        if incoming['from_path'] == '':
            object_type = 'bucket'
        else:
            object_type = 'file'
            try:
                s3 = boto3.resource('s3', aws_access_key_id=incoming['aws_access_key_id'],
                                    aws_secret_access_key=incoming['aws_secret_access_key'])
                bucket = s3.Bucket(incoming['bucket'])
                prefix = incoming['from_path'] + '/'
                #bucket.Object(prefix).load()
                for obj in bucket.objects.filter(Prefix=prefix):
                    object_type = 'dir'
                    break
            except botocore.exceptions.ClientError as e:
                print '!!!!!!!!!!!boto'
                print e
                print vars(e)
                pass

        incoming['object_type'] = object_type

        # record creds
        if not os.path.isdir(cred_file_dir):
            os.mkdir(cred_file_dir, 0700)

        cred_file_name = 's3_transfer_key_{0}_{1}'.format(username, identifier)
        cred_file_loc = os.path.join(cred_file_dir, cred_file_name)

        with open(cred_file_loc, 'w') as f:
            if 'aws_access_key_id' in incoming:
                f.write('{0}\n'.format(incoming['aws_access_key_id']))
                incoming['aws_access_key_id'] = '**'
            else:
                f.write('\n')
            if 'aws_secret_access_key' in incoming:
                f.write('{0}'.format(incoming['aws_secret_access_key']))
                incoming['aws_secret_access_key'] = '**'
            else:
                f.write('')

        incoming['cred_file'] = cred_file_name

        transfer_info = {
            'bucket': incoming['bucket'],
            'from_path': incoming['from_path'],
            'to_path': incoming['to_path'],
            'cred_file': incoming['cred_file'],
            'object_type': incoming['object_type'],
            'include_subdir': incoming['include_subdir'],
            'tags': incoming['tags']
        }
        incoming['transfer_info'] = transfer_info

        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)

        hda = output["out_file1"]
        dataset = hda.dataset

        if incoming['object_type'] in ['dir', 'bucket']:
            if not os.path.exists(incoming['to_path']):
                os.makedirs(incoming['to_path'])
                symlink_dir = str(dataset.get_file_name())[0:-4] + "_files"
                os.symlink(incoming['to_path'], symlink_dir)
            ext = "txt"
        else:
            if "fastq.gz" in dataset.file_name:
                ext = "fastqsanger"
            else:
                ext = os.path.splitext(incoming["to_path"])[-1]
                if ext == ".gz":
                    ext = "fastqsanger"
                else:
                    ext = ext[1:] # remove the leading dot.
            dataset.set_file_name(incoming["to_path"])

        hda.extension = ext
        hda_name = os.path.basename(incoming["to_path"])
        if hda_name == identifier:
            hda.name = ''
        else:
            hda.name = hda_name
        info = "S3 transfer\n from: {0}\n path: {1}\n".format(incoming['bucket'], incoming['from_path'])
        if incoming['object_type'] == 'dir':
            info = info + "Directory Object\n"
        elif incoming['object_type'] == 'bucket':
            info = info + "Bucket Object\n"
        hda.info = info

        trans.sa_session.flush()

        return job, output

class BatchSubmit(DefaultToolAction):
    """This tool basically just works like a globus client."""
    GOAUTH_TOKEN_KEY = 'X-GLOBUS-TOKEN'
    REMOTE_USER_KEY = 'X-GLOBUS-USER'

    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):

        galaxy_config = trans.app.config.config_dict

        import socket
        incoming["userapi"] = trans.user.api_keys[0].key
        incoming["url"] = "https://%s" % socket.gethostname()
        username = trans.user.username
        if "username" in incoming and incoming['username'] not in ['**', '', None]:
            if self.REMOTE_USER_KEY in trans.environ:
                incoming['username'] = trans.environ[self.REMOTE_USER_KEY]
                log.debug('Using remote user %s' % incoming['username'])
            else:
                incoming["username"] = username

        if 'goauth' in incoming and incoming['goauth'] not in ['**', '', None]:
            goauth_token = incoming['goauth']
        elif hasattr(trans, 'environ') and 'X-GLOBUS-TOKEN' in trans.environ:
            goauth_token = trans.environ['X-GLOBUS-TOKEN']
        elif os.path.isfile(os.path.join(galaxy_config['globus_dir'], 'tokens', username)):
            with open(os.path.join(galaxy_config['globus_dir'], 'tokens', username), 'r') as f:
                goauth_token = f.readline().strip()
        else:
            goauth_token = None
        incoming['goauth'] = goauth_token

        if "historyid" in incoming:
            if incoming["historyid"] is None or incoming["historyid"] == "**":
                incoming["historyid"] = trans.security.encode_id( trans.history.id )
                print "HISTORY_ID: %s" %  trans.security.encode_id( trans.history.id )
            else:
                print "HISTORY_ID exists: %s" % incoming["historyid"]

        # queue the job run:
        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)
    
        return job, output


class HistoryManagement(DefaultToolAction):
    """This tool basically just works like a globus client."""

    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):
        import socket
        username = trans.user.username
        incoming["userapi"] = trans.user.api_keys[0].key
        incoming["userkey"] = trans.user.api_keys[0].key
        incoming["url"] = "https://%s" % socket.gethostname()
        galaxy_config = trans.app.config.config_dict
        if "historyid" in incoming:
            if incoming["historyid"] is None or incoming["historyid"] == "**":
                incoming["historyid"] = trans.security.encode_id( trans.history.id )

        if 'goauth_token' in incoming and incoming['goauth_token'] not in ['**', '', None]:
            goauth_token = incoming['goauth_token']
        #elif os.path.isfile(os.path.join(galaxy_config['globus_dir'], 'tokens-auth', username)):
        #    with open(os.path.join(galaxy_config['globus_dir'], 'tokens-auth', username), 'r') as f:
        #        goauth_token = f.readline().strip()
        elif os.path.isfile(os.path.join(galaxy_config['globus_dir'], 'tokens-identifiers', username)):
            with open(os.path.join(galaxy_config['globus_dir'], 'tokens-identifiers', username), 'r') as f:
                goauth_token = f.readline().strip()
                print "GOATH_TOKEN: %s" % goauth_token

        elif hasattr(trans, 'environ') and 'X-GLOBUS-TOKEN-AUTH' in trans.environ:
            goauth_token = trans.environ['X-GLOBUS-TOKEN-AUTH']
        else:
            goauth_token = None

        #if goauth_token != None:
        #    goauth_token = check_globus_token(goauth_token, galaxy_config['globus_cred_file'], galaxy_config['globus_dir'], username, 'auth')

        print "GOATH_TOKEN-NOW: %s" % goauth_token
        #if 'identifier_token' in incoming and incoming['identifier_token'] not in ['**', '', None]:
        #    identifier_token = incoming['identifier_token']
        #elif os.path.isfile(os.path.join(galaxy_config['globus_dir'], 'tokens-identifiers', username)):
        #    with open(os.path.join(galaxy_config['globus_dir'], 'tokens-identifiers', username), 'r') as f:
        #        identifier_token = f.readline().strip()
        #elif hasattr(trans, 'environ') and 'X-GLOBUS-TOKEN-AUTH' in trans.environ:
        #    identifier_token = trans.environ['X-GLOBUS-TOKEN-AUTH']
        #else:
        #    identifier_token = None

        #if identifier_token != None:
        #    identifier_token = check_globus_token(goauth_token, galaxy_config['globus_cred_file'], galaxy_config['globus_dir'], username, 'identifier')


        incoming['goauth_token'] = goauth_token
        print "GOATH_TOKEN FINAL: %s" % incoming['goauth_token']
        #incoming['identifier_token'] = identifier_token

        # queue the job run:
        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)

        return job, output

class GlobusOptimizedWorkflows(DefaultToolAction):
    """This tool basically just works like a globus client."""
    GOAUTH_TOKEN_KEY = 'X-GLOBUS-TOKEN'
    REMOTE_USER_KEY = 'X-GLOBUS-USER'

    def execute(self, tool, trans, incoming=None, return_job=False, set_output_hid=True, history=None, job_params=None, rerun_remap_job_id=None, execution_cache=None, dataset_collection_elements=None, completed_job=None, collection_info=None):
        """Before executing the optimized workflow command rename the output 
           to the sample name variable by getting the common name of the inputs 
        """

        #is this paired or single
        read1 = incoming['paired']['input1']
        read2 = None
        if incoming['paired']['sPaired'] == "paired":  
            read2 = incoming['paired']['input2']

        # get the common sample name.
        common = []
        if read2 is not None:
            # get the common string of the two reads
            for count, i in enumerate(read1.name):
                if i == read2.name[count]:
                    common.append(i)
                else:
                    break
        else:
            for i in read1.name:
                if i != ".":
                    common.append(i)
                else:
                    break
        name = "".join(common)

        # get the common sample name.
        incoming['hidden_name'] = name

        incoming['username'] = trans.galaxy_session.user.username
        incoming['goauth'] = trans.environ[self.GOAUTH_TOKEN_KEY]
        galaxy_config = trans.app.config.config_dict
        incoming['from_endpoint'] = galaxy_config['globus_endpoint']

        # queue the job run:
        tool_action_instance = DefaultToolAction()
        job, output = tool_action_instance.execute(tool, trans, incoming=incoming, return_job=return_job, set_output_hid=set_output_hid, history=history, job_params=job_params, rerun_remap_job_id=rerun_remap_job_id, execution_cache=execution_cache, dataset_collection_elements=dataset_collection_elements, completed_job=completed_job, collection_info=collection_info)

        return job, output
