#!/usr/bin/python
import boto3
import botocore
import os
import sys
import ast

from optparse import OptionParser
args=None
parser = OptionParser()

parser.add_option("--transfer-info", dest="transfer_info", help="Detailed Transfer Information")

options, args = parser.parse_args(args)

transfer_info = ast.literal_eval(options.transfer_info)

cred_file_dir = '/home/galaxy/.globusgenomics/tmp'
cred_file_loc = os.path.join(cred_file_dir, transfer_info['cred_file']) 

with open(cred_file_loc, 'r') as f:
    lines = f.readlines()
    try:
        aws_access_key_id = lines[0].strip()
    except:
        aws_access_key_id = None
    try:
        aws_secret_access_key = lines[1].strip()
    except:
        aws_secret_access_key = None
os.remove(cred_file_loc)

if aws_access_key_id in [None, '']:
    sys.exit("Could not get aws_access_key_id.")
if aws_secret_access_key in [None, '']:
    sys.exit("Could not get aws_secret_access_key.")

bucket_name = transfer_info['bucket'].strip().strip('/').strip('\\')
if transfer_info['sse'] in [True, False]:
    server_side_encryption = transfer_info['sse']
else:
    sys.exit("Failed to read SSE value.")
tags = transfer_info['tags']
# check tags
if tags != []:
    tags_to_put = []
    for tag in tags:
        if tag['value'] == '' or tag['key'] == '':
            sys.exit("Empty tags inputs. Please make sure input both Key and Value.")
        tmp = {'Key': str(tag['key']), 'Value': str(tag['value'])}
        tags_to_put.append(tmp)

# Transfer
s3 = boto3.resource('s3', aws_access_key_id=aws_access_key_id,
                        aws_secret_access_key=aws_secret_access_key)
bucket = s3.Bucket(bucket_name)

def upload_file(from_loc, to_loc):
    try:
        if server_side_encryption:
            bucket.upload_file(from_loc, to_loc, ExtraArgs={'ServerSideEncryption': 'AES256'})
        else:
            bucket.upload_file(from_loc, to_loc)
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] in ['404', 'NoSuchKey']:
            sys.exit("The object does not exist.")
        elif e.response['Error']['Code'] in ['403', 'InvalidAccessKeyId', 'SignatureDoesNotMatch']:
            sys.exit("Forbidden, no pemission to upload to ({0}), please check you credentials.".format(to_loc))
        else:
            raise
    if tags != []:
        # tag
        s3_client = boto3.client('s3', aws_access_key_id=aws_access_key_id,
                                    aws_secret_access_key=aws_secret_access_key)
        try:
            s3_client.put_object_tagging(Bucket=bucket_name, Key=to_loc, Tagging={'TagSet':tags_to_put})
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] in ['404', 'NoSuchKey']:
                sys.exit("Problem when tagging, the object({0}) does not exist.".format(to_loc))
            elif e.response['Error']['Code'] in ['403', 'InvalidAccessKeyId', 'SignatureDoesNotMatch']:
                sys.exit("Problem when tagging, Forbidden to tag object({0}), check your credentials.".format(to_loc))
            else:
                raise

for job in transfer_info["jobs"]:
    object_type = job['object_type']
    from_path = job['from_path'].strip().rstrip('/').rstrip('\\')
    to_path = job['to_path'].strip().strip('/').strip('\\')
    if job['rename'] != '':
        rename = job['rename'].strip().strip('/').strip('\\')
    else:
        rename = None

    # Check if the to_path is a existing file
    try:
        objs = list(bucket.objects.filter(Prefix=to_path))
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] in ['404', 'NoSuchKey']:
            sys.exit("The bucket does not exist.")
        elif e.response['Error']['Code'] in ['403', 'InvalidAccessKeyId', 'SignatureDoesNotMatch']:
            sys.exit("Forbidden, check you credentials")
        else:
            raise
    if len(objs) > 0 and objs[0].key == to_path:
        sys.exit("{0} is a existing file in your bucket.".format(to_path))

    if object_type == 'dir':
        if os.path.isdir(from_path):
            for root, dirs, files in os.walk(from_path):
                if rename != None:
                    dir_name = os.path.join(rename, root[len(from_path):].strip('/'))
                else:
                    dir_name = os.path.join(os.path.basename(from_path), root[len(from_path):].strip('/'))
                for item in files:
                    upload_file(os.path.join(root, item),
                                os.path.join(to_path, dir_name, item))
        else:
            sys.exit("{0} does not exist.".format(from_path))

    elif object_type == 'file':
        if os.path.isfile(from_path):
            if rename != None:
                file_name = rename
            else:
                file_name = os.path.basename(from_path)
            upload_file(from_path, os.path.join(to_path, file_name))
        else:
            sys.exit("{0} does not exist.".format(from_path))
    else:
        sys.exit("Object type not supported.")

