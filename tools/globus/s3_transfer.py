#!/usr/bin/python
import boto3
import botocore
import os
import sys
import errno
import ast

from optparse import OptionParser
args=None
parser = OptionParser()

parser.add_option("--transfer-info", dest="transfer_info", help="Detailed Transfer Information")
parser.add_option("--out-file", dest="out_file", help="Output file")

options, args = parser.parse_args(args)

transfer_info = ast.literal_eval(options.transfer_info)

object_size_limit = 100*1000000000 #B, 100GB
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

object_type = transfer_info['object_type']
bucket_name = transfer_info['bucket'].strip().strip('/').strip('\\')
from_path = transfer_info['from_path'].strip().lstrip('/').lstrip('\\')
to_path = transfer_info['to_path'].strip()
if transfer_info['include_subdir'] in [True, False]:
    include_subdir = transfer_info['include_subdir']
else:
    sys.exit("Failed to read include_subdir value.")

if object_type == 'dir':
    prefix = from_path + '/'
elif object_type == 'bucket':
    prefix = None
elif object_type == 'file':
    prefix = from_path
else:
    sys.exit("Object type not supported.")

tags = transfer_info['tags']

# check tags
if tags != []:
    for tag in tags:
        if tag['value'] == '' or tag['key'] == '':
            sys.exit("Empty tags inputs. Please make sure input both Key and Value.")

# Transfer
job_info = []
overall_object_size = 0

def decide_transfer(key):
    if '/' in key:
        if include_subdir:
            return True
        else:
            return False
    else:
        return True

s3_client = boto3.client('s3', aws_access_key_id=aws_access_key_id,
                        aws_secret_access_key=aws_secret_access_key)

def check_tags(key):
    if tags == [] or key.endswith('/'):
        return True
    else:
        try:
            object_tags = s3_client.get_object_tagging(Bucket=bucket_name, Key=key)['TagSet']
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] in ['404', 'NoSuchKey']:
                sys.exit("Getting object tags, the object({0}) does not exist.".format(key))
            elif e.response['Error']['Code'] in ['403', 'InvalidAccessKeyId', 'SignatureDoesNotMatch']:
                sys.exit("Getting object({0}) tags. Forbidden, check your credentials.".format(key))
            else:
                raise
        
        new_tags = []
        for tag in tags:
            tmp = (str(tag['key']), str(tag['value']))
            new_tags.append(tmp)

        new_object_tags = []
        for object_tag in object_tags:
            tmp = (str(object_tag['Key']), str(object_tag['Value']))
            new_object_tags.append(tmp)

        if set(new_tags).issubset(set(new_object_tags)):
            return True
        else:
            return False

        

s3 = boto3.resource('s3', aws_access_key_id=aws_access_key_id,
                        aws_secret_access_key=aws_secret_access_key)

try:
    bucket = s3.Bucket(bucket_name)
    if object_type == 'bucket':
        for obj in bucket.objects.all():
            if decide_transfer(obj.key) and check_tags(obj.key):
                job_info.append([obj.key, os.path.join(to_path, obj.key)])
                overall_object_size += int(obj.size)
    elif object_type == 'dir':
        for obj in bucket.objects.filter(Prefix=prefix):
            if obj.key != prefix:
                shortened_key = obj.key.replace(prefix,'',1)
                if decide_transfer(shortened_key) and check_tags(obj.key):
                    job_info.append([obj.key, os.path.join(to_path, shortened_key)])
                    overall_object_size += int(obj.size)
    elif object_type == 'file':
        if check_tags(from_path):
            job_info.append([from_path, to_path])
            overall_object_size += int(bucket.Object(prefix).content_length)
        else:
            sys.exit("Object({0}) doesn't meet the tags requirement.".format(from_path))
except botocore.exceptions.ClientError as e:
    if e.response['Error']['Code'] in ['404', 'NoSuchKey']:
        sys.exit("The object does not exist.")
    elif e.response['Error']['Code'] in ['403', 'InvalidAccessKeyId', 'SignatureDoesNotMatch']:
        sys.exit("Forbidden, check your credentials")
    else:
        raise

if overall_object_size > object_size_limit:
    sys.exit("The object size exceeds 100GB limit, please contact administrator.")

def check_empty_or_all_dir(ls):
    if ls == []:
        return True
    for i in ls:
        if not i[1].endswith('/'):
            return False
    return True

if check_empty_or_all_dir(job_info):
    sys.exit("Nothing to transfer, check whether there is data under the object to transfer, or check your options such as tags and include subdir.")

for item in job_info:
    if item[1].endswith('/'):
        try:
            os.makedirs(item[1])
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    else:
        if not os.path.exists(os.path.dirname(item[1])):
            try:
                os.makedirs(os.path.dirname(item[1]))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        try:
            bucket.download_file(item[0], item[1])
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] in ['404', 'NoSuchKey']:
                sys.exit("The object ({0}) does not exist.".format(item[0]))
            elif e.response['Error']['Code'] in ['403', 'InvalidAccessKeyId', 'SignatureDoesNotMatch']:
                sys.exit("Forbidden, no pemission to download ({0}), check you credentials.".format(item[0]))
            else:
                raise

def list_files(startpath):
    to_return = ''
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        to_return = to_return + '{}{}/\n'.format(indent, os.path.basename(root))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            to_return = to_return + '{}{}\n'.format(subindent, f)
    return to_return

if object_type in ['dir', 'bucket']:
    to_write = list_files(transfer_info['to_path'])
    with open(options.out_file, 'w') as f:
        f.write(to_write)
