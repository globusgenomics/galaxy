import urllib
import argparse

import requests
import boto
from boto.s3.connection import OrdinaryCallingFormat


def parse_cli():
    description = 'Make a transfrer from Amazon S3'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--session', required=True)
    parser.add_argument('--token', required=True)
    parser.add_argument('--host', required=True)
    parser.add_argument('--ext', required=False)
    parser.add_argument('bucket', nargs=1)
    parser.add_argument('path', nargs=1)
    parser.add_argument('output_path', nargs=1)
    return parser.parse_args()

def fetch_credentials(host, session, token):
    url = urllib.basejoin(host, '/globusonline/s3/browser/credentials')
    response = requests.get(url, params={'token': token},
                            cookies={'galaxysession': session}, verify=False)
    if response.ok:
        return {k: v.encode() for k, v in response.json.items()}
    else:
        raise Exception('Unable to fetch AWS credentials.')

def make_transfer(keyid, sakey, bucket, path, output):
    conn = boto.connect_s3(keyid, sakey, calling_format=OrdinaryCallingFormat())
    bucket = conn.get_bucket(bucket)
    key = bucket.get_key(path)
    with open(output, 'w') as outfile:
        key.get_contents_to_file(outfile)

def get_real_path(output_path, ext):
    if not ext:
        return output_path
    else:
        if output_path.endswith('.data'):
            return output_path[:-4] + ext
        else:
            return output_path

def main():
    args = parse_cli()
    bucket = args.bucket[0]
    path = args.path[0]
    output_path = args.output_path[0]
    cred = fetch_credentials(args.host, args.session, args.token)
    make_transfer(cred['aws_access_key_id'], cred['aws_secret_access_key'],
                  bucket, path, get_real_path(output_path, args.ext))

if __name__ == '__main__':
    main()
    
