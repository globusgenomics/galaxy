import os
import subprocess
from optparse import OptionParser
args=None
parser = OptionParser()

parser.add_option("--user", dest="user", help="user email")
parser.add_option("--globus-cred-file", dest="globus_cred_file", help="globus cred file path")
parser.add_option("--from-endpoint", dest="from_endpoint", help="from endpoint")
parser.add_option("--to-endpoint", dest="to_endpoint", help="to endpoint")
parser.add_option("--from-dataset", dest="from_dataset", help="from dataset")
parser.add_option("--to-path", dest="to_path", help="to_path")
parser.add_option("--extra-source-path", dest="extra_source_path", help="extra source path", default=None)

options, args = parser.parse_args(args)

#print options

username = options.user[0:options.user.find('@')].strip()

dataset_path = options.from_dataset
# Check whether it is a Directory object thansferred through Globus Genomics
dir_link = dataset_path[0:-4] + "_files"
if os.path.islink(dir_link):
    dir_realpath = os.path.realpath(dir_link)
    dataset_path = dir_realpath
elif os.path.isdir(dir_link):
    dataset_path = dir_link

dataset_path = dataset_path.rstrip('/').rstrip('\\')

to_path = options.to_path.strip().rstrip('/').rstrip('\\').replace(' ', '_').replace(':\\', '/').replace('\\', '/')

transfer_item = "{0}::to::{1}".format(dataset_path, to_path)

access_token_file = os.path.join(os.path.join(os.path.dirname(options.globus_cred_file), "tokens"), username)
with open(access_token_file, 'r') as f:
    access_token = f.readlines()[0].strip()

script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "globus_transfer_v2.py")

# construct the command
command = "python {0} --username {1} --globus-cred-file {2} --from-endpoint {3} --to-endpoint {4} -i {5} --goauth-token {6} --transfer-direction out".format(script_path, username, options.globus_cred_file, options.from_endpoint, options.to_endpoint, transfer_item, access_token)
if options.extra_source_path != None:
    command = command + " --extra-source-path {0}".format(options.extra_source_path)

#print command

subprocess.call(command, shell=True)
