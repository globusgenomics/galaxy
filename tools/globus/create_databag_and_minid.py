import sys, os, shutil
import optparse
import requests
import json

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '-o', '--output', dest='output_path', action='store', type='string', help='output file' )
    parser.add_option( '-p', '--bag_path', dest='bag_path', action='store', type='string', help='physical location where bag will be stored' )
    parser.add_option( '-d', '--dataset', dest='datasets', action='append', type="string", nargs=4, help='"file_path" "extra_files_path" "dataset_name" "datatype"' )
    parser.add_option( '-t', '--token', dest='token', action='store', type='string', help='Globus auth token' )
    parser.add_option( '-u', '--minid-user', dest='minid_user', action='store', type='string', help='Minid user' )
    parser.add_option( '-e', '--minid-email', dest='minid_email', action='store', type='string', help='Minid email' )
    parser.add_option( '-m', '--minid-title', dest='minid_title', action='store', type='string', help='Minid title' )
    (options, args) = parser.parse_args()

    try:
        os.makedirs(options.bag_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    minid_post = {
                   "minid_user": options.minid_user,
                   "minid_email": options.minid_email,
                   "access_token": options.token,
                   "minid_title": options.minid_title,
                 }

    endpoint_id = "d1c150fa-b43d-11e7-b0a7-22000a92523b"
    if options.datasets:
        for ( dataset, extra_files_path, name, ext ) in options.datasets:
            target_filepath = "%s/%s.%s" % (options.bag_path, name.replace(" ", "_").replace(":", "_").replace("(", "_").replace(")", "_"), ext.lower())
            shutil.copy2(dataset, target_filepath)

            # copy extra_files_path directory if it exists
            if os.path.isdir(extra_files_path):
                dst = "%s/%s.directory" % (options.bag_path, name.replace(" ", "_").replace(":", "_"))
                shutil.copytree(extra_files_path, dst)

    # build the BDBag and minid
    list_of_files = []
    for (dirpath, dirnames, filenames) in os.walk(options.bag_path):
        for filename in filenames:
            file_dict = {}
            file_dict['url'] = "globus://%s:%s" % (endpoint_id, os.sep.join([dirpath, filename]))
            file_dict['filename'] = filename
            list_of_files.append(file_dict)

    minid_post["remote_files_manifest"] = list_of_files

    # POST minid_post to create BDbag and minid
    url = "https://portal.sc17.nick.globuscs.info/concierge/api/bags/"
    headers = {'content-type': 'application/json'}
    r = requests.post(url, data=json.dumps(minid_post), headers=headers)

    minid_id = None
    if r.status_code == 201:
        minid_id = r.json()['minid_id']
        fh = open(options.output_path, "w")
        fh.write("%s\n" % minid_id)
        fh.write("%s\n" % r.json())
        fh.close()

if __name__=="__main__": __main__()
