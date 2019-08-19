import os
import functools
import urllib
from argparse import ArgumentParser

import requests



class CatalogAPI(object):
    #AUTH_SERVER = 'graph.api.test.globuscs.info'
    AUTH_SERVER = 'graph.api.globusonline.org'
    AUTH_URL = '/goauth/token?grant_type=client_credentials'
    DS_SERVER  = 'datasets.globuscs.info'
    DS_BASE_URI = '/tagfiler/catalog'
    PROTOCOL = 'https'


    def __init__(self, gouser, gopasswd='', gotoken=None, catalog_id=1,
                 auth_server=None, ds_server=None, verify_certs=False):
        self.user = gouser
        self.passwd = gopasswd
        self.catalog_id = catalog_id
        self.gotoken = gotoken
        self.verify_certs = verify_certs
        if auth_server is None:
            self.auth_server = self.AUTH_SERVER
        if ds_server is None:
            self.ds_server = self.DS_SERVER
        self.auth_server_url = '{0}://{1}{2}'.format(
            self.PROTOCOL,self.auth_server, self.AUTH_URL)
        self.ds_server_url = '{0}://{1}{2}/{3}'.format(
            self.PROTOCOL, self.ds_server, self.DS_BASE_URI,
            self.catalog_id)

    def _login(self):
        auth = requests.auth.HTTPBasicAuth(self.user, self.passwd)
        response = requests.get(self.auth_server_url,
                                auth=auth, verify=self.verify_certs)
        if response.ok:
            self.gotoken= response.json['access_token']
            return  True
        else:
            return False

    @property
    def auth_request(self):
        if self.gotoken is None:
            raise Exception('Unable to get authenticated request, '
                            'the instance does not have an access code.')
        else:
            return functools.partial(requests.get,  verify=self.verify_certs, 
                                     headers={
                                         'Authorization': 
                                         'Globus-Goauthtoken %s' % self.gotoken})

    def auth_get(self, url_part, is_part=True):
        if is_part:
            url = self.ds_server_url + url_part
        else:
            url = url_part
        return self.auth_request(url)

    def login(self, validate=True):
        if validate:
            if self.gotoken is None:
                return self._login()
            else:
                return True
        else:
            return self._login()

    def query(self,  attributes):
        resp = self.auth_get('/subject/{0}'.format(attributes))
        if resp.ok:
            return resp.json
        else:
            raise Exception('Invalid request %s' % resp.text)

    def files_of_dataset(self, id_):
        return self.query('dataset_reference=%s(http_name;http_path;http_size;id)' % id_)



class BasespaceAPI(object):

    def __init__(self, access_token):
        self.atoken = access_token


    def get_file(self, location, outpath):
        response = requests.get(location + '/content',
                                prefetch=False,
                                headers={'x-access-token': self.atoken})
        with open(outpath, 'w') as output:
            for chunk in response.iter_content(2048):
                if not chunk:
                    break
                output.write(chunk)


def parse_cli():
    description = 'Fetch datasets from the Globus Online Catalog'
    parser = ArgumentParser(description=description)
    parser.add_argument('--session', required=True)
    parser.add_argument('--token', required=True)
    parser.add_argument('--host', required=True)
    parser.add_argument('--files', default='')
    parser.add_argument('catalog', nargs=1)
    parser.add_argument('dataset_ref', nargs=1)
    parser.add_argument('output_primary', nargs=1)
    parser.add_argument('output_id', nargs=1)
    parser.add_argument('output_dir', nargs=1)
    return parser.parse_args()



def fetch_credentials(host, session, token):
    url = urllib.basejoin(host, '/globusonline/catalog/browser/credentials')
    response = requests.get(url, params={'token': token},
                            cookies={'galaxysession': session})
    if response.ok:
        return {k: v.encode() for k, v in response.json.items()}
    else:
        raise Exception('Unable to fetch the GO Catalog credentials.')





def make_transfer(gouser, gotoken, bsaccess, catalog, dataset_id,
                  outlog, out_id, out_dir, selected_files=()):
    capi = CatalogAPI(gouser, gotoken=gotoken, catalog_id=int(catalog))
    bapi = BasespaceAPI(bsaccess)
    files = capi.files_of_dataset(dataset_id)
    if selected_files:
        files = [f for f in files if f['id'] in selected_files]
    if not files:
        outlog('There is no file matching that criteria.\n')
        raise Exception('There is no file matching that criteria.')
    outlog('Fetching %s %s:\n' % (len(files), 
                                  'files' if  len(files) > 1 else 'file'))
    for f in files:
        outlog('\t- %s [%s bytes]\n' % (f['http_name'], f['http_size']))
        name = f['http_name'].replace('_', '-')
        location = f['http_path']
        out_name = "%s_%s_%s_%s_%s" % \
                    ('primary', out_id, name, 'visible', 'txt')
        bapi.get_file(location, os.path.join(out_dir, out_name))
    outlog('\n')


def main():
    args = parse_cli()
    catalog = args.catalog[0]
    dataset_ref = args.dataset_ref[0]
    output_log = args.output_primary[0]
    output_id = args.output_id[0]
    output_dir = args.output_dir[0]
    cred = fetch_credentials(args.host, args.session, args.token)
    if args.files.strip():
        try:
            files = map(int, args.files.strip().split(','))
        except ValueError:
            args.error('Invalid format of file id.')
    else:
        files = None
    # we are going to ignore catalog for the demo purposes.
    log_msg = ('Globus Online Catalog Tool Information\n\n' 
               'User: %s\n'
               'Catalog: %s\n'
               'Dataset: %s\n'
               'Selected files: %s\n') % \
               (cred['gouser'], 'BIRN Catalog', dataset_ref, files)
    with open(output_log, 'w') as outlog:
        outlog.write(log_msg)
        outlog.write(('-' * 80) + '\n')
        make_transfer(cred['gouser'], cred['gotoken'], cred['bsaccess'],
                      catalog, dataset_ref, outlog.write, output_id,
                      output_dir, files)


if __name__ == '__main__':
    main()
