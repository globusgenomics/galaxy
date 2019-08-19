import logging
import time
import socket
from Cookie import BaseCookie
from cgi import parse_qs
import globus_sdk
from cStringIO import StringIO
import json
import ConfigParser
import re
import os
import ast
import requests
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()

log = logging.getLogger(__name__)

AUTH_PAGE = """
<html lang="en">
<head>
<title>Globus Genomics : Login</title>
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous"> 
<script src=https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>
<style type="text/css">
body {
    min-width: 500px;
    text-align: center;
}
.navbar-bg {
    background-color: #1D5BA9;
    padding-top: 20;
    font-size: 19;
    font-color: white;
    height: 80;
}
.brand-name {
    color: white;
}

.logo {
    height: 100%%;
}

.gg-panel {
    width: 500;
    margin: auto;
}
</style>
</head>
<body>
    <nav class="navbar navbar-default" style="border: 0px;">
            <div class="container-fluid navbar-bg">
            <div class="navbar-header brand-name">
                <img class="logo" height="280%%" src="https://www.globus.org/sites/default/files/logo_globus_genomics_0.png" alt="Home" /> globus <b>genomics</b>
            </div>
            </div>
     </nav>
    <div class="gg-panel">

            <div class="panel panel-primary">
                <div class="panel-body">
                <h4>Log in using your Globus account:</h4>
                <hr>
                <a href="%(url)s" class="btn btn-primary btn-lg active" role="button">Authenticate using Globus</a>
                <!-- <div class="alert alert-danger" role="alert"><h6>%(message)s</h6></div> -->
                <h6><p class="bg-danger">%(message)s</p></h6>
                <hr>
                <h4>
                   <small>Note: Globus Genomics uses Globus groups to manage access to the instance. To get access to your instance, please
                         contact your administrator or reach us at <a href="mailto:support@globus.org?subject=Access to Globus Genomics instance">support@globus.org</a>
                   </small>
                </h4>
                <hr>
                <h4>
                <small>Don't have a Globus account? Create an account
                        using <a href="https://www.globusid.org/create">Globus ID</a>
                </small>
                </h4>
            </div>
           </div>
    </div>
</body>
</html>

"""

class SessionStorage():
    def __init__(self, paste_session):
        self.paste_session = paste_session

    def _get_in_memory_storage(self):
        try:
            ses_storage = self.paste_session['globus.storage']
        except KeyError:
            ses_storage = self.paste_session['globus.storage'] = {}
        finally:
            return ses_storage
        
    def store(self, key, value):
        ses_storage = self._get_in_memory_storage()
        ses_storage[key] = value

    def fetch(self, key):
        ses_storage = self._get_in_memory_storage()
        return ses_storage[key]

    def pop(self, key):
        ses_storage = self._get_in_memory_storage()
        return ses_storage.pop(key)


class UserAuthentication(object):
    LOGOUT = 'LOGOUT'
    LOGIN = 'LOGIN'
    AUTHENTICATED = 'AUTHENTICATED'
    NOT_AUTHORIZED = 'NOT_AUTHORIZED'
    NOT_ALLOWED = 'NOT_ALLOWED'
    IMPERSONATE = 'IMPERSONATE' 
    API_CALL = 'API_CALL'
    DISPLAY_SERVER = 'DISPLAY_SERVER'
    # Redundant names just to avoid any string typos.

    def __init__(self, app, display_servers=None):
        log.debug('init')
        self.display_servers = display_servers or []
        self.app = app

        self.handlers = {
            self.LOGOUT: Logout(app),
            self.LOGIN:  Login(app),
            self.AUTHENTICATED: Authenticate(app),
            self.NOT_AUTHORIZED: NotAuthorized(app),
            self.NOT_ALLOWED: NotAllowed(app),
            self.IMPERSONATE: Impersonate(app),
            self.API_CALL: APICall(app),
            self.DISPLAY_SERVER: DisplayServer(app)
        }

    def __call__(self, environ, start_response):
        session = self._session(environ)
        self.paste_session = environ['paste.session.factory']()
        qs = parse_qs(environ['QUERY_STRING'])
        name, params = self._detect_type_of_request(session, environ, qs)
        handler = self.handlers[name]
        return handler(session, start_response, environ, *params)

    def _req_cookie(self, environ):
        cookie = BaseCookie()
        cookie.load(environ.get('HTTP_COOKIE', ''))
        print ""
        return cookie

    def _session(self, environ):
        cookie = self._req_cookie(environ)
        if cookie:
            try:
                session = cookie['galaxysession']
            except KeyError:
                return None
            else:
                return session.value

    def _detect_type_of_request(self, session, env, qs):
        """Always return a type of request and the arguments
        required for the handler.
        
        :returns: handler_name, (params)
        """
        # The "globusonline" part of this validation is kind of a hack
        # just to allow the "display_servers" part that the remote_user
        # middleware support and still work with the custom tools that
        # needs to make request with the user.
        path_info = env['PATH_INFO']
        #print "PATH_INFO"
        #print path_info
        referer = env.get('HTTP_REFERER', '')
        #print "Display Servers"
        #print self.display_servers
        # Validate that the request comes from a display server
        #if (self.display_servers and self._is_a_valid_display_server(env)):
        if (self.display_servers and 
                not path_info.startswith('/globus.org/') and
                self._is_a_valid_display_server(env)):
            print "returning display server"
            return self.DISPLAY_SERVER, ()
        #if (path_info.startswith('/display_application/')):
        #    print "this should be true now"
        #    return self.__normal_flow_on_dtype_of_request(session)
            #return self.DISPLAY_SERVER, () 
        # Handle OAuth
        if path_info == '/':
            result = self._handle_oauth_flow(session, qs)
            if result is not None: # this is an OAuth flow.
                #print "PRINTING RESULT"
                #print result
                return result
        if path_info.startswith('/user/logout'):
            return self.LOGOUT, ()
        # Handle the api calls
        if path_info.startswith('/api/'):
            if 'HTTP_COOKIE' in env:
            #if referer.endswith('/history') or referer.endswith('/genomes') or referer.startswith('/api/datatypes') or referer.startswith('/api/tools'): # the request comes from the browser.
                return self.__normal_flow_on_dtype_of_request(session)
            else:
                return self.API_CALL, ()
        # Handle the impersonate case.
        if (env['REQUEST_METHOD'] == 'POST' and 
              referer.endswith('/admin/impersonate')):
                return self.IMPERSONATE, ()
        # If none of the special cases are used then use the
        # default handlers AUTHENTICATED or NOT_AUTHENTICATED.
        return self.__normal_flow_on_dtype_of_request(session)

    def _is_a_valid_display_server(self, env):
        try:
            host = socket.gethostbyaddr(env['REMOTE_ADDR'])[0]
            #print "HOST in display_server"
            #print host
        except (socket.error, socket.herror,
                socket.gaierror, socket.timeout, KeyError):
            return False
        else:
            return host in self.display_servers

    def _handle_oauth_flow(self, session, qs):
        if 'action' in qs:
            action = qs['action'][0]
            if action == 'logout' and \
                    session in self.paste_session:
                return self.LOGOUT, ()
            elif action == 'user_not_allowed':
                user_not_allowed = None
                if 'user' in qs:
                    user_not_allowed = qs['user'][0]
                    return self.NOT_ALLOWED, (user_not_allowed,)
        if 'code' in qs:
            code = qs['code'][0]
            return self.LOGIN, (code,)
        return None

    def __normal_flow_on_dtype_of_request(self, session):
        #print "I am here in __normal_flow_on_dtype_of_request"
        user_n_token = self._user_in_session(session)
        if user_n_token is not None:
            return self.AUTHENTICATED, user_n_token
        return self.NOT_AUTHORIZED, ()

    def _user_in_session(self, session):
        if session is not None and \
               session in self.paste_session:
            user, tokens = self.paste_session[session]
            return user, tokens
        return None


class ActionHandler(object):
    required_args = set([])

    def __init__(self, app, **kwargs):
        self.app = app
        self.galaxy_config = self.app.app.application.api_controllers["configuration"].app.config.config_dict
        if self.required_args - set(kwargs):
            raise Exception("Missing required arguments %s" % self.required_args)
        else:
            for name, value in kwargs.items():
                setattr(self, name, value)

    def get_globus_auth_url(self):
        self.load_auth_client()
        auth_url = self.client.oauth2_get_authorize_url()
        return auth_url

    def load_auth_client(self):
        Config = ConfigParser.ConfigParser()
        Config.read(self.galaxy_config['globus_cred_file'])
        client_id = Config.get('globus', 'client_id')
        client_secret = Config.get('globus', 'client_secret')
        redirect_uri = Config.get('globus', 'redirect_uri')

        self.client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
        self.client.oauth2_start_flow(redirect_uri, refresh_tokens=True)

class BypassHandler(object):
    """Just a transparent class to pass the request to Galaxy
    """

    def __init__(self, app):
        self.app = app

    def __call__(self, session, start_response, environ):
        return self.app(environ, start_response)

class APICall(BypassHandler):
    """
    The API calls are authenticated by the api keys, just pass the request.
    """

class DisplayServer(BypassHandler):
    """
    If the request is a trusted display_server this will pass
    the request as it came to the the remote_user middleware.
    """

class Authenticate(ActionHandler):

    def __call__(self, session, start_response, environ, user, atoken):
        # check whether the token is valid
        self.load_auth_client()
        if atoken and not self.client.oauth2_validate_token(atoken)['active']:
            # if not valid, then log out
            try:
                paste_session = environ['paste.session.factory']()
                del paste_session[session]
            except:
                pass
            try:
                del paste_session['globus.storage']
            except KeyError:
                pass
            start_response( '303 See other', [('Content-type', 'text/html'),
                                              ('Location', 'https://auth.globus.org/v2/web/logout')])
            return ''

        #  HTTP_REMOTE_USER will be modified by galaxy.
        environ[ self.app.remote_user_header ] = user
        environ['X-GLOBUS-USER'] = user
        environ['X-GLOBUS-TOKEN'] = atoken
        self.record_token(user, atoken)

        environ['globus.storage'] = SessionStorage(environ['paste.session.factory']())
        return self.app(environ, start_response)

    def record_token(self, user, atoken):
        if '@' in user:
            file_name = user[0:user.find('@')]
        else:
            file_name = user
        
        dir_name =  os.path.join(self.galaxy_config['globus_dir'], 'tokens')
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name, 0700)

        record_file = os.path.join(dir_name, file_name)

        with open(record_file, 'w') as write_file:
            write_file.write(atoken)


class Login(ActionHandler):
    def __call__(self, session, start_response, environ, code):
        self.load_auth_client()
        tokens = self.client.oauth2_exchange_code_for_tokens(code)
        transfer_token = tokens.by_resource_server['transfer.api.globus.org']
        transfer_access_token = transfer_token['access_token']
        transfer_token_info = self.client.oauth2_token_introspect(transfer_access_token)
        
        username = transfer_token_info['username']
        userid = transfer_token_info['sub']

        if self.galaxy_config['globus_use_group'] and not self._is_user_in_group(username):
            return self._user_not_in_group(username, start_response)
        environ[ self.app.remote_user_header ] = username

        def custom_start_response(status, headers, exc_info=None):
            """This is closure to bind Galaxy Session with Globus User.
            Look for the Set-Cookie header in the next middleware.
            start_response( '303 See other', [('Content-type', 'text/html'),
            ('Location', 'https://auth.globus.org/v2/web/logout')])
            """
            #print "in custom_start_response"
            if session is None:
                for name, value in headers:
                    if name == 'Set-Cookie':
                        cookie = BaseCookie()
                        cookie.load(value)
                        if 'galaxysession' in cookie:
                            sess = cookie['galaxysession'].value
                            self._create_user(sess, environ, username, userid, transfer_access_token)
                        break
            else:
                self._create_user(session, environ, username, userid, transfer_access_token)
            headers  += [('Content-type', 'text/html'),
                         ('Location', '/')]
            #print "RETURNING HERE"
            return start_response('303 See Other', headers, exc_info)
        return self.app(environ, custom_start_response)

    def _create_user(self, _session, _environ, _user, _userid, _tokens):
        # share endpoint
        self._share_endpoint(_user, _userid)

        paste_session = _environ['paste.session.factory']()
        paste_session[_session] = (_user, _tokens)

    def _share_endpoint(self, username, userid):
        """
        For shared endpoint: Share the endpoint with a user
        - check whether a user's dir exits in the shared dir (/scratch/shared), if not then share
        - determine whether the endpoint in use is shared endpoint
        - check whether the user is in the rule list
        - create user's dir under shared dir and share the endpoint with the user

        Requirement: globus_sdk, CRED_FILE
        """
        if '@' in username:
            username = username[0:username.find('@')]
        # Config
        SHARED_DIR = self.galaxy_config["ftp_upload_dir"].strip()
        GLOBUS_DIR = self.galaxy_config["globus_dir"].strip()
        ENDPOINT_NAME = self.galaxy_config["globus_endpoint"].strip().replace('"', '').strip()
        GLOBUS_USER_ID = userid
        GLOBUS_USER_NAME = username

        # Check if the endpoint is shared with the user
        shared_user_list_file = os.path.join(GLOBUS_DIR, "shared_endpoint_shared_user_list")
        if os.path.isfile(shared_user_list_file):
            with open(shared_user_list_file) as f:
                users = [line.strip() for line in f]
            if GLOBUS_USER_NAME in users:
                msg = "{0} is shared with {1}".format(ENDPOINT_NAME, GLOBUS_USER_NAME)
                log.info(msg)
                return

        path = os.path.join(SHARED_DIR, GLOBUS_USER_NAME)

        CRED_FILE = self.galaxy_config["globus_transfer_cred_file"]

        client_id = None
        transfer_rt = None
        transfer_at = None
        expires_at_s = None
        with open(CRED_FILE) as f:
            for line in f:
                if line.startswith("client_id:"):
                    client_id = re.sub("client_id:","",line,count=1).strip()
                if line.startswith("refresh_token:"):
                    transfer_rt = re.sub("refresh_token:","",line,count=1).strip()
                if line.startswith("access_token:"):
                    transfer_at = re.sub("access_token:","",line,count=1).strip()
                if line.startswith("expires_at_seconds:"):
                    expires_at_s = int(re.sub("expires_at_seconds:","",line,count=1).strip())

        if None in [client_id, transfer_rt, transfer_at, expires_at_s]:
            msg = "Cannot get the creds needed for globus_sdk"
            raise Exception(msg)

        client = globus_sdk.NativeAppAuthClient(client_id)

        authorizer = globus_sdk.RefreshTokenAuthorizer(
                        transfer_rt, client, access_token=transfer_at, expires_at=expires_at_s)

        tc = globus_sdk.TransferClient(authorizer=authorizer)

        endpoint = None
        # Search for the shared endpoint
        for ep in tc.endpoint_search(filter_fulltext=ENDPOINT_NAME, filter_scope='shared-by-me'):
            if ep["display_name"] == ENDPOINT_NAME:
                endpoint = ep
                break
        # make sure it is a shared endpoint:
        if endpoint == None or endpoint["host_endpoint_id"] == None:
            msg = "Cannot find shared endpoint {0} or it is not a shared endpoint".format(ENDPOINT_NAME)
            log.info(msg)
        else:
            # Share with a user
            endpoint_id = endpoint["id"]
            acl_list = tc.endpoint_acl_list(endpoint_id)
            # check if the user exists in the rules:
            proceed = True
            for item in acl_list["DATA"]:
                if item["principal"] == GLOBUS_USER_ID:
                    msg = "user {0} already exist in the access rule".format(GLOBUS_USER_NAME)
                    log.info(msg)
                    proceed = False
                    break
            if proceed:
                rule_data = {
                               "DATA_TYPE": "access",
                               "principal_type": "identity",
                               "principal": GLOBUS_USER_ID,
                               "path": "/shared/{0}/".format(GLOBUS_USER_NAME),
                               "permissions": "rw",
                             }
                tc.add_endpoint_acl_rule(endpoint_id, rule_data)
                rule_data = {
                               "DATA_TYPE": "access",
                               "principal_type": "identity",
                               "principal": GLOBUS_USER_ID,
                               "path": "/galaxy/files/",
                               "permissions": "r",
                             }
                tc.add_endpoint_acl_rule(endpoint_id, rule_data)

        # create user's dir
        if not os.path.isdir(path):
            os.mkdir(path, 0755)
        # record the user has been shated with
        with open(shared_user_list_file, "a") as f:
            to_write = "{0}\n".format(GLOBUS_USER_NAME)
            f.write(to_write)  

    def _user_not_in_group(self, user, start_response):
        tmp = '/?action=user_not_allowed&user={0}'.format(user)
        headers = [('Content-type', 'text/html'),
                   ('Location', tmp)]
        start_response('303 See other', headers)
        return ''

    def _is_user_in_group(self, user):
        if '@' in user:
            user = user[0:user.find('@')]
        users = []

        tmp_user_list = '/home/galaxy/.globusgenomics/globus_user_granted'
        if os.path.isfile(tmp_user_list):
            with open(tmp_user_list) as f:
                users = [line.strip() for line in f]
        else:
            globus_group_id = self.galaxy_config['globus_group_id']
            Config = ConfigParser.ConfigParser()
            Config.read(self.galaxy_config['globus_cred_file'])
            refresh_token = Config.get('globus_group', 'refresh_token')
            client_id = Config.get('globus_group', 'client_id')
            client_secret = Config.get('globus_group', 'client_secret')
            client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
            token_response = client.oauth2_refresh_token(refresh_token)
            globus_data = token_response.by_resource_server['nexus.api.globus.org']
            access_token = globus_data['access_token']

            url = 'https://nexus.api.globusonline.org/groups/{0}/members'.format(globus_group_id)
            headers = {'authorization': 'Bearer {0}'.format(access_token)}
            r = requests.get(url, headers=headers)
            response = r.json()
            users = []
            for i in response['members']:
                if i['status'] == 'active':
                    users.append(i['username'])
                    
        if user in users:
            return True
        else:
            return False

class Logout(ActionHandler):

    def __call__(self, session, start_response, environ, *args):
        try:
            paste_session = environ['paste.session.factory']()
            del paste_session[session]
        except:
            pass
        try:
            del paste_session['globus.storage']
        except KeyError:
            pass
        start_response( '303 See other', [('Content-type', 'text/html'),
                                          ('Location', 'https://auth.globus.org/v2/web/logout')])
        return ''


class NotAuthorized(ActionHandler):
    def __call__(self, session, start_response, environ):
        #print "I am in NotAuthorized"
        start_response( '303 See Other', [('Content-type', 'text/html'), 
            ('Location', self.get_globus_auth_url())] )
        
        return ''


class NotAllowed(ActionHandler):
    def __call__(self, session, start_response, environ, user):
        """Handle the case when the user is not allowed in this instance
        for the necessary group permissions.
        """
        if user is None:
            message = 'The user is not allowed'
        else:
            message = ('The user <strong>%s</strong> is not '
                       'allowed in this instance.' % user)
        start_response( '403 Forbidden', [('Content-type', 'text/html')])
        return AUTH_PAGE % {'url': self.get_globus_auth_url(),
                            'message': message}


class Impersonate(ActionHandler):
    def __call__(self, session, start_response, environ):
        impuser = self._catch_user_for_impersonate_request(session, environ)

        def custom_start_response(status, headers, exc_info=None):
            self._impersonate_session(session, environ, impuser)
            return start_response('200 OK', headers, exc_info)

        self._set_current_user(session, environ)
        return self.app(environ, custom_start_response)

    def _set_current_user(self, session, environ):
        paste_session = environ['paste.session.factory']()
        user, tokens = paste_session[session]
        environ[ self.app.remote_user_header ] = user
        environ['X-GLOBUS-USER'] = user
        environ['X-GLOBUS-TOKEN'] = tokens
        environ['globus.storage'] = SessionStorage(environ['paste.session.factory']())

    def _impersonate_session(self, session, environ, email):
        user_token = ({'username': email.split('@')[0]}, {'access': ''})
        paste_session = environ['paste.session.factory']()
        paste_session[session] = user_token

    def _catch_user_for_impersonate_request(self, session, env):
        body = env['wsgi.input'].read()
        params = parse_qs(body)
        env['wsgi.input'] = StringIO(body)
        return params['email'][0]

