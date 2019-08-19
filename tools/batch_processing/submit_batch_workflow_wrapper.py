#!/usr/bin/python
#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse, os, shutil, subprocess, sys, tempfile, psycopg2, socket, re
#from galaxy.web import url_for
import requests
requests.packages.urllib3.disable_warnings()

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()


parser = optparse.OptionParser()
parser.add_option( '--user', dest='userid', help='The user id from the postgresql table associated with the user' )
parser.add_option( '--input', dest='table_file', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--out', dest='log', help='Log file' )
parser.add_option( '--username', dest='username', help='The user name from the postgresql table associated with the user' )
parser.add_option( '--goauth-token', dest='goauth', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--rootdir', dest='rootdir', help='The root directory for galaxy installation' )
parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
(options, args) = parser.parse_args()

# get the database, user, host, password, url
#print os.getcwd()
universe_file = options.rootdir + "/universe_wsgi.ini"
ufile = open( universe_file, "r" )
for line in ufile:
    if "database_connection = " in line:
        matches = re.search(r"database_connection = postgres:\/\/\/(.*)\?user=(.*)\&password=(.*)\&host=(.*)", line)
        database = matches.group(1)
        user = matches.group(2)
        host = matches.group(4)
        password = matches.group(3)

    #elif "display_servers = " in line: 
        #matches = re.search(r"display_servers = (.*)", line)
        #url_host = matches.group(1)
    elif "globus_endpoint = " in line:
        matches = re.search(r"globus_endpoint = (.*)#(.*)", line)
        url_host = matches.group(2)
        
ufile.close()
#print "URL_UNIVERSE: %s" % url

# comment next line if socket.gethostname is not giving correct hostname
url_host = socket.gethostname()
#print "URL_HOST: %s" % url_host
url = "https://%s" % url_host

con = None
try:
    con = psycopg2.connect(database=database, user=user, host=host, password=password) 
    cur = con.cursor()

    # get the user key
    #key_query_sql = "SELECT  api_keys.key FROM api_keys INNER JOIN galaxy_user ON galaxy_user.id=api_keys.user_id where galaxy_user.username='%s'" % (userid)
    #key_query_sql = "SELECT  api_keys.key FROM api_keys WHERE api_keys.user_id='%s'" % (options.userid)
    key_query_sql = "SELECT  api_keys.key FROM api_keys WHERE api_keys.user_id='%s' ORDER BY api_keys.create_time DESC LIMIT 1" % (options.userid)
    #print key_query_sql
    cur.execute(key_query_sql)
    key = cur.fetchone()

    try:
        if len(key) > 0:
           api_key = key[0]
           #print api_key
 
           # submit the API batch job with the key
           cmd = "python %s/tools/batch_processing/batch_submit_globus_genomics.py -t '%s' -k %s -u %s -i %s > %s" % (options.rootdir, options.goauth, api_key, url, options.table_file, options.log)
           print cmd

        # Run
        try:
            tmp_out = tempfile.NamedTemporaryFile().name
            tmp_stdout = open( tmp_out, 'wb' )
            tmp_err = tempfile.NamedTemporaryFile().name
            tmp_stderr = open( tmp_err, 'wb' )
            proc = subprocess.Popen( args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr )
            returncode = proc.wait()
            tmp_stderr.close()
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp_err, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stdout.close()
            tmp_stderr.close()
            if returncode != 0:
                raise Exception, stderr
            
            # TODO: look for errors in program output.
        except Exception, e:
            stop_err( 'Error in batch submission:\n' + str( e ) ) 

    except Exception, e:
        print "An error was encounteres: %s" % str(e)
        sys.exit(1)

except psycopg2.DatabaseError, e:
    print 'Error %s' % e    

