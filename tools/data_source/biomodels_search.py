#!/usr/bin/python

import cgi
import logging
from suds.client import Client
import suds
import cgi, cgitb 


print """content-type: text/html

<html><head><title>BioModels Search</title></head>
<body><h1>Search BioModels Database</h1>
This tool lets you search and import SBML models from the <a href="http://www.ebi.ac.uk/biomodels-main/" target="_blank">BioModels</a> database.<br>
<hr>
<form method=GET>
Search: <input name=query autofocus="true"> 
<select name="qtype">
  <option value="chebi">ChEBI (term)</option>
  <option value="chebiid">ChEBI (ID)</option>
  <option value="go">Gene Ontology (GO) (term)</option>
  <option value="goid">Gene Ontology (GO) (ID)</option>
  <option value="name" selected="true">Model Name</option>
  <option value="person">Author Name</option>
  <option value="publication">Publication (PMID, DOI, title/abstract)</option>
  <option value="taxonomy">Taxonomy (text)</option>
  <option value="taxonomyid">Taxonomy (ID)</option>
  <option value="uniprot">UniProt (text)</option>
  <option value="uniprotid">UniProt (ID)</option>
  <option value="uniprotids">UniProt Multiple (comma-separated ID list)</option>
</select> 

<input type=submit>
</form>

"""

#Retrieve results using web services

galaxy_url = 'https://cvrg2.galaxycloud.org'

url = 'http://www.ebi.ac.uk/biomodels-main/services/BioModelsWebServices?wsdl'

logging.basicConfig(level=logging.INFO)
logging.getLogger('suds.client').setLevel(logging.DEBUG)

try:
  client = Client(url)

  output = ""

  # Create instance of FieldStorage 
  form = cgi.FieldStorage() 

  # Get data from fields
  query = form.getvalue('query')
  qtype = form.getvalue('qtype')

  if query:
    #Get information
    if qtype=='chebi':
      result = client.service.getModelsIdByChEBI(query)
    elif qtype=='chebiid':
      result = client.service.getModelsIdByChEBIId(query)
    elif qtype=='go':
      result = client.service.getModelsIdByGO(query)
    elif qtype=='goid':
      result = client.service.getModelsIdByGOId(query)
    elif qtype=='name':
      result = client.service.getModelsIdByName(query)
    elif qtype=='person':
      result = client.service.getModelsIdByPerson(query)
    elif qtype=='publication':
      result = client.service.getModelsIdByPublication(query)
    elif qtype=='taxonomy':
      result = client.service.getModelsIdByTaxonomy(query)
    elif qtype=='taxonomyid':
      result = client.service.getModelsIdByTaxonomyId(query)
    elif qtype=='uniprot':
      result = client.service.getModelsIdByUniprot(query)
    elif qtype=='uniprotid':
      result = client.service.getModelsIdByUniprotId(query)
    elif qtype=='uniprotids':
      result = client.service.getModelsIdByUniprotIds(query.split(','))
    else:
      result = ''
      print "Error: invalid query type.<br>"
    
    if result:
      output = output + "<h2>Results for \"{0}\" ({1}):</h2>\n".format(query,len(result))
      for r in result:
        info_link = 'http://www.ebi.ac.uk/biomodels-main/' + r
        dl_link = 'http://www.ebi.ac.uk/biomodels-main/download?mid=' + r
        gal_link = galaxy_url + '/tool_runner?tool_id=biomodels&URL=' + dl_link
        output = output + r + " (<a href=\"{0}\" target=\"_blank\">info</a>, <a href=\"{1}\">import to Galaxy</a>)<br>\n".format(info_link,gal_link)
    else:
      output = output + "<br>No results for <b>\"%s\"</b><br>\n" % query
except suds.WebFault, e:
  print "An error occurred while querying the server.<br>"
  print e

# Generate the HTML response

#######################################################

print output

print """
</body></html>
"""
