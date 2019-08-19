def get_analysis_status = {

	return 'completed'

}

username='admin'
password = 'pass'

//Replace server and IP with accurate information.
processURI = 'http://ec2-50-19-133-33.compute-1.amazonaws.com:8080/api/v2/processes/ATT-SA2-120925-24-3909'

processNode = GLSGeneusRestApiUtils.httpGET(processURI, username, password)

Analysis_Status_UDF = processNode.'udf:field'.find { it.@name== 'Analysis Status' }

Analysis_Status = get_analysis_status()

//if (Analysis_Status_UDF != null) {

Analysis_Status_UDF.setValue(Analysis_Status)

//} 


returnNode = GLSGeneusRestApiUtils.httpPUT(processNode, processNode.@uri, username, password)

println GLSGeneusRestApiUtils.nodeToXmlString(returnNode)

