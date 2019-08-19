import argparse, time
from common import submit, get2

def main( options ):
    """Collect all user data and installs the tools via the Galaxy API."""
    start_time = time.time() #creates time stamp for entire job
    tools=getTools(options.tool_shed_url)
    errors = [] #creates errors list

    if options.tool_shed_url == "": #returns errors if there are inputs missing
        print "Error: no tool shed URL input."
        return
    if len(options.local_url) == 0:
        print "Error: no galaxy URL input."
        return
    if len(options.api) == 0:
        print "Error: no galaxy API key input."
        return
    if len(options.name) == 0 and len(options.shedlist) == 0 and len(options.testshedlist) == 0:
        print "Error: no tool input."
        return
    if len(options.api) != len(options.local_url):
        print "Error: too few galaxy API keys or galaxy URLs input."
        return
    if len(options.name) != len(options.owner):
        print "Error: too few owners or tool names."
        return

    if options.shedlist: #runs installToolShedList method to install tools via Galaxy dynamic lists
        errors = installToolShedList(options.shedlist, 'http://toolshed.g2.bx.psu.edu/', errors)
    if options.testshedlist:
        errors = installToolShedList(options.testshedlist, 'http://testtoolshed.g2.bx.psu.edu/', errors)
    if options.privateshedlist:
        errors = installToolShedList(options.privateshedlist, options.tool_shed_url, errors)

    for i in range(len(options.local_url)): #loops through all instances
        for x in range(len(options.name)): #loops through name arguments
            found = False #creates boolean to determine if tool is found in tool shed
            for dictionary in tools: #loops through all the tools in the tool shed
                if options.name[x] == dictionary['name']:
                    found = True
                    if dictionary['owner'] == options.owner[x] or options.owner[x] == "none":
                        dicts=getDict(dictionary)
                        if dicts['changeset_revision'] == "":
                            errors.append(dicts['name'] + " by " + dicts['owner'] + " failed to install on instance: " + options.local_url[i] + ": no repo revision number found.")
                            continue #returns error if there is no revision number
                        errors.append(installTool(i, dicts)) #installs tool and adds errors to errors list
                    else:
                        errors.append("The tool " + options.name[x] + " does not have owner " + options.owner[x] + " (or multiple owners for same tool).")
            if found == False:
                errors.append(options.name[x] + " tool not found in toolshed: " + options.tool_shed_url + ".")
                    
    print "---------------------------------------------"
    print 'REPORT:' #prints error report from errors list with time stamps
    for error in errors:
        print "\n" + error
    print "\nTool runtime: " + str(round((time.time() - start_time),3)) + " seconds."
    print "---------------------------------------------" 

def getTools(toolshedurl="http://toolshed.g2.bx.psu.edu/"): #retrieves list of dictionaries with tool name and owners from input URl
    tools=[]
    toollists = get2(toolshedurl.strip('/') + "/api/repositories")
    for dictionary in toollists:
        dicts = {}
        dicts['name'] = dictionary['name']
        dicts['owner'] = dictionary['owner']
        tools.append(dicts)
    return tools

def getDict(dictionary): #gets rest of required tool info needed to install
    dicts={}
    dicts['name'] = dictionary['name']
    dicts['owner'] = dictionary['owner']
    dicts['changeset_revision'] = getRepo(dicts)
    if options.install_tool_dependencies:
        dicts["install_tool_dependencies"] = "True"
    if options.install_repository_dependencies:
        dicts['install_repository_dependencies'] = "True"
    if options.new_tool_panel_section_label != "":
        dicts['new_tool_panel_section_label'] = options.new_tool_panel_section_label
    dicts['tool_shed_url'] = options.tool_shed_url
    return dicts

def getRepo(dicts): #gets latest repository install number
    repos =  get2(options.tool_shed_url.strip('/') + "/api/repositories/get_ordered_installable_revisions?name=" + dicts['name'] + "&owner=" + dicts['owner'])[-1:]
    if len(repos) == 0:
        return ""
    else:
        return repos[0]

def installToolShedList(arguments, URL, errors): #installs tools from Galaxy dynamic lists
    tools = getTools(URL)
    inputs= arguments
    new_inputs=[]
    names = []
    owners = []
    errors = errors
    for element in inputs: #removes "by" and commas from argument list
        if element == "by":
            continue
        new_element = element.replace(",",  "")
        new_inputs.append(new_element)
    for i in range(len(new_inputs)): #separates arguments into names and owners
        if i % 2 == 0:
            names.append(new_inputs[i])
        else:
            owners.append(new_inputs[i])
    for i in range(len(options.local_url)): #see main, getDicts, and getRepo comments
        for x in range(len(names)):
            for dictionary in tools:
                if dictionary['name'] == names[x] and dictionary['owner'] == owners[x]:
                    dicts={}
                    dicts['name'] = dictionary['name']
                    dicts['owner'] = dictionary['owner']
                    if options.install_tool_dependencies:
                        dicts["install_tool_dependencies"] = "True"
                    if options.install_repository_dependencies:
                        dicts['install_repository_dependencies'] = "True"
                    if options.new_tool_panel_section_label != "":
                        dicts['new_tool_panel_section_label'] = options.new_tool_panel_section_label
                    dicts['tool_shed_url'] = URL
                    repo = get2(URL.strip('/') + "/api/repositories/get_ordered_installable_revisions?name=" + dicts['name'] + "&owner=" + dicts['owner'])[-1:]
                    if len(repo) == 0:
                        errors.append(dicts['name'] + " by " + dicts['owner'] + " failed to install on instance: " + options.local_url[i] + ": no repo revision number found.")
                        continue
                    dicts['changeset_revision'] = repo[0]
                    errors.append(installTool(i, dicts))
    return errors

def installTool(instance, toolinfo): #makes submit request to install tool
    try:
        tool_time=time.time() #creates time stamp for specific tool
        submit(options.api[instance], '%s%s' % (options.local_url[instance].strip('/'), '/api/tool_shed_repositories/new/install_repository_revision'), toolinfo) #makes actual submit request to install tool
        if options.new_tool_panel_section_label: #print extra information
            print "  tool_panel: " + options.new_tool_panel_section_label
        else:
            print "  tool_panel: None"
        if options.install_tool_dependencies:
            print "  installed_tool_dependencies: True"
        else:
            print "  installed_tool_dependencies: False"
        if options.install_repository_dependencies:
            print "  installed_repository_dependencies: True"
        else: 
            print "  installed_repository_dependencies: False"
        print "  instance_installed_on: " + options.local_url[instance]
        return toolinfo['name'] +  " by " + toolinfo['owner'] + " successfully installed on instance: " + options.local_url[instance] + " in " + str(round((time.time() - tool_time),3)) + " seconds"
    except:
        return toolinfo['name'] + " failed to install on instance: " + options.local_url[instance] + ": error or tool already installed"

if __name__ == '__main__': #options parser
    parser = argparse.ArgumentParser( description='Installation of tool shed repositories via the Galaxy API.' )
    parser.add_argument( "-u", "--url", dest="tool_shed_url", required=True, help="Tool Shed URL" )
    parser.add_argument( "-a", "--api", nargs='*', dest="api", required=True, help="API Key" )
    parser.add_argument( "-l", "--local", nargs="*", dest="local_url", required=True, help="URL of the galaxy instance." )
    parser.add_argument( "-n", "--name", nargs="*", required=True, help="Repository name." )
    parser.add_argument( "-o", "--owner", nargs="*",help="Repository owner." )
    parser.add_argument( "-s", "--shedlist", nargs="*", dest="shedlist", help="For comprehensive input list (Galaxy tool shed) from galaxy ONLY.")
    parser.add_argument( "-t", "--testshedlist", nargs="*", dest="testshedlist", help="For comprehensive input list (Galaxy test tool shed) from galaxy ONLY.")
    parser.add_argument( "-p", "--privateshedlist", nargs="*", dest="privateshedlist", help="For comprehensive input list (Galaxy tool shed) from galaxy ONLY.")
    parser.add_argument( "--panel-section-name", dest="new_tool_panel_section_label", help="New tool panel section label. If specified a new tool section will be created." )
    parser.add_argument( "--repository-deps", dest="install_repository_dependencies", action="store_true", default=False, help="Install repository dependencies. [False]")
    parser.add_argument( "--tool-deps", dest="install_tool_dependencies", action="store_true", default=False, help="Install tool dependencies. [False]" )
    options = parser.parse_args()
    main( options )