from common import get2
import sys

def main(address):
    createFile(address)

def createFile(address):
    try:
        tools = getTools(address)
        valid = []
        for dicts in tools:
            if len(get2(address.strip('/') + "/api/repositories/get_ordered_installable_revisions?name=" + dicts['name'] + "&owner=" + dicts['owner'])) == 0:
                continue
            else:
                valid.append(dicts)
        opfile= open("output.txt", 'w')
        if address == "http://toolshed.g2.bx.psu.edu/":
            opfile.write("Success! Tool info exported to toolshedinfo.txt in tool-data folder: \n")
        elif address == "http://testtoolshed.g2.bx.psu.edu/":
            opfile.write("Success! Tool info exported to testtoolshedinfo.txt in tool-data folder: \n")
        else:
            opfile.write("Success! Tool info exported to privatetoolshedinfo.txt in tool-data folder: \n")
        opfile.close()
        for dictionary in valid:
            print dictionary['name'] + " by " + dictionary['owner'] + " " + "\n"
    except:
        print "Invalid tool shed URL: " + address

def getTools(toolshedurl):
    return get2(toolshedurl.strip('/') + "/api/repositories")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()