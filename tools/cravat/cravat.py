import requests
import json
import sys
from __builtin__ import False
#Gets the input and output from galaxy
input_filename = sys.argv[1]
output_filename = sys.argv[2]

#opens each file, in to read, out to write
in_file = open(input_filename, "r")
out_file = open(output_filename, "w")


#sets replacements to replace each space in genomic coordinates with an underscore to run with the query
replacements = {' ':'_'}
#so we only print out the Keys once
write_header = True 

#loops through the input file line by line
for line in in_file:
    #strips the input line of \n (new line) and replaces every space with an underscore   
    line = "_".join( line.split()[1:6] )
    print line
    #gets request from CRAVAT server with the inputed mutation line
    call = requests.get('http://staging.cravat.us/CRAVAT/rest/service/query', params={'mutation': line} )
    #puts the string of data into a json dictionary
    try:
        json_data = json.loads(call.text)
    except:
        continue
    #manually sets the order of the Keys to the same as CRAVAT Server
    keys = ["Chromosome","Position","Strand","Reference base(s)","Alternate base(s)","HUGO symbol",
            "Sequence ontology transcript","Sequence ontology protein change","Sequence ontology",
            "Sequence ontology all transcripts","ExAC total allele frequency",
            "ExAC allele frequency (African/African American)","ExAC allele frequency (Latino)",
            "ExAC allele frequency (East Asian)","ExAC allele frequency (Finnish)",
            "ExAC allele frequency (Non-Finnish European)","ExAC allele frequency (Other)",
            "ExAC allele frequency (South Asian)", "1000 Genomes allele frequency",
            "ESP6500 allele frequency (European American)","ESP6500 allele frequency (African American)",
            "Transcript in COSMIC","Protein sequence change in COSMIC",
            "Occurrences in COSMIC [exact nucleotide change]","Mappability Warning","Driver Genes",
            "TARGET","dbSNP","MuPIT Link"]
    print json_data
    
    #Spit out first 8 or 9, then loop through rest and print out
    for key in json_data:
        if key not in keys:
            keys.append(key)
    
    #used so we only print out the Keys once
    if write_header == True:
        #writes out the keys of the dictionary
        out_file.write('\t'.join(keys) + '\n')
        write_header = False
        
    #sets value to the first value in the first key
    value = json_data[keys[0]]
    #actually writes out the value
    out_file.write(value)
    #print "key[" + key[0] + "] value[" + str(value) + "]"
    #sets value to the second key
    value = json_data[keys[1]]
    out_file.write('\t' + value)
    #print "key[" + key[1] + "] value[" + str(value) + "]"
    
    
    #loops through all other values for each key
    for key in keys[2:]:
        #strips the value
        if key in json_data:
            value = json_data[key].strip()
            #another try, except statement to convert the rest of the values to floats, and then round them to four decimals
            try:
                value = float(value)
                value = '%.4f'%value
            except:
                pass
            #writes out the value with a tab after for galaxy formatting 
            out_file.write("\t" + str(value))
        #print for debugging
        print "key[" + key + "] value[" + str(value) + "]"
    #creates a new line for the next set of values
    out_file.write('\n')
    

#closes both files    
in_file.close()
out_file.close()
