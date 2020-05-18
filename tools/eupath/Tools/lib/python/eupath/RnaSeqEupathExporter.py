import EupathExporter
import ReferenceGenome
import sys
import json
import os
import re

class RnaSeqExport(EupathExporter.Export):

    # Constants
    TYPE = "RnaSeq"
    VERSION = "1.0"

    def __init__(self, args):

        # print >> sys.stderr, "arguments:"
        # for i in range(0, len(args)):
        #   print >> sys.stderr, "args[" + str(i) + "] = " + args[i]
        # print >> sys.stderr, "end of arguments"

        EupathExporter.Export.__init__(self,
                                       RnaSeqExport.TYPE,
                                       RnaSeqExport.VERSION,
                                       "validateRnaSeq",
                                       args)

        if len(args) < 10:
            raise EupathExporter.ValidationException("The tool was passed too few arguments.")

        # grab first dataset provided ref genome
        self._initial_refGenome = args[10]
        # . . . and strandedness -- we're passed strandednessParam and we write strandedness
        strandednessParam = args[7]
        strandedness = strandednessParam

        self._datasetInfos = []

        # open manifest file
        manifestPath = "/tmp/manifest." + str(os.getpid()) + ".txt"
        manifest = open(manifestPath, "w+")

        # process variable number of [dataset refgenome] pairs.
        fileNumber = 0
        for i in range(8, len(args), 4):   # start on args[8], increment by 4
            # print >> sys.stderr, "args[" + str(i) + "] = " + args[i]
            samplename = args[i+1]
            suffix = args[i+3]
            filename = self.clean_file_name(re.sub(r"\s+", "_", samplename) + "." + suffix)

            fileNumber += 1
            if strandednessParam == "stranded":
                if filename[-4:] == ".txt":
                    filename = re.sub("forward", "one", re.sub("reverse", "two", filename))
                    samplename = re.sub("forward", "one", re.sub("reverse", "two", samplename))
                    strandedness = "sense" if (fileNumber % 2) == 1 else "antisense"
                else:
                    strandedness =  "firststrand" if (fileNumber % 2) == 1 else "secondstrand"

            self._datasetInfos.append({"name": filename, "path": args[i]})
            print >> manifest, samplename + "\t" + filename + "\t" + strandedness

        manifest.close()
        self._datasetInfos.append({"name": "manifest.txt", "path": manifestPath})

        self._refGenome = ReferenceGenome.Genome(args[10])

        # print >> sys.stderr, "datasetInfos: " + json.dumps(self._datasetInfos) + "<<- END OF datasetInfos"

    def identify_dependencies(self):
        """
        The appropriate dependency(ies) will be determined by the reference genome selected - only one for now
        """
        return [{
            "resourceIdentifier": self._refGenome.identifier,
            "resourceVersion": self._refGenome.version,
            "resourceDisplayName": self._refGenome.display_name
        }]

    def identify_projects(self):
        return [self._refGenome.project]

    def identify_dataset_files(self):
        """
        :return: A list containing the dataset files accompanied by their VEuPathDB designation.
        """
        return self._datasetInfos
