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
        EupathExporter.Export.__init__(self,
                                       RnaSeqExport.TYPE,
                                       RnaSeqExport.VERSION,
                                       "validateRnaSeq",
                                       args)

        # beyond the standard 7 params, this exporter requires one or more pairs of args: dataset1 dataset1.refGenome
        # dataset2...
        if len(args) < 10:
            raise EupathExporter.ValidationException("The tool was passed too few arguments.")

        # grab first dataset provided ref genome
        self._initial_refGenome = args[9]

        self._datasetInfos = []

        # open manifest file
        manifestPath = "/tmp/manifest." + str(os.getpid()) + ".txt"
        manifest = open(manifestPath, "w+")

        # process variable number of [dataset refgenome] pairs.
        for i in range(7, len(args), 4):   # start on 8th arg, increment by 4
            # print >> sys.stderr, "args[" + str(i) + "] = " + args[i]
            filename = re.sub(r"\s+", "_", args[i+1]) + "." + args[i+3]
            self._datasetInfos.append({"name": filename, "path": args[i]})
            print >> manifest, args[i+1] + "\t" + filename + "\tunstranded"

        manifest.close()
        self._datasetInfos.append({"name": "manifest.txt", "path": manifestPath})

        # now override the dataset provided ref genome with the one obtained from the form assuming it is correctly
        # selected.  Otherwise throw an error.
        # if len(args[6].strip()) == 0:
        #     raise EupathExporter.ValidationException("A reference genome must be selected.")
        self._refGenome = ReferenceGenome.Genome(args[9])

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
        :return: A list containing the dataset files accompanied by their EuPathDB designation.
        """
        return self._datasetInfos
