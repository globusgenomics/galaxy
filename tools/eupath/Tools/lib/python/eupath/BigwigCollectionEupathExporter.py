import EupathExporter
import ReferenceGenome
import sys
import json


class BigwigCollectionExport(EupathExporter.Export):

    # Constants
    TYPE = "BigWig collection"
    VERSION = "1.0"

    def __init__(self, args):
        EupathExporter.Export.__init__(self,
                                       BigwigCollectionExport.TYPE,
                                       BigwigCollectionExport.VERSION,
                                       None,
                                       args)

        # print >> sys.stderr, "more debugging stuff!!"

        # beyond the standard 7 params, this exporter requires one or more pairs of args: dataset1 dataset1.refGenome
        # dataset2...
        if len(args) < 10:
            raise EupathExporter.ValidationException("The tool was passed too few arguments.")

        # grab first dataset provided ref genome
        self._initial_refGenome = args[9]

        self._datasetInfos = []
        
        # process variable number of [dataset refgenome] pairs.
        # confirm that all dataset provided ref genomes are identical.
        for i in range(7, len(args), 3):   # start on 8th arg, increment by 3
            # if args[i+2] != self._initial_refGenome:
            #     raise EupathExporter.ValidationException("All provided bigwig datasets must have the same reference genome.  Found " + self._initial_refGenome + " and " + args[i+2])
            self._datasetInfos.append({"name": args[i+1], "path": args[i]})
            print >> sys.stderr, "name: " + args[i+1]
            print >> sys.stderr, "path: " + args[i]

        # now override the dataset provided ref genome with the one obtained from the form assuming it is correctly
        # selected.  Otherwise throw an error.
        if len(args[6].strip()) == 0:
            raise EupathExporter.ValidationException("A reference genome must be selected.")
        self._refGenome = ReferenceGenome.Genome(args[6])

        print >> sys.stderr, "datasetInfos: " + json.dumps(self._datasetInfos) + "<<- END OF datasetInfos"

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
