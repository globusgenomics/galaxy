import EupathExporter
import ReferenceGenome
import sys


class BigwigFilesExport(EupathExporter.Export):

    # Constants
    TYPE = "BigwigFiles"
    VERSION = "1.0"

    def __init__(self, args):
        EupathExporter.Export.__init__(self,
                                       BigwigFilesExport.TYPE,
                                       BigwigFilesExport.VERSION,
                                       None,
                                       args)

        # beyond the standard 7 params, this exporter requires one or more pairs of args: dataset1 dataset1.refGenome
        # dataset2...
        if len(args) < 10:
            raise EupathExporter.ValidationException("The tool was passed an insufficient numbers of arguments.")

        # grab first dataset provided ref genome
        self._initial_refGenome = args[9]

        self._datasetInfos = []
        
        # process variable number of [dataset refgenome] pairs.
        # confirm that all dataset provided ref genomes are identical.
        for i in range(7, len(args), 3):   # start on 8th arg, increment by 3
            #if not args[i+1].endswith(".bigwig") and not args[i+1].endswith(".bw"):
            #    raise EupathExporter.ValidationException("All datafiles must have either the .bigwig or .bw extension.")
            if args[i+2] != self._initial_refGenome:
                raise EupathExporter.ValidationException("All provided bigwig datasets must have the same reference genome.  Found " + self._initial_refGenome + " and " + args[i+2])
            self._datasetInfos.append({"name": args[i+1], "path": args[i]})

        # now override the dataset provided ref genome with the one obtained from the form assuming it is correctly
        # selected.  Otherwise throw an error.
        if len(args[6].strip()) == 0:
            raise EupathExporter.ValidationException("A reference genome must be selected.")
        self._refGenome = ReferenceGenome.Genome(args[6])


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
