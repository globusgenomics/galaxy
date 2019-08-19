#!/usr/bin/python

import EupathExporter
import ReferenceGenome


class GeneListExport(EupathExporter.Export):
    """
    This class is a specialized version of the Galaxy to EuPathDB dataset export tool.  This tool's
    specialty is furnishing user's gene list data to EuPathDB.  As with all specialty export tools, this
    tool implements 3 abstract classes.
    """

    # Name given to this type of dataset and to the expected file
    GENE_LIST_TYPE = "GeneList"
    GENE_LIST_VERSION = "1.0"
    GENE_LIST_FILE = "genelist.txt"

    # The validation script to be applied to the dataset files.  A failed validation should
    # return in a system exit status of other than 0.
    GENE_LIST_VALIDATION_SCRIPT = "validateGeneList"

    def __init__(self, args):
        """
        Initializes the gene list export class with the parameters needed to accomplish the particular
        type of export.
        :param args: parameters provided from tool form
        """
        EupathExporter.Export.__init__(self,
                                       GeneListExport.GENE_LIST_TYPE,
                                       GeneListExport.GENE_LIST_VERSION,
                                       GeneListExport.GENE_LIST_VALIDATION_SCRIPT,
                                       args)

        # For the gene list export, three parameters beyond generic 7 are required.
        if len(args) < 10:
            raise EupathExporter.ValidationException("The tool was passed an insufficient numbers of arguments.")

        # Data for the input given by the user
        self._dataset_file_path = args[7]

        # Overriding the dataset genome reference with that provided via the form.
        if len(args[6].strip()) == 0:
            raise EupathExporter.ValidationException("A reference genome must be selected.")
        self._genome = ReferenceGenome.Genome(args[6])
        self._datatype = args[9]

    def identify_dependencies(self):
        """
        The appropriate dependency(ies) will be determined by the reference genome selected - only one for now
        The EuPathDB reference genomes will have a project id, a EuPath release number, and a genome description
        all separated by a dash in the first instance and an underscore in the second instance.
        :return: list containing the single dependency with the component parts parsed out (only one for now)
        """
        return [{
            "resourceIdentifier": self._genome.identifier,
            "resourceVersion": self._genome.version,
            "resourceDisplayName": self._genome.display_name
        }]

    def identify_projects(self):
        """
        The appropriate project(s) will be determined by the reference genome selected - only one for now
        The project name must be listed in the SUPPORTED_PROJECTS array.  Failure to find it will be
        regarded as a validation exception.
        :return: list containing the single relevant EuPath project (only one for now)
        """
        return [self._genome.project]

    def identify_dataset_files(self):
        """
        The user provided gene list file is combined with the name EuPathDB expects
        for such a file
        :return: A list containing the single dataset file accompanied by its EuPathDB designation.
        """
        return [{"name": self.GENE_LIST_FILE, "path": self._dataset_file_path}]
