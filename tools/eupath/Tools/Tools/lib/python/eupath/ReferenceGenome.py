#!/usr/bin/python

import EupathExporter
import re


class Genome:

    def __init__(self, reference_genome):
        """
        Teases out from the user's parameter, the reference genome information used in the construction of dependency
        data.  The reference genome parameter should be of the form: ProjectId-EupathBuildNumber_Strain_Genome
        :param reference_genome: the reference genome parameter provided by the user
        """

        # Insure the the reference genome matches the pattern for Eupath originated reference genomes.
        if not reference_genome and not re.match(r'^.+-\d+_.+_Genome$', reference_genome, flags=0):
            raise EupathExporter.ValidationException(
                "A syntactically correct reference genome is required for exports to EuPathDB.")
        self._identifier = reference_genome
        self._project = reference_genome[0:reference_genome.index("-")]
        sans_project = reference_genome[reference_genome.index("-") + 1:]
        components = sans_project.split("_")
        self._version = components[0]
        self._display_name = components[1] + " Genome"

    @property
    def project(self):
        return self._project

    @property
    def version(self):
        return self._version

    @property
    def display_name(self):
        return self._display_name

    @property
    def identifier(self):
        return self._identifier
