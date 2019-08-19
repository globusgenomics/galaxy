#!/usr/bin/python
import sys
import json
import re
import os

class ValidationException(Exception):
  pass

def __main__():
    """
    Simple validation of a gene list dataset.  The validation checks:
    1.  that exactly one dataset file is provided and
    2.  that no lines in the file contain embedded whitespace
    If any validation issue is discovered, a validation error is returned.
    """
    error = None
    dataset_file_json = sys.stdin.read()
    dataset_file = json.loads(dataset_file_json)
    if len(dataset_file) != 1:
        error = "A gene list dataset should have only one file."
    else:
        with open(dataset_file[0]['path'],'r') as source_file:
            line = source_file.readline().strip()
            if re.search(r"\s", line):
                error = "No lines in a gene list file should contain embedded whitespace."
    if error is not None:
        print >> sys.stderr, "Error: " + error
        sys.exit(1)

if __name__ == "__main__": __main__()