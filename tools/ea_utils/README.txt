== ea-utils Galaxy Wrapper ==

This is a Galaxy wrapper for some ea-utils tools, fastq-join and sam-stats.

** Installation **

Installation from a tool shed provides the necessary tool dependencies.

Otherwise, make sure fastq-join and sam-stats are in the path.
Move the test data files to your galaxy root test-data.
Move the xml files to a subdirectory of your tools directory and add lines in tool_conf.xml to point to them.
Restart the Galaxy server.

** Attribution **

The ea-utils package and associated documentation can be found at: http://code.google.com/p/ea-utils/

The galaxy wrapper code was written by Lance Parsons (lparsons@princeton.edu), Lewis-Sigler Institute for Integrative Genomics, Princeton University.
The code is housed on BitBucket at: https://bitbucket.org/lance_parsons/ea_utils_galaxy_wrapper
