#!/usr/bin/env python
# Alex Rodriguez

"""
A wrapper script for running Tabix.
"""

import sys, optparse, os, tempfile, subprocess, shutil

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-P', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to tabix, without any modification.' )
    parser.add_option( '-o', '--output', dest='output', type="string", help='output file' )
    parser.add_option( '-i', '--input', dest='input', type="string", help='input file' )
    (options, args) = parser.parse_args()

    tmp_dir = tempfile.mkdtemp( prefix='tabix-temp' )
    tmp_tabix_input = tmp_dir + "/tabix_input"
    tmp_tabix_output = tmp_tabix_input  + ".tbi"
    shutil.copyfile(options.input, tmp_tabix_input)

    cmd = "tabix " + tmp_tabix_input + " "
    if options.pass_through_options:
        cmd += '  '.join( options.pass_through_options )
    else:
        cmd += ' '

    print cmd
    proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir )
    return_code = proc.wait()

    shutil.copyfile(tmp_tabix_output, options.output)

if __name__=="__main__": __main__()
