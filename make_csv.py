#!/usr/bin/python

import sys
from optparse import OptionParser


def make_csv(result_file, output_file):
    ''' Reads prediction results and writes formatted csv to output_file in order to display several paralogs in one cell
    '''
    in_handle = open(result_file, "rU")
    lines = in_handle.readlines()
    in_handle.close()
    i = 1
    out_handle = open(output_file, "w")
    instances = []
    while i < len(lines):
        l = lines[i]
        l_prev = lines[i - 1]
        splitted = l.rstrip().split("\t")
        splitted_prev = l_prev.rstrip().split("\t")
        if len(splitted) > 1:
            splitted[1] = splitted[1].split("|")[1]
        if len(splitted_prev) > 1:
            splitted_prev[1] = splitted_prev[1].split("|")[1]
        instances.append(splitted_prev)
        if splitted[0] == splitted_prev[0] and i == len(lines) - 1:
        	instances.append(splitted)
        if splitted[0] == splitted_prev[0] and i != len(lines) - 1:
            i += 1
        else:
            if len(instances) > 1:
                out_handle.write( instances[0][0] + ",\"" + "\n".join([p[1] for p in instances]) + "\",\"" + "\n".join([p[2] for p in instances]) + "\",\"" + "\n".join([p[3] for p in instances]) + "\"\n")
            else:
                out_handle.write( ",".join(instances[0]) + "\n" )
            i += 1
            instances = []
    out_handle.close()
    return 0


parser = OptionParser()
parser.add_option("-i", "--in-file", dest="result_file", default=None,
                  help="file with homolog predictions")
parser.add_option("-o", "--out-file", dest="output_file", default=None,
                  help="formatted output csv file")

# read options and arguments passed to the script
(opts, args) = parser.parse_args()
result_file = opts.result_file
output_file = opts.output_file

make_csv(result_file, output_file)