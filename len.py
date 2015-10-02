#!/usr/bin/python
''' This script computes the length of one or more sequences in the provided file.
	It returns sequence name and length.
'''

from Bio import SeqIO
from optparse import OptionParser
import os.path


def seq_len(filename):
    if os.path.isfile(filename):
        format = filename.split(".")[-1]
        if format == "fasta" or format == "gb":
            in_handle = open(filename, "rU")
            records = list(SeqIO.parse(in_handle, format))
            in_handle.close()
            return([[r.id, len(r.seq)] for r in records])
        else:  
            print "Unknown file extension"
    else:
        print "File ", filename, " not found"

parser = OptionParser()

(opts, args) = parser.parse_args()
for filename in args:
    # print filename
    for l in seq_len(filename):
    	print l[0], " \t", str(l[1])
