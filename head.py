#!/usr/bin/python
''' This script extracts protein head, starting from the N-terminus to the first known domain of the annotated *.gb file
'''
import os
import sys
from Bio import SeqIO

try:
	filename = sys.argv[1]
	if os.path.isfile(filename):
		in_handle = open(filename, "rU")
		records = list(SeqIO.parse(in_handle, "genbank"))
		in_handle.close()
		for record in records:
			i = 0
			for f in record.features:
				if i == 0:
					first_feature_start = f.location.start
				if int(f.location.start) < first_feature_start:
					first_feature_start = f.location.start
				i += 1
			print ">" + record.name + "/1-" + str(first_feature_start) + "OS=" + record.annotations['organism']
			print record.seq[0:first_feature_start + 1]
except:
	print "This script extracts protein head, starting from the N-terminus to the first known domain of the annotated *.gb file"
