#!/usr/bin/python
''' This script extracts protein tail, starting from the end of last known domain from annotated *.gb file
'''
import os
import sys
from Bio import SeqIO

cutoff = 25
try:
	filename = sys.argv[1]
	if os.path.isfile(filename):
		in_handle = open(filename, "rU")
		records = list(SeqIO.parse(in_handle, "genbank"))
		in_handle.close()
		for record in records:
			last_feature_end = 0
			for f in record.features:
				if f.location.end > last_feature_end and (f.location.end - f.location.start) >= cutoff:
					last_feature_end = f.location.end + 1
			print ">" + record.name + "/" + str(last_feature_end) + "-" + str(len(record.seq)) + "OS=" + record.annotations['organism']
			print record.seq[last_feature_end - 1:len(record.seq)]
except:
	print "This script extracts protein tail, starting from the end of last known domain from annotated *.gb file"
