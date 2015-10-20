#!/usr/bin/python
''' This generates a bunch of images to represent domain/motif structure using PSImage.
	For convenience, it is possible to produce single pdf file with convert.
	The script requires annotated Genbank file and Internet connection.
'''

from Bio import SeqIO
import sys
import os
from urllib2 import urlopen, URLError, HTTPError
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--sequence-file", dest="seq_file", default=None,
                  help="annotated genbank file, default: None")
parser.add_option("-o", "--output", dest="output_file",  default="output.pdf",
                  help="output pdf file")


(opts, args) = parser.parse_args()
seq_file = opts.seq_file
output_file = opts.output_file

domain_types = {
# first number: shape
# 1 - square, 2 - circle, 3 - right pentagonal, 4 - left pentagonal, 5 - hexagonal, 6 - rotated hexagonal
# second number: color
# 1 - orange, 2 - green, 3 - blue, 4 - gray
	"CCR": ["1", "1"],
	"EF-hand_4": ["3", "2"],
	"FCH": ["2", "3"],
	"MHD": ["3", "4"],
    "PID": ["3", "2"],
    "PID_2": ["3", "2"],
    "PH_13": ["1", "4"],
    "SH3": ["1", "2"],
    "BAR": ["1", "4"],
    "ANTH": ["5", "1"],
    "ENTH": ["5", "1"],
    "AP2_binding": ["1", "1"],
    "NPF": ["1", "2"],
    "Clathrin_binding": ["1", "4"],
    "Arrestin_N": ["2", "2"],
    "G_DYNAMIN_2": ["2", "5"],
    "GED": ["1", "1"],

}


def dlfile(url, filename):
    # Open the url
    try:
        f = urlopen(url)
        print "Downloading " + filename

        # Open our local file for writing
        with open(os.path.basename(filename), "wb") as local_file:
            local_file.write(f.read())

    # Handle errors
    except HTTPError, e:
        print "HTTP Error: ", e.code, url
    except URLError, e:
        print "URL Error: ", e.reason, url

in_handle = open(seq_file, "rU")
records = list(SeqIO.parse(in_handle, "genbank"))
in_handle.close()

print "We're cleaning up..."
os.system("rm *.png")

i = 0
for r in records:
	features = []
	for f in r.features:
		f_name = f.qualifiers['name'][0]
		if f_name in domain_types.keys():
			f_type = domain_types[f_name][0]
			f_color = domain_types[f_name][1]
		else:
			f_type = "2"
			f_color = "3"
		features.append(",".join([str(f.location.start), str(f.location.end), f_type + "_" + f_color, f_name]))
	f_string = "+".join(features)
	prosite_url = "http://prosite.expasy.org/cgi-bin/prosite/PSImage.cgi?hit=" + f_string + "&len=" + str(len(r.seq)) + "&hscale=0.8"
	dlfile(prosite_url, str(i).zfill(3) + "_" + r.name + ".png")
	i += 1
os.system("convert *.png -gravity North -annotate 0 '%f' " + output_file)