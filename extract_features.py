#!/usr/bin/python
''' This script extracts sequence regions from annotated Genbank file by the exact name of the region.
	We use PFAM domain names.
	To check what domains are available, use -l flag.
'''

from optparse import OptionParser
from Bio import SeqIO

parser = OptionParser()
parser.add_option("-f", "--sequence-file", dest="seq_file", default=None,
                  help="annotated genbank file, default: None")
parser.add_option("-n", "--names", dest="domain_names",  default="",
                  help="QUOTED string containing domain names")
# parser.add_option("-i", "--domain-number", dest="domain_number", default=None,
#                   help="domain_number, default: None")
parser.add_option("-c", "--coords-file", dest="coords_file", default=None,
                  help="coordinate file, default: None")
parser.add_option("-l", "--list", action="store_true", dest="is_listed", default=False,
				  help="print all domain names in the file")

(opts, args) = parser.parse_args()
seq_file = opts.seq_file
domain_names = opts.domain_names.split()
coords_file = opts.coords_file
is_listed = opts.is_listed


in_handle = open(seq_file, "rU")
records = list(SeqIO.parse(in_handle, "genbank"))
in_handle.close()

coords_file_string = ""

selected_domains = []
unique_domains = []
for r in records:
	for d in r.features:
		if d.type == "Region":
			if 'name' in d.qualifiers:
				if d.qualifiers['name'][0] not in unique_domains:
					unique_domains.extend(d.qualifiers['name'])
				coords_file_string += " ".join([ d.qualifiers['name'][0], str(d.location.start), str(d.location.end), r.name ]) + "\n"
				for dn in domain_names:
					if d.qualifiers['name'][0] == dn:
						selected_domains.append(d)
						if not is_listed:
							print ">" + r.id + "/" + str(d.location.start) + "-" + str(d.location.end) + " " + r.description + "\n" + d.extract(r.seq)
if is_listed:
	print "Unique domains in the file:\n", "\n".join(unique_domains)

if coords_file:
	coords_handle = open(coords_file, "w")
	coords_handle.write(coords_file_string)
	coords_handle.close()