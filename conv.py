#!/usr/bin/python
''' This script convert between sequence formats.
	Ideally, it should convert everything to everything.
    Now it converts Uniprot fasta to genbank and back.
'''

from optparse import OptionParser
from Bio import SeqIO
from Bio.Alphabet import generic_protein
import re

parser = OptionParser()
(opts, args) = parser.parse_args()

for arg in args:
    input_handle = open(arg, "rU")
    ext = arg.split(".")[-1]
    name = ".".join(arg.split(".")[:-1])
    if ext == "gb":
        sequences = list(SeqIO.parse(input_handle, "genbank"))
        output_handle = open(name + ".fasta", "w")
        count = SeqIO.write(sequences, output_handle, "fasta")
        output_handle.close()
    elif ext == "fasta":
        records = list(SeqIO.parse(input_handle, "fasta"))
        count = 0
        converted_records = []
        for r in records:
            if len(r.id.split("|")) > 1:
                r.annotations['accession'] = r.id.split("|")[1]
                r.id = r.id.split("|")[2]
                r.name = r.id
                r.seq.alphabet = generic_protein
            converted_records.append(r)
            organism_regex = re.compile(r'OS=([\w\s\/\(\)\.\d\-]+)\w\w=')
            for match in organism_regex.finditer(r.description):
                r.annotations['organism'] = match.group(1)
                r.annotations['source'] = match.group(1)
        output_handle = open(name + ".gb", "w")
        count = SeqIO.write(converted_records, output_handle, "genbank")
        output_handle.close()
    else:
        print "Unknown file extension"

    print "Converted %i records" % count
    input_handle.close()
