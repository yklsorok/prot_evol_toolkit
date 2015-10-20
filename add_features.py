#!/usr/bin/python
''' This script adds the feature from the coordinate file.
    Format of the coordinate file:
    feature_name start end source
    motif 12 20 A0A088A7Y6
'''

import os, sys
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def make_protein_feature(feature_name, feature_start, feature_end, feature_type):
    ''' Returns sequence feature, using start, end, name and type as input
    '''
    feature = SeqFeature(FeatureLocation(int(feature_start), int(feature_end)), type=feature_type)
    if feature_type == "Region":
        feature.qualifiers = {'name': [feature_name]}
    return feature


parser = OptionParser()
parser.add_option("-f", "--file", dest="seq_file", default=None,
                  help="genbank file")
parser.add_option("-c", "--coord-file", dest="coord_file", default=None,
                  help="coordinate file")
parser.add_option("-r", "--range", dest="seq_range", default="1..10",
                  help="range in format 1..100, full by default")
parser.add_option("-s", "--source", dest="seq_source", default="A0A088A7Y6",
                  help="source name for sequence")
parser.add_option("-n", "--feature_name", dest="feature_name", default="motif",
                  help="new source name for sequence")

# read options and arguments passed to the script
(opts, args) = parser.parse_args()
seq_file = opts.seq_file
coord_file = opts.coord_file
seq_range = opts.seq_range.split("..")
feature_name = opts.feature_name
seq_source = opts.seq_source

new_features = {}
if coord_file:
    in_handle = open(coord_file, "rU")
    lines = in_handle.readlines()
    in_handle.close()
    for l in lines:
        splitted = l.rstrip().split(" ")
        f_name = splitted[0]
        f_source = splitted[3]
        f_start = int(splitted[1]) - 1
        f_end = int(splitted[2])
        d = make_protein_feature(f_name, f_start, f_end, "Region")
        if f_source not in new_features.keys():
            new_features[f_source] = []
        new_features[f_source].append(d)

in_handle = open(seq_file, "rU")
records = list(SeqIO.parse(in_handle, "genbank"))
in_handle.close()

if seq_range and seq_source and feature_name:
    d = make_protein_feature(feature_name, seq_range[0], seq_range[1], "Region")

for r in records:
    for s, f in new_features.items():
        if r.name == s:
            r.features.extend(f)
out_handle = open(seq_file, "w")
SeqIO.write(records, out_handle, "genbank")
out_handle.close()

