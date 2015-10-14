#!/usr/bin/python
''' This script helps to fix errors of protein annotation.
	The hypothesis behind this script is that orthologous multidomain proteins existed long long ago.
	If the protein lacks several domains that exist in its vertebrate ortholog, 
	the script searches protein domains in the adjacent genome region.
	Protein sequence, genome region sequence and domain name are used as input.
	Currently the sequence is sent to SMART for prediction
'''
# 1. find orfs in sequence
# 2. prepare big fasta file with proteins
# 3. run SMART prediction
# 4. parse results

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
from optparse import OptionParser
from StringIO import StringIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML

parser = OptionParser()
parser.add_option("-d", "--domain", dest="domain", default="PH",
                  help="SMART domain name to search for")
parser.add_option("-f", "--genomic-file", dest="filename",
                  help="fasta file with genomic sequence")
parser.add_option("-r", "--range", dest="range", default="1...",
                  help="range in format 1..10000, full by default")
parser.add_option("-l", "--min-length", dest="min_pro_len", default=50,
                  help="minimal protein length in ORF")


(opts, args) = parser.parse_args()
##########################################
# initial parameters
##########################################
# domain name we are searching for
domain_name = opts.domain
# file with genomic sequence in fasta format
genomic_filename = opts.filename
# range
seq_range = opts.range.split("..")
# minimal protein length
min_pro_len = int(opts.min_pro_len)
# filename to store translated proteins combined in one
protein_filename = "protein.fasta"
# directory to store SMART results
smart_results_dir = "../Results/SMART_results/"
# part of protein allowed to use for SMART prediction
n = 30000
# threshold for SMART prediction
E_VALUE_THRESH = 10

def parse_smart_domains(file_name):
    in_handle = open(file_name, "rU")
    lines = in_handle.readlines()
    domains = []
    domain_status = False
    domain_type = False
    is_domain = False
    for line in lines:
        line = line.rstrip()
        if len(line) > 1:
            pairs = line.split("=")
            is_domain = True
            if len(pairs) == 2:
                if pairs[0] == "DOMAIN":
                    domain_name = pairs[1]
                elif pairs[0] == "START":
                    domain_start = int(pairs[1])
                elif pairs[0] == "END":
                    domain_end = int(pairs[1])
                elif pairs[0] == "TYPE":
                    if pairs[1] != "PFAM":
                        domain_type = True
                    else:
                        domain_type = False
                elif pairs[0] == "STATUS":
                    if pairs[1] == "visible|OK":
                        domain_status = True
                    #False
                    else: domain_status = True
                else:
                    is_domain = False
        else:
            if is_domain & domain_type & domain_status:
                d = SeqFeature(FeatureLocation(domain_start, domain_end), type="Region")
                d.qualifiers = {'region_name': [domain_name]}
                if domain_name != 'low_complexity_region':
                    domains.append(d)
                is_domain = False
                domain_type = False
                domain_status = False
    in_handle.close()
    return domains

in_handle = open(genomic_filename, "rU")
record = SeqIO.read(in_handle, "fasta")
in_handle.close()

if seq_range[1] != ".":
    seq = record.seq[int(seq_range[0]):int(seq_range[1])+1]
else:
    seq = record.seq

table = 1
i = 0
protein_obj = SeqRecord("")
protein_obj.id = "proteinfull"
protein_obj.description = "translated_from_" + record.id

print "finding domains of minimal length "+str(min_pro_len)+"..."
print "=============================="

for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
    for frame in range(3):
        for pro in nuc[frame:].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                # print "%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame)
                protein_obj.seq += pro

if len(protein_obj.seq) == 0:
    print "No ORFs found in sequence or you have problems with it"
    sys.exit()
elif len(protein_obj.seq) > n:
    print "Translation of found ORFs in selected genomic region is larger than\t", n
    print "Please, narrow it"
    sys.exit()
else:
    out_handle = open(protein_filename, "w")
    SeqIO.write(protein_obj, out_handle, "fasta")
    out_handle.close()

try:
    os.remove(smart_results_dir + protein_obj.id + "_SMART_results.txt")
except OSError:
    pass

print "submitting to SMART..."
print "=============================="

os.system("perl SMART_batch.pl --outputDirectory " + smart_results_dir + " --inputFile " + protein_filename)

print "Aligning SMART domains to initial genomic sequence..."
print "=============================="
domains = parse_smart_domains(smart_results_dir + protein_obj.id + "_SMART_results.txt")

found = False
for d in domains:
    if d.qualifiers['region_name'][0] == domain_name:
        found = True
        pro = SeqRecord(d.extract(protein_obj.seq))
        print domain_name, "\t", d.extract(protein_obj.seq)
        out_handle = open("domain.fasta", "w")
        SeqIO.write(pro, out_handle, "fasta")
        out_handle.close()
        # blast the sequences
        output = NcbitblastnCommandline(query="domain.fasta", subject=genomic_filename, outfmt=5)()[0]
        # parse output
        blast_result_record = NCBIXML.read(StringIO(output)) 
        query_start_list = []
        if len(blast_result_record.alignments) > 0:
            alignment = blast_result_record.alignments[0]
            hsps = []
            for hsp in alignment.hsps:
                # filter garbage
                if hsp.expect < E_VALUE_THRESH:
                    hsps.append(hsp)
                    print hsp
if not found:
    print "No domains found. Try again:)"