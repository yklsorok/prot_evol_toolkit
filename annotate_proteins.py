#!/usr/bin/python
''' This script predicts protein domains and motifs.
    It takes fasta as input and produces GenBank file with features corresponding to the hits.
    You could use Prosite prediction (via ps_scan) or HMMER prediction using locally installed HMMER and pfam_scan.
    UIMs are predicted only if using ps_scan.
    Coiled-coil regions are predicted with locally installed ncoils.
'''

from optparse import OptionParser
from Bio import SeqIO
import os
import re

parser = OptionParser()
parser.add_option("-f", "--file", dest="seq_file", default=None,
                  help="fasta file with protein sequence/s")
parser.add_option("-c", "--coord-file", dest="coord_file", default=None,
                  help="output coordinate file")
parser.add_option("-d", "--pfam-dir", dest="pfam_dir", default=None,
                  help="directory with Pfam-A database")

# read options and arguments passed to the script
(opts, args) = parser.parse_args()
seq_file = opts.seq_file
coord_file = opts.coord_file
pfam_dir = opts.pfam_dir

# compile regular expressions to match motifs in the sequence
# compile regular expression to match coiled coil in the sequence x+
motifs = {
    'AP2_binding': re.compile(r'(DP[WF]|F.D.F|WV.F|F..F..[FL])'),
    'Clathrin_binding': re.compile(r'(L[FILMV].[FILMV][DE]|L[FILMV].[DE][FILMV])'),
    'NPF': re.compile(r'NPF'),
    'CCR': re.compile(r'x+'),
    }

in_handle = open(seq_file, "rU")
records = list(SeqIO.parse(in_handle, "fasta"))
in_handle.close()
out_coord_handle = open(coord_file, "w")
found_features = 0
for r in records:
    # deal with Uniprot sequence id
    if len(r.id.split("|")) > 1:
        accession = r.id.split("|")[1]
        r.id = r.id.split("|")[2]
        r.name = r.id
        r.accession = accession
    print r.id
    helper_fasta_file = "seq.fasta"
    ncoils_result_file = "ncoils_res.fasta"
    # write fasta sequence of the record to the intermediate file
    helper_handle = open(helper_fasta_file, "w")
    SeqIO.write(r, helper_handle, "fasta")
    helper_handle.close()
    # do ncoils prediction
    ncoils_command = "ncoils -f -win 28 -w < " + helper_fasta_file + "> " + ncoils_result_file
    os.system(ncoils_command)
    # parse ncoils results
    ncoils_result_handle = open(ncoils_result_file, "rU")
    ncoils_prediction = SeqIO.read(ncoils_result_handle, "fasta")
    ncoils_result_handle.close()
    for match in motifs["CCR"].finditer(str(ncoils_prediction.seq)):
        out_coord_handle.write("CCR" + " " + str(match.start()) + " " + str(match.end()) + " " + r.id + "\n")
        found_features += 1
    
    helper_fasta_handle = open(helper_fasta_file, "rU")
    helper_fasta_record = SeqIO.read(helper_fasta_handle, "fasta")
    helper_fasta_handle.close()

    # predict UIM
    prosite_result_file = "prosite_res.pff"
    ps_scan_command = "./ps_scan.pl --pfscan ./pfscan -d data/custom.dat -o pff " + helper_fasta_file + "> " + prosite_result_file
    os.system(ps_scan_command)
    prosite_prediction_handle = open(prosite_result_file, "rU")
    for l in prosite_prediction_handle.readlines():
        splitted = l.split()
        out_coord_handle.write("UIM" + " " + splitted[1] + " " + splitted[2] + " " + splitted[0] + "\n")
        found_features += 1
    prosite_prediction_handle.close()
 
    # predict motifs with regular expressions
    for motif_name, regex in motifs.iteritems():
        for match in regex.finditer(str(helper_fasta_record.seq)):
            out_coord_handle.write(motif_name + " " + str(match.start()) + " " + str(match.end()) + " " + r.id + "\n")
            found_features += 1

    # predict domains with pfam_scan
    pfam_scan_command = "./pfam_scan -fasta " + helper_fasta_file + " -dir " + pfam_dir + " > res"
    os.system(pfam_scan_command)
    res_handle = open("res", "rU")
    lines = res_handle.readlines()
    res_handle.close()
    for l in lines:
        if l[0] != "#" and len(l) > 1:
            stripped = " ".join(l.split())
            # print stripped
            splitted = stripped.split(" ")
            out_coord_handle.write(splitted[6] + " " + splitted[3] + " " + splitted[4] + " " + splitted[0] + "\n")
            found_features += 1
out_coord_handle.close()
print "Found " + str(found_features) + " features"