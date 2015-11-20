#!/usr/bin/python
''' This script counts how every protein feature is presented in the phylum.
    The script uses custom organism mapping to phyla.
    The script uses annotated genbank file as the input.
'''

from optparse import OptionParser
from Bio import SeqIO
from ProteinEvolution import eukaryotic_map
import sys
from ProteinEvolution import parse_coord_file

organisms = [
                "LEIBR", "LEIIN", "LEIMA", "TRYB2", "TRYCI", "TRYCC", "TRYVY", "THETR", 
                "FONAL", "TRIVA", "NAEGR", "GIAIB", "PLAF7", "THEPA", "CRYMR", "EIMAC", 
                "EIMMA", "TOXGO", "GRENI", "TETTS", "PARTE", "PERM5", "THEAN", "BABBO", 
                "NEOCL", "PHATC", "THAOC", "THAPS", "BLAHO", "HYAAE", "PHYIT", "PHYRM", 
                "PHYSP", #"SAPPC", 
                "AURAN", "ECTSI", "PYTUL", "EMIHU", "RETFI", "CHLRE", 
                "CHLVA", "ARATH", "AUXPR", "SOLTU", "VITVI", "SOLLC", "SOYBN", "GALSU", 
                "RICCO", "OSTTA", "WHEAT", "BRAOL", "EUTSA", "ORYSJ", "AMBTC", "SELML", 
                "MEDTR", "OSTLU", "BRADI", "VOLCA", "ENTHI", "ENTDS", "DICPU", "DICDI", 
                "ENTNP", #"ACACA", 
                "POLPA", "EDHAE", "VITCO", "SPRLO", "VAVCU", "BATDJ", 
                "ZYMTI", "YARLI", "MAGO7", "GROCL", "OGAPD", "ARTOC", "GIBM7", "FUSO4", 
                "BEAB2", "PARBP", "ASPNA", "CHATD", "NECH7", "TALSN", "COLSU", "BOTF1", 
                "PENMQ", "PICPG", "CANAL", "YEAST", "MIXOS", "PHACS", "RHOT1", "RHOG2", 
                "CRYNJ", "FIBRA", "WALI9", "CRYNH", "USTMA", "MALS4", "GLOTA", "RHIO9", 
                "MUCC1", "CAPO3", "SPHAR", #"AMOPA", 
                "SALR5", "MONBE", "CLOSI", "SCHMA", 
                "TRIAD", "ONCVO", "TRISP", "PRIPA", "LOALO", "TRITR", "CAEEL", "PEDHC", 
                "APIME", "DROME", "BOMMO", "CULQU", "SOLIN", "DENPD", "NASVI", "STRPU", 
                "AMPQE", "THEKT", "NEMVE", "LOTGI", "CRAGI", "CIOIN", "BRAFL", #"CIOSA", 
                "OIKDI", "ANAPL", "HETGA", "LEPOC", "TETNG", "ANOCA", "ORYLA", "ASTMX",
                "ORENI", "GASAC", "CRIGR", "CHEMY", "DANRE", "IXOSC", "FICAL", "XENTR", 
                "CHICK", "MOUSE", 
                "HUMAN",
            ]

def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
        return None
    if len(lst) %2 == 1:
        return lst[((len(lst)+1)/2)-1]
    else:
        return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0


def protein_tail(record, cutoff):
    ''' Determines protein tail. The protein tail starts from the end of the last "big" domain
    '''
    last_feature_end = 0
    for f in record.features:
        if f.location.end > last_feature_end and (f.location.end - f.location.start) >= cutoff:
            last_feature_end = f.location.end + 1
    print ">" + record.name + "/" + str(last_feature_end) + "-" + str(len(record.seq)) + "OS=" + record.annotations['organism']
    print record.seq[last_feature_end - 1:len(record.seq)]
    return len(record.seq[last_feature_end - 1:len(record.seq)])


parser = OptionParser()
parser.add_option("-f", "--file", dest="seq_file", default=None,
                  help="fasta file with protein sequence/s")
parser.add_option("-c", "--coords-file", dest="coords_file", default=None,
                  help="coordinate file with feature/s")
# read options and arguments passed to the script
(opts, args) = parser.parse_args()
seq_file = opts.seq_file
coords_file = opts.coords_file

if coords_file and seq_file:
    print "Select sequence file or coordinate file, but not both"
    sys.exit()

phyla_count = {}
phyla_org_count = {}
for phylum in eukaryotic_map.keys():
    phyla_count[phylum] = {}
    phyla_org_count[phylum] = []
    for o in organisms:
        if o in eukaryotic_map[phylum]:
            phyla_org_count[phylum].append(o)

# if seq_file is set
if seq_file:
    seq_file_handle = open(seq_file, "rU")
    records = list(SeqIO.parse(seq_file_handle, "genbank"))
    seq_file_handle.close()
    for r in records:
        organism = r.annotations['source']
        organism_is_found = False
        # define to which phylum organism belongs
        for phylum, orgs in eukaryotic_map.iteritems():
            if organism in orgs.values():
                organism_is_found = True
                # phyla_org_count[phylum].append(organism)
                for f in r.features:
                    f_name = f.qualifiers['name'][0]
                    if f_name not in phyla_count[phylum].keys():
                        phyla_count[phylum][f_name] = []
                        phyla_count[phylum][f_name].append(r.id)
                    elif r.id not in phyla_count[phylum][f_name]:
                        phyla_count[phylum][f_name].append(r.id)

        if not organism_is_found:
            print "Organism " + organism + " of the protein " + r.id + " was not found in the mapping"

# if coords_file is set
if coords_file:
    features = parse_coord_file(coords_file, "coord")
    for f in features:
        try:
            organism = f[3].split("_")[1]
        except:
            print "Failed to guess organism from " + f[3]
            organism = ""
        organism_is_found = False
        f_name = f[0]
        seq_name = f[3]

        # define to which phylum organism belongs
        for phylum, orgs in eukaryotic_map.iteritems():
            if organism in orgs.keys():
                organism_is_found = True
                if f_name not in phyla_count[phylum].keys():
                    if f_name != "Tail":
                        phyla_count[phylum][f_name] = {}
                    elif f_name == "Tail":
                        phyla_count[phylum][f_name] = []
                if seq_name not in phyla_count[phylum][f_name]:
                    if f_name != "Tail":
                        phyla_count[phylum][f_name][seq_name] = 0
                if f_name != "Tail":
                    phyla_count[phylum][f_name][seq_name] += 1
                else:
                    phyla_count[phylum][f_name].append(int(f[2]) - int(f[1]))
        if not organism_is_found:
            print "Organism " + organism + " of the protein " + seq_name + " was not found in the mapping"

organism_order = ["Euglenozoa", "Excavata", "Chromista", "Archaeplastida", "Amoebozoa", "Fungi", "Choanozoa", "Vertebrata", "Non-vertebrata", "Metazoa",
]
print "\nDomain statistics\n"
for phylum in organism_order:
    l = phyla_count[phylum]
    for m, ll in l.iteritems():
        if m != "Tail":
            org_set = set([p.split("_")[1] for p in ll])
            domain_average = float(len(org_set)) / len(phyla_org_count[phylum])
            domain_occurences = [p for p in ll.values()]
            print "\t".join([phylum, m, "%.2f" % domain_average, str(median(domain_occurences))])

        else:
            average_tail_len = float(sum(ll)) / len(ll)
            print phylum, m, "%.2f" % average_tail_len
