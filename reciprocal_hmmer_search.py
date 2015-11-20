#!/usr/bin/python -W ignore
''' This script predicts detects homologs using reciprocal HMMER search.
    The script generates three layers of prediction:
    a) just reciprocal hit
    b) reciprocal hit with domain prediction, each domain has different importance to regard this protein as ortholog
    c) reciprocal hit with exact domain structure
    The algorithm is the following:
    1) We have known proteins from H. sapiens, C. intestinalis, D. melanogaster and C. elegans. Using these sequences, we build an alignment and hmm. Using this hmm, we search for hits in all proteomes.
    2) We perform reciprocal search to ensure that every protein returns hit used for hmm creation. It could be not the best hit but should have score of the same log e-value as first hit in the list.
    3) We then perform domain matching. We predict all domains in these proteins using Pfam or Prosite models.
    4) Then we draw domains and perform manual check.
        4a) for every protein we have a list of the most important functional domains. To make future work faster, we select matches for this domain into the corresponding file for every proteome;
        4b) inside the matches, we search for other principal features. If all features are found, we report the protein;
        4c) if only part of the features is found, we report the protein;
        4d) in these uncertain cases, we search for additional domains in the adjacent genome regions;
        4e) we report hmm score for the search againt the proteome of every representative organism of every phylum;
        4f) if nothing is found, we report the best reciprocal hit found according to Wall et al., 2003;

'''

from Bio import SeqIO
from Bio import SearchIO
import os, sys
from optparse import OptionParser
import re


def make_hmm(fastafile, modelname):
    ''' Takes the file name with fasta sequences and make hmm, returns 0
    '''
    align_command = "clustalo --outfmt=clu -i " + fastafile + " > " + modelname + ".clustal"
    print "Aligning sequences with ClustalO..."
    os.system(align_command)
    make_hmm_command = "hmmbuild -n " + modelname + " " + modelname + ".hmm " + modelname + ".clustal"
    print "Building hmm..."
    os.system(make_hmm_command)


def hmm_search(modelname, organism, out_filename = None):
    ''' Takes hmm name and performs search through proteomes. The directory with proteomes is ../proteomes
    '''
    if out_filename == None:
        out_filename = modelname + ".hits"
    hmm_search_command = "hmmsearch --noali --tblout " + out_filename + " " + modelname + ".hmm ../proteomes/all/" + organism + ".fasta"
    print "Searching in all proteomes with built hmm ..."
    os.system(hmm_search_command)


def filter_forw_hmm_hits(modelname, forw_inc_bitscore_percentage, unique=True):
    ''' Reads hmmer hits from standart hmmer output.
    Returns hits with bitscore more than provided max bitscore percentage.
    By default returns only hits that from unique loci.
    '''
    gene_regex = re.compile(r'GN=([^=]+)')
    try:
        hits = SearchIO.read(modelname + ".hits", "hmmer3-tab")
    except:
        hits = []
    # select the longest isoform
    filtered_hits = []
    try:
        max_bitscore = hits[0].bitscore
    except:
        max_bitscore = 0
    unique_genes = []
    na_gene_count = 0
    for h in hits:
        # try to guess gene from hit description
        gene = None
        for m in gene_regex.finditer(h.description):
            gene = m.group(1)[:-3]
        if not gene:
            gene = "gene" + str(na_gene_count)
            na_gene_count += 1
        #print h.id, gene, na_gene_count

        if gene not in unique_genes:
            unique_genes.append(gene)
            if h.bitscore > max_bitscore * forw_inc_bitscore_percentage:
                filtered_hits.append(h)
        if manual_mode:
            print unique_genes
            # raw_input("...")
    return filtered_hits


def seq_fetch(full_seq_name, organism=None):
    ''' Fetches the sequence from the proteome file.
    Needs full sequence name in Uniprot format like sp|Q15811|ITSN1_HUMAN
    '''
    seq_id = full_seq_name.split("|")[2]
    if not organism:
        # if organism is not set, try to guess from seq name
        organism = full_seq_name.split("|")[2].split("_")[1]
    seq_filename = seq_id + ".fasta"
    fetch_command = "esl-sfetch ../proteomes/all/" + organism + ".fasta \"" + full_seq_name + "\" > " + seq_filename
    try:
        os.system(fetch_command)
        return seq_filename
    except:
        print "Sequence was not found in the proteome of " + organism
        return None


def make_pfam_prediction(seq_filename, out_filename):
    ''' Runs pfam_scan
    '''
    pfam_command = "./pfam_scan -fasta " + seq_filename + " -dir data > " + out_filename
    # print pfam_command
    try:
        os.system(pfam_command)
        return 0
    except:
        return 1


def make_prosite_prediction(seq_filename, out_filename):
    ''' Runs prosite prediction with ps_scan
    '''
    prosite_command = "./ps_scan.pl --pfscan ./pfscan -d data/custom.dat -o pff " + seq_filename + " > " + out_filename
    # print prosite_command
    try:
        os.system(prosite_command)
        return 0
    except:
        return 1


def make_ncoils_prediction(seq_filename, out_filename):
    ''' Predicts coiled coil region with ncoils.
    Currently results does not correspond to SMART, which also uses ncoils prediction.
    On coils server they use ncoils with -nw parameter which is not described and not working both for 1999 and 2002 ver of ncoils.
    Giving up as for now.
    '''
    ncoils_result_file = "ncoils.fasta"
    ccr_regex = re.compile(r'x+')
    ncoils_command = "ncoils -f -win 28 < " + seq_filename + " > " + ncoils_result_file
    os.system(ncoils_command)
    # parse ncoils results
    ncoils_result_handle = open(ncoils_result_file, "rU")
    ncoils_prediction = SeqIO.read(ncoils_result_handle, "fasta")
    ncoils_result_handle.close()
    out_handle = open(out_filename, "w")
    for match in ccr_regex.finditer(str(ncoils_prediction.seq)):
        out_handle.write("CCR" + " " + str(match.start()) + " " + str(match.end()) + " " + ncoils_prediction.id + "\n")
    out_handle.close()


def make_regex_prediction(seq_filename, out_filename):
    ''' Predicts short motives using regular expressions
    '''
    motifs = {
        'AP2Binding': re.compile(r'(DP[WF]|F.D.F|WV.F|F..F..[FL])'),
        'ClathrinBinding': re.compile(r'(L[FILMV].[FILMV][DE]|L[FILMV].[DE][FILMV])'),
        'NPF': re.compile(r'(NPF)'),
    }
    in_handle = open(seq_filename, "rU")
    records = list(SeqIO.parse(in_handle, "fasta"))
    in_handle.close()
    out_handle = open(out_filename, "w")

    for r in records:
        for motif_name, motif_regex in motifs.iteritems():
            for match in motif_regex.finditer(str(r.seq)):
                out_handle.write(motif_name + " " + str(match.start()) + " " + str(match.end()) + " " + r.id + "\n")
    out_handle.close()


def parse_coord_file(coord_file, type):
    ''' Parses the file with coordinates of the predicted features
    Type could be ncoils, prosite, pfam, regex.
    '''
    features_array = []
    in_handle = open(coord_file, "rU")
    for l in in_handle.readlines():
        # if the line is not commented
        if l[0] != "#" and len(l) > 1:
            splitted = l.split()
            if type == "pfam":
                seq_name_index = 0
                f_start_index = 3
                f_end_index = 4
                f_name_index = 6
            elif type == "prosite":
                f_name_index = 3
                f_start_index = 1
                f_end_index = 2
                seq_name_index = 0
            else:
                # currently no differences with prosite
                f_name_index = 0
                f_start_index = 1
                f_end_index = 2
                seq_name_index = 3
            try:
                seq_name = splitted[seq_name_index].split("|")[2]
            except:
                seq_name = splitted[seq_name_index]
            f_start = splitted[f_start_index]
            f_end = splitted[f_end_index]
            f_name = splitted[f_name_index]
            features_array.append([f_name, f_start, f_end, seq_name])
    in_handle.close()
    return features_array    


def domain_match(seq_filename, domains):
    ''' Performs match of the sequence to the domain set.
    Domains is the dictionary { "SH3": 1, "EF-hand": 1, "PH": 0} where domain has 1 if it is crucial for the protein and less otherwise. Domain names should correspond to the PFAM hmm. Motifs also could be included.
    '''
    pfam_coords_file = "pfam_coords"
    prosite_coords_file = "prosite_coords"
    ncoils_coords_file = "ncoils_coords"
    regex_coords_file = "regex_coords"
    make_pfam_prediction(seq_filename, pfam_coords_file)
    make_prosite_prediction(seq_filename, prosite_coords_file)
    make_ncoils_prediction(seq_filename, ncoils_coords_file)
    make_regex_prediction(seq_filename, regex_coords_file)
    features = []
    features.extend(parse_coord_file(pfam_coords_file, "pfam"))
    features.extend(parse_coord_file(prosite_coords_file, "prosite"))
    features.extend(parse_coord_file(ncoils_coords_file, "ncoils"))
    features.extend(parse_coord_file(regex_coords_file, "regex"))
    f_dict = {}
    for f in features:
        if f[3] not in f_dict.keys():
            f_dict[f[3]] = []
        f_dict[f[3]].append(f[:3])
    final_score = 0
    print f_dict
    print f_dict.keys()

#    raw_input("...")
    for seq_name, f in f_dict.iteritems():
        not_found = ""
        found = ""
        seq_features = [ff[0].split("_")[0] for ff in f]
        seq_features_set = set(seq_features)
        score = 0
        for d in domains:
            if d in seq_features_set:
                score += domains[d]
                found += " " + str(score) + " " + d
            else:
                not_found += " " + d
        if manual_mode:
            print "Domains found in " + seq_filename + ": " + found
            raw_input("Check the results of domain match...\n")
        final_score = score / sum(domains.values())
    return final_score


def reciprocal_hmm_search(modelname, modelname_regex, filename, organism, rev_inc_bitscore_percentage, out_filename = None):
    ''' Performs reciprocal hmmer search against the proteome of given organism
    '''
    print "# Reciprocal search..."
    is_found = False
    if out_filename == None:
        out_filename = modelname + ".rechits_" + organism
    reciprocal_search_command = "phmmer --noali --tblout " + out_filename + " " + filename + " ../proteomes/all/" + organism + ".fasta > hmmer_res"
    os.system(reciprocal_search_command)
    try:
        hits = SearchIO.read(out_filename, "hmmer3-tab")
        max_bitscore = hits[0].bitscore
    except:
        hits = []
        max_bitscore = 0
    if len(hits) > 0:
        if re.search(modelname_regex, hits[0].description):
            is_found = True
    # for h in hits:
    #     if h.bitscore > rev_inc_bitscore_percentage * max_bitscore:
    #         if manual_mode:
    #             print modelname_regex, h.description, re.search(modelname_regex, h.description)
    #             raw_input("...")
    #         if re.search(modelname_regex, h.description):
    #             is_found = True
    if manual_mode:
        print filename, is_found
        raw_input("Check reciprocal search results...")
    return is_found


#domains = { "EF-hand": 1, "SH3": 1, }
manual_mode = False
# domains = { "PID": 1, "NumbF": 1 }
domains = { "Adaptin": 1, }
# domains = { "Adap": 1, "PS51070": 1}
modelname = "AP2A"
# modelname_regex = r'EPN|CLINT|Epsin'
#modelname_regex = r'AGFG'
# modelname_regex = r'LDLRAP'
modelname_regex = r'AP2A'
forw_inc_bitscore_percentage = 0.3
rev_inc_bitscore_percentage = 0.5
# organisms = [
#                 "LEIBR", "LEIIN", "LEIMA", "TRYB2", "TRYCI", "TRYCC", "TRYVY", "THETR", 
#                 "FONAL", "TRIVA", "NAEGR", "GIAIB", "PLAF7", "THEPA", "CRYMR", "EIMAC", 
#                 "EIMMA", "TOXGO", "GRENI", "TETTS", "PARTE", "PERM5", "THEAN", "BABBO", 
#                 "NEOCL", "PHATC", "THAOC", "THAPS", "BLAHO", "HYAAE", "PHYIT", "PHYRM", 
#                 "PHYSP", #"SAPPC", 
#                 "AURAN", "ECTSI", "PYTUL", "EMIHU", "RETFI", "CHLRE", 
#                 "CHLVA", "ARATH", "AUXPR", "SOLTU", "VITVI", "SOLLC", "SOYBN", "GALSU", 
#                 "RICCO", "OSTTA", "WHEAT", "BRAOL", "EUTSA", "ORYSJ", "AMBTC", "SELML", 
#                 "MEDTR", "OSTLU", "BRADI", "VOLCA", "ENTHI", "ENTDS", "DICPU", "DICDI", 
#                 "ENTNP", #"ACACA", 
#                 "POLPA", "EDHAE", "VITCO", "SPRLO", "VAVCU", "BATDJ", 
#                 "ZYMTI", "YARLI", "MAGO7", "GROCL", "OGAPD", "ARTOC", "GIBM7", "FUSO4", 
#                 "BEAB2", "PARBP", "ASPNA", "CHATD", "NECH7", "TALSN", "COLSU", "BOTF1", 
#                 "PENMQ", "PICPG", "CANAL", "YEAST", "MIXOS", "PHACS", "RHOT1", "RHOG2", 
#                 "CRYNJ", "FIBRA", "WALI9", "CRYNH", "USTMA", "MALS4", "GLOTA", "RHIO9", 
#                 "MUCC1", "CAPO3", "SPHAR", #"AMOPA", 
#                 "SALR5", "MONBE", "CLOSI", "SCHMA", 
#                 "TRIAD", "ONCVO", "TRISP", "PRIPA", "LOALO", "TRITR", "CAEEL", "PEDHC", 
#                 "APIME", "DROME", "BOMMO", "CULQU", "SOLIN", "DENPD", "NASVI", "STRPU", 
#                 "AMPQE", "THEKT", "NEMVE", "LOTGI", "CRAGI", "CIOIN", "BRAFL", #"CIOSA", 
#                 "OIKDI", "ANAPL", "HETGA", "LEPOC", "TETNG", "ANOCA", "ORYLA", "ASTMX",
#                 "ORENI", "GASAC", "CRIGR", "CHEMY", "DANRE", "IXOSC", "FICAL", "XENTR", 
#                 "CHICK", "MOUSE", 
#                 "HUMAN",
#             ]

organisms = [
                # "LEIBR", "LEIIN", "LEIMA", "TRYB2", "TRYCI", "TRYCC", "TRYVY", "THETR", 
                # "FONAL", "TRIVA", "NAEGR", 
                "GIAIB", 
                # "PLAF7", "THEPA", "CRYMR", "EIMAC", 
                # "EIMMA", "TOXGO", "GRENI", "TETTS", "PARTE", "PERM5", "THEAN", "BABBO", 
                # "NEOCL", "PHATC", "THAOC", "THAPS", "BLAHO", "HYAAE", "PHYIT", "PHYRM", 
                # "PHYSP", 
                # "AURAN", "ECTSI", "PYTUL", "EMIHU", "RETFI", "CHLRE", "CHLVA", "ARATH", 
                # "AUXPR", "SOLTU", "VITVI", "SOLLC", "SOYBN", "GALSU", "RICCO", "OSTTA", 
                # "WHEAT", "BRAOL", "EUTSA", "ORYSJ", "AMBTC", "SELML", "MEDTR", "OSTLU", 
                # "BRADI", "VOLCA", 
                # "ENTHI", "ENTDS", "DICPU", "DICDI", "ENTNP", "POLPA", "EDHAE", 
                # "VITCO", "SPRLO", "VAVCU", 
                # "BATDJ", "ZYMTI", "YARLI", "MAGO7", "GROCL", "OGAPD", 
                # "ARTOC", "GIBM7", "FUSO4", "BEAB2", "PARBP", "ASPNA", "CHATD", "NECH7", 
                # "TALSN", "COLSU", "BOTF1", "PENMQ", "PICPG", "CANAL", "YEAST", "MIXOS", 
                # "PHACS", "RHOT1", "RHOG2", "CRYNJ", "FIBRA", "WALI9", "CRYNH", "USTMA", 
                # "MALS4", "GLOTA", "RHIO9", "MUCC1", "CAPO3", "SPHAR", "SALR5", "MONBE", 
                # "CLOSI", "SCHMA", "TRIAD", "ONCVO", "TRISP", "PRIPA", "LOALO", "TRITR", 
                # "CAEEL", "PEDHC", "APIME", "DROME", "BOMMO", "CULQU", "SOLIN", "DENPD", 
                # "NASVI", "IXOSC", "STRPU", "AMPQE", "THEKT", "NEMVE", "LOTGI", "CRAGI", 
                # "CIOIN", "BRAFL", "OIKDI", "DANRE", "HETGA", "LEPOC", "TETNG", "ORYLA", 
                # "ASTMX", "ORENI", "GASAC", "CHEMY", "ANOCA", "ANAPL", "FICAL", "XENTR", 
                # "CHICK", "CRIGR", "MOUSE", "HUMAN"
]
result_handle = open(modelname + "_res", "w")
garbage_handle = open(modelname + "_garbage", "w")
for organism in organisms:
    hmm_search(modelname, organism)
    hits = filter_forw_hmm_hits(modelname, forw_inc_bitscore_percentage)
    is_something_found = False
    found_hits = []
    for hit in hits:
        is_hit_found = False
        hit_id = hit.id
        filename = seq_fetch(hit_id)
        score = domain_match(filename, domains)
        if score == 1:
            is_hit_found = reciprocal_hmm_search(modelname, modelname_regex, filename, "HUMAN", rev_inc_bitscore_percentage)
        if is_hit_found:
            found_hits.append(hit)
            is_something_found = True

    if is_something_found:
        garbage_handle.write(organism + "\n")
        for h in found_hits:
            result_handle.write(organism + "\t" + h.id + "\t" + str(h.bitscore) + "\t" + str(h.evalue) + "\n")
    else:
        result_handle.write(organism + "\n")
        if len(hits) > 0:
            garbage_handle.write(organism + "\t" + hits[0].id + "\t" + str(hits[0].bitscore) + "\t" + str(hits[0].evalue) + "\n")
        else:
            garbage_handle.write(organism + "\n")
result_handle.close()
garbage_handle.close()