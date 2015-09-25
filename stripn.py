#!/usr/bin/python
''' This script strips \n from formatted FASTA file/s.
    If the directory is given, it found all *.fasta files and strips \n from them.
    The script does not modify existing files.
    The script uses awk, so it's fast but works on Linux only.
'''

import os
import sys

try:
    purpose = sys.argv[1]
    if os.path.isfile(purpose):
        os.system("awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }' " + sys.argv[1])
    elif os.path.isdir(purpose):
        for root, dirs, files in os.walk("."):
            if not os.path.isfile(root) and root != ".":
                os.chdir(root)
                for filename in files:
                    splitted = filename.split(".")
                    if splitted[1] == "fasta":
                        # awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' file.fasta
                        os.system("awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }' " + filename + " > " + filename +".fas")
                        print "Processed " + filename
                os.chdir("..")
except:
    print "This script strips \\n from formatted FASTA file/s.\nIf the directory is given, it found all *.fasta files and strips \\n from them.\nIt uses awk so it works on Linux only."