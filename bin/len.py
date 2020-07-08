#!/usr/bin/env python3

from Bio import SeqIO
import sys

arg1=sys.argv[1]

filename=arg1
for record in SeqIO.parse(filename, "fasta"):
    print("ID = %s, length %i" % (record.id, len(record.seq)+1))
