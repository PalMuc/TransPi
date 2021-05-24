#!/usr/bin/env python3

from Bio import SeqIO
import sys

arg1=sys.argv[1]

filename=arg1
for record in SeqIO.parse(filename, "fasta"):
    print(record.id,"\t", len(record.seq)+1)
