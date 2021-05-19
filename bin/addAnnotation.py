#!/usr/bin/env python3

import os
import glob
import sys
import argparse

parser = argparse.ArgumentParser(usage='xls.py -trinotate_file FILENAME', description='')
parser.add_argument('-trinotateFile', dest='trinotate_file', required=True)
parser.add_argument('-db', dest='db_name', required=False, default="swissprot", help="DB to use for header: uniprot or swissprot")
parser.add_argument('-type', dest='db_type', required=False, default="prot", help="Type of DB to use: nucl or prot")
parser.add_argument('-combine', dest='db_combine', required=False, default="false", help='Use two DBs in headers')
args = parser.parse_args()

swissProtCount=0
uniProtCount=0
for line in open(args.trinotate_file, 'r'):
    line = line.strip()
    lineSplit = line.split("\t")
    if args.db_name == "swissprot" and args.db_type == "nucl":
        if lineSplit[8] != ".":
            print(">" + lineSplit[0] + " SwissProt_Blastx: " + lineSplit[2].split("^")[0])
            uniProtCount += 1
    elif args.db_name == "swissprot" and args.db_type == "prot":
        if lineSplit[6] != ".":
            print(">" + lineSplit[0] + " SwissProt_Blastp: " + lineSplit[6].split("^")[0])
            swissProtCount += 1
    elif args.db_name == "uniprot" and args.db_type == "nucl":
        if lineSplit[8] != ".":
            print(">" + lineSplit[0] + " UniProt_Blastx: " + lineSplit[7].split("^")[0])
            uniProtCount += 1
    elif args.db_name == "uniprot" and args.db_type == "prot":
        if lineSplit[8] != ".":
            print(">" + lineSplit[0] + " UniProt_Blastp: " + lineSplit[8].split("^")[0])
            uniProtCount += 1
