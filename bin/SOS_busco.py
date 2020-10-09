#!/usr/bin/env python3

import pandas as pd
import Bio
from Bio import SeqIO
import argparse
from functools import reduce
import numpy as np
import math
import sys
from collections import Counter

parser = argparse.ArgumentParser(usage='', description='')
parser.add_argument('-input_file_busco', dest='input_file_busco',required=True)
parser.add_argument('-input_file_fasta', dest='input_file_fasta',required=True)
parser.add_argument('-min', dest='min_num_assembler', type=float, required=True)
parser.add_argument('-kmers',dest='kmers',required=True)

args = parser.parse_args()

assemblers_names = ['SOAP','SPADES','TransABySS','Velvet']

all_missing_list = []
list_of_databases = []
final_list = []

Busco_to_save = []

with open(args.input_file_busco) as input_busco_file:

    kmers_list = args.kmers.strip().split(',')
    nr_of_kmers = (len(kmers_list)*4+2)
    column_names = [(assembler + '_' + kmer) for assembler,kmer in zip(assemblers_names,kmers_list) for kmer in kmers_list]
    column_names.insert(3*len(kmers_list) ,'Trinity')
    column_names.insert(len(column_names),'Transpi')
    column_names.insert(0,'Busco ID')

    busco_df = pd.read_csv(input_busco_file, sep=',',header=0,names=['Busco_id','Status','Sequence','Score','Length'])
    busco_unique = busco_df.groupby((busco_df['Busco_id'] !=busco_df['Busco_id'].shift()).cumsum().values).first()

    busco_tables = np.array_split(busco_unique, nr_of_kmers)
    transpi_table = busco_tables[nr_of_kmers-1]

    for table in busco_tables:
        busco_missing = table[table.Status.eq('Missing')].iloc[:,0].tolist()
        all_missing_list.extend(busco_missing)
        missing_Busco = list(dict.fromkeys(all_missing_list))

    for table in busco_tables:
        final_df = table[table['Busco_id'].isin(missing_Busco)].iloc[:, 0:2]
        final_list.append(final_df)

    comparison_table = reduce(lambda left,right: pd.merge(left,right,on='Busco_id'), final_list)
    comparison_table.columns = column_names
    transpi_table = comparison_table[(comparison_table['Transpi'] == 'Missing')]

    comparison_table.to_csv('Complete_comparison_table',sep='\t',index=False)
    transpi_table.to_csv('Transpi_comparison_table',sep='\t',index=False)

    BUSCO_to_rescue =  transpi_table[(transpi_table == 'Complete').any(axis=1)].iloc[:,0].tolist()

    if len(BUSCO_to_rescue) == 0:
        sys.exit(0)
    elif len(BUSCO_to_rescue) != 0:
        for table in busco_tables[:-1]:
            for i in BUSCO_to_rescue:
                seqs = (i,table['Sequence'].loc[table['Busco_id'] == i].values[0],table['Score'].loc[table['Busco_id'] == i].values[0])
                Busco_to_save.append(seqs)

potential_seqs = [t for t in Busco_to_save if not any(isinstance(n, float) and math.isnan(n) for n in t)]
flat_list = [i[0] for i in potential_seqs]
busco_count = Counter(flat_list)

min_number = nr_of_kmers * args.min_num_assembler
busco_to_save = [k for k, v in busco_count.items() if v >= min_number]

seqs_to_save = [item for item in potential_seqs if item[0] in busco_to_save]

seqs_to_save.sort(key= lambda x: x[2], reverse=True)

checked = set()
unique_seqs_list = []

for busco_id, sequence, score in seqs_to_save:
    if not busco_id in checked:
         checked.add(busco_id)
         unique_seqs_list.append((busco_id,sequence))

#The fasta file is parsed with Biopython SeqIO.parse. And target sequences are extracted.
sequences_IDs_to_rescue = [ x[1] for x in unique_seqs_list]
fasta_to_extract = []

for seqrecord in SeqIO.parse(args.input_file_fasta, 'fasta'):
    if seqrecord.id in sequences_IDs_to_rescue:
        fasta_to_extract.append(seqrecord)

#Output files are written.
with open('sequences_to_add.fasta','w') as outputh:
    SeqIO.write(fasta_to_extract,outputh,'fasta')
