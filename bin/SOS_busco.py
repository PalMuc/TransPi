

############################################## SOS BUSCO #####################################################


#Import required stuff

import pandas as pd
import Bio
from Bio import SeqIO
from functools import reduce
import argparse
import csv
from collections import Counter

#Three arguments are required:

#   1) a table with busco scores. First row has to be the header in the 5 following columns (Busco id  Status  sequence legnth  score)
#   2) a fasta file from which to extract the sequences. The sequence name has to match with the sequence name in --input_busco_table.
#   3) an integer n. This is the number of BUSCO genes that are present in the database used for the analysis.

#   4) OPTIONAL: to be 'rescued' a busco hit has to be - by default - present in at least 3 assemblers outputs.
#      It is possible however to change this value (between 1 and 5) by sepcifying the -min option. 

parser = argparse.ArgumentParser(usage='', description='')
parser.add_argument('-input_busco_table', dest='input_busco_table', required = True)
parser.add_argument('-input_fasta', dest='input_fasta', required = True)
parser.add_argument('-n_genes', dest='genes_number', type=int, required = True)
parser.add_argument('-min', dest='min_assemblers', type=int, required = False)

args = parser.parse_args()



#Stores the number of genes present in the Busco database used for the analysis.

n = args.genes_number

# The busco scores table is imported and converted to a pandas dataframe.

input_database = pd.read_csv(args.input_busco_table, sep='\t', header=0)

# For duplicated entries only the first match is kept. This is done to have each database contain the same number (n) of entries.

unique_database = input_database.groupby((input_database['Busco id'] !=input_database['Busco id'].shift()).cumsum().values).first()


# Split the input busco table into 'assembler-specific' dataframes. From each of these, missing Busco IDs are extracted.
# IDs for Busco genes which are missing is converted to list.

SOAP_database = unique_database[:n]
SOAP_missing = SOAP_database[SOAP_database.Status.eq('Missing')].iloc[:,0].tolist()


SPADES_database = unique_database[n:n*2]
SPADES_missing = SPADES_database[SPADES_database.Status.eq('Missing')].iloc[:,0].tolist()


TransABySS_database = unique_database[n*2:n*3]
TransABySS_missing = TransABySS_database[TransABySS_database.Status.eq('Missing')].iloc[:,0].tolist()


Trinity_database = unique_database[n*3:n*4]
Trinity_missing = Trinity_database[Trinity_database.Status.eq('Missing')].iloc[:,0].tolist()


Velvet_database = unique_database[n*4:n*5]
Velvet_missing = Velvet_database[Velvet_database.Status.eq('Missing')].iloc[:,0].tolist()

TransPi_database = unique_database[n*5:n*6]
TransPi_missing = TransPi_database[TransPi_database.Status.eq('Missing')].iloc[:,0].tolist()


#Missing Busco ID from each assembler output are combined in a single list and duplicates are removed.

all_missing_Busco = SOAP_missing + SPADES_missing + TransABySS_missing + Trinity_missing + Velvet_missing + TransPi_missing
missing_Busco = list(dict.fromkeys(all_missing_Busco))

# For each Busco entry missing in one assembler, the 'status' (e.g., 'complete', 'fragmented',etc) is determined for other assemblers.
# Column name is changed from 'status' to assembler name. This is done to produce  'Busco_comparison_table.csv'

SOAP_final = SOAP_database[SOAP_database['Busco id'].isin(missing_Busco)].iloc[:, 0:2]
SOAP_final.columns = ['SOAP' if x=='Status' else x for x in SOAP_final.columns]

SPADES_final = SPADES_database[SPADES_database['Busco id'].isin(missing_Busco)].iloc[:, 0:2]
SPADES_final.columns = ['SPADES' if x=='Status' else x for x in SPADES_final.columns]

TransABySS_final = TransABySS_database[TransABySS_database['Busco id'].isin(missing_Busco)].iloc[:, 0:2]
TransABySS_final.columns = ['TransABySS' if x=='Status' else x for x in TransABySS_final.columns]

Trinity_final = Trinity_database[Trinity_database['Busco id'].isin(missing_Busco)].iloc[:, 0:2]
Trinity_final.columns = ['Trinity' if x=='Status' else x for x in Trinity_final.columns]

Velvet_final = Velvet_database[Velvet_database['Busco id'].isin(missing_Busco)].iloc[:, 0:2]
Velvet_final.columns = ['Velvet' if x=='Status' else x for x in Velvet_final.columns]

TransPi_final = TransPi_database[TransPi_database['Busco id'].isin(missing_Busco)].iloc[:, 0:2]
TransPi_final.columns = ['TransPi' if x=='Status' else x for x in TransPi_final.columns]


#The output 'Busco_comparison_table' is computed.

dataframes = [SOAP_final, SPADES_final, TransABySS_final, Trinity_final, Velvet_final, TransPi_final]
comparison_table = reduce(lambda left,right: pd.merge(left,right,on='Busco id'), dataframes)

# Determines if Busco entries - missing in Transpi - have a corresponding 'complete' match in other assemblers outputs.
# Matches are converted to list

rescued_from_SOAP = comparison_table[(comparison_table.TransPi =='Missing') & (comparison_table.SOAP == 'Complete')].iloc[:,0].tolist()
rescued_from_SPADES = comparison_table[(comparison_table.TransPi =='Missing') & (comparison_table.SPADES == 'Complete')].iloc[:,0].tolist()
rescued_from_TransABySS = comparison_table[(comparison_table.TransPi =='Missing') & (comparison_table.TransABySS == 'Complete')].iloc[:,0].tolist()
rescued_from_trinity = comparison_table[(comparison_table.TransPi =='Missing') & (comparison_table.Trinity == 'Complete')].iloc[:,0].tolist()
rescued_from_Velvet = comparison_table[(comparison_table.TransPi =='Missing') & (comparison_table.Velvet == 'Complete')].iloc[:,0].tolist()

all_seqs_to_save = []

# If matches are present in each assembler 'list' store (for each match) the Busco id, sequence id, score in list 'all_seqs_to_save'

if len(rescued_from_SOAP) != 0:
    for i in rescued_from_SOAP:
        SOAP_seqs = (i,SOAP_database[SOAP_database['Busco id'] == i].iloc[:,2])
        all_seqs_to_save.append(SOAP_seqs)

if len(rescued_from_SPADES) != 0:
        for i in rescued_from_trinity:
            SPADES_seqs = (i,SPADES_database['Sequence'].loc[SPADES_database['Busco id'] == i].values[0],SPADES_database['Score'].loc[SPADES_database['Busco id'] == i].values[0])
            all_seqs_to_save.append(SPADES_seqs)

if len(rescued_from_TransABySS) != 0:
        for i in rescued_from_TransABySS:
            TransABySS_seqs = (i,TransABySS_database['Sequence'].loc[TransABySS_database['Busco id'] == i].values[0],TransABySS_database['Score'].loc[TransABySS_database['Busco id'] == i].values[0])
            all_seqs_to_save.append(TransABySS_seqs)


if len(rescued_from_trinity) != 0:
        for i in rescued_from_trinity:
            trinity_seqs = (i,Trinity_database['Sequence'].loc[Trinity_database['Busco id'] == i].values[0],Trinity_database['Score'].loc[Trinity_database['Busco id'] == i].values[0])
            all_seqs_to_save.append(trinity_seqs)

if len(rescued_from_Velvet) != 0:
    for i in rescued_from_Velvet:
        Velvet_seqs = (i,Velvet_database['Sequence'].loc[Velvet_database['Busco id'] == i].values[0],Velvet_database['Score'].loc[Velvet_database['Busco id'] == i].values[0])
        all_seqs_to_save.append(Velvet_seqs)


#The next section makes sure that only Busco with matches in at least 3 different assemblers are rescued.
#Matches that satisfy this criteria are appended to the list 'seqs_to_save'

flat_list = [i[0] for i in all_seqs_to_save]
busco_count = Counter(flat_list)

if args.min_assemblers:

    busco_to_save = [k for k, v in busco_count.items() if v >= args.min_assemblers]

else:

    busco_to_save = [k for k, v in busco_count.items() if v >= 3]


seqs_to_save = [item for item in all_seqs_to_save if item[0] in busco_to_save]

#Data is sorted based on the busco score and duplicates are removed. This makes sure that, if a sequence can be retrieved from multiple
#assemblers output, the one with the highest score is considered.

seqs_to_save.sort(key= lambda x: x[2], reverse=True)

checked = set()
unique_seqs_list = []

for busco_id, sequence, score in seqs_to_save:
    if not busco_id in checked:
         checked.add(busco_id)
         unique_seqs_list.append((busco_id,sequence))

#Finally, the IDs for the sequences to rescue are inserted in the list 'sequences_IDs_to_rescue'

sequences_IDs_to_rescue = [ x[1] for x in unique_seqs_list]


#The fasta file is parsed with Biopython SeqIO.parse. And target sequences are extracted.

fasta_to_extract = []

for seqrecord in SeqIO.parse(args.input_fasta, 'fasta'):
    if seqrecord.id in sequences_IDs_to_rescue:
        fasta_to_extract.append(seqrecord)

#Output files are written.

comparison_table.to_csv('Busco_comparison_table.csv',index=False)

with open('sequences_to_add.fasta','w') as outputh:
    SeqIO.write(fasta_to_extract,outputh,'fasta')
