#!/usr/bin/python

# USAGE EXAMPLE

from reg_phos_reader import get_kinase_group, get_substrates
from get_windows import get_windows
from fasta_tools import get_relevant_db
from calculate_alignment_scores import cross_score, cross_score_local

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import motifs

import pandas as pd

from ggplot import *

my_kinases = get_kinase_group("./regPhos/RegPhos_kinase_human.txt",
                              "CMGC")
my_substrates = get_substrates("./regPhos/RegPhos_Phos_human.txt",
                               my_kinases)

fasta_db = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                       "fasta",
                       IUPAC.extended_protein)

my_windows = []

for i in (my_substrates['substrates'].tolist()):
    fasta_db = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                           "fasta", IUPAC.extended_protein)
    relevant_db = get_relevant_db(fasta_db, i['AC'])
    my_windows.append(
        get_windows(
            relevant_db,
            i['AC'],
            i['position']))

my_motifs = [[] if len(window['window']) == 0 else
             motifs.create(window['window']) for
             window in my_windows]

# PWM is acronym for position weight matrices

my_pwm = [[] if (len(m) == 0) else m.counts.normalize(pseudocounts=1) for
          m in my_motifs]

# pssm is acronym for position specific storin matrices

my_pssm = [[] if (len(pwm) == 0) else pwm.log_odds() for
           pwm in my_pwm]

# Scoring all elements of a given list
model = "./ModelOrganisms/UP000000625_83333.fasta"
score_lists = cross_score(my_pssm, model, start=1, end=100)
score_lists[0].head()["scores"]
score_lists[0].head()["id"]

# score_list = cross_score(my_pssm, "./ModelOrganisms/UP000000625_83333.fasta")

# convert to nested data frames

my_data_frame = pd.DataFrame()
my_data_frame['kinase'] = [None if isinstance(i, list) else
                           str(i) for i in my_kinases]
my_data_frame['matches'] = [None if isinstance(i, list) else
                            i for i in score_lists]

# Concatenation of dataframes

concat = pd.concat(my_data_frame['matches'].tolist(),
                   keys=my_data_frame['kinase'])
concat.reset_index(level=0, inplace=True)
concat = concat[concat['scores'].notnull()]

ggplot(concat, aes(x='scores', color='kinase')) + geom_density()




# Query aggains the original substrate database

fasta_db1 = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                        "fasta", IUPAC.extended_protein)
relevant_db1 = get_relevant_db(fasta_db1,
                               pd.concat([i['ID'] for
                                          i in my_substrates['substrates']]))

print(len(relevant_db1))
relevant_db1 = dict(zip([i.id for i in relevant_db1],
                        [i for i in relevant_db1]))

score_lists1 = cross_score_local(my_pssm, relevant_db1, end=100)

# convert to nested data frames

my_data_frame1 = pd.DataFrame()
my_data_frame1['kinase'] = [None if isinstance(i, list) else
                            str(i) for i in my_kinases]
my_data_frame1['matches'] = [None if isinstance(i, list) else
                             i for i in score_lists1]

# Concatenation of dataframes

concat1 = pd.concat(my_data_frame1['matches'].tolist(),
                    keys=my_data_frame1['kinase'])
concat1.reset_index(level=0, inplace=True)
concat1 = concat1[concat1['scores'].notnull()]

ggplot(concat1, aes(x='scores', color='kinase')) + geom_density()
