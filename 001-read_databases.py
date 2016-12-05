#!/usr/bin/python

# USAGE EXAMPLE

from reg_phos_reader import get_kinase_group, get_substrates
from get_windows import get_windows
from fasta_tools import get_relevant_db

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import motifs

import pandas as pd

my_kinases = get_kinase_group("./regPhos/RegPhos_kinase_human.txt",
                              "CMGC")

low_memory=False

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

#TODO make a motif per central aminoacid ??

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
model= "./ModelOrganisms/UP000000625_83333.fasta"
score_lists=cross_score(my_pssm, model, start=1, end=100)
score_lists[0].head()["scores"]
score_lists[0].head()["id"]

# convert to nested data frames

my_data_frame = pd.DataFrame()
my_data_frame['kinase'] = my_kinases
my_data_frame['matches'] = score_lists
