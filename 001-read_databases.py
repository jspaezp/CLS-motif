#!/usr/bin/python

# USAGE EXAMPLE

from reg_phos_reader import get_kinase_group, get_substrates
from get_windows import get_windows
from fasta_tools import get_relevant_db
from calculate_alignment_scores import calculate_alignment_scores

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import motifs

my_kinases = get_kinase_group("./regPhos/RegPhos_kinase_human.txt", "CMGC")
my_substrates = get_substrates("./regPhos/RegPhos_Phos_human.txt", my_kinases)

fasta_db = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                       "fasta", IUPAC.extended_protein)

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

# Scoring

my_scores = [[] if (len(pssm) == 0) else
             calculate_alignment_scores(pssm, 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA') for
             pssm in my_pssm]

##################################################################

