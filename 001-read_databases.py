#!/usr/bin/python

import pandas as pd

kinasedb = pd.read_csv("./regPhos/RegPhos_kinase_human.txt", sep = "\t")
print(kinasedb.head())

ppidb = pd.read_csv("./regPhos/RegPhos_kinase_PPI_human.txt", sep = "\t")
print(ppidb.head())

phosdb = pd.read_csv("./regPhos/RegPhos_Phos_human.txt", sep = "\t")
print(phosdb.head())


phosdb = phosdb.fillna("UNKNOWN")

DMPKfamily = kinasedb[kinasedb['Family'].str.contains("DMPK")]

CMGCfamily = kinasedb[kinasedb['Group'] == "CMGC"]

print (DMPKfamily[DMPKfamily["Kinase"].str.contains("MAP")])

# Generates a list of the kinases

CMGC_list = CMGCfamily["Gene_Symbol"].tolist()
CMGC_list = [str(elem).replace(' ', '') for elem in CMGC_list]
print(CMGC_list)


mapsites = phosdb[phosdb['catalytic kinase'].str.contains("MAP")]

print(mapsites.head())

# match all interactions that have the associated prots

ppidb[ppidb['GENE_a'].isin(CMGC_list)]

# match all sites that have as a catalytic kinase somethin in the list


matching_sites = phosdb[phosdb['catalytic kinase'].isin(CMGC_list)]
print(matching_sites.head())


# DEFINITIONS

# USAGE EXAMPLE

from get_kinase_group import get_kinase_group
from get_substrates import get_substrates
from get_windows import get_windows
from get_relevant_db import get_relevant_db
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Alphabet import IUPAC

my_kinases = get_kinase_group("./regPhos/RegPhos_kinase_human.txt", "CMGC")
my_substrates = get_substrates("./regPhos/RegPhos_Phos_human.txt", my_kinases)

fasta_db = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                       "fasta", IUPAC.extended_protein)

my_windows = []

for i in list(my_substrates.values()):
    fasta_db = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                           "fasta", IUPAC.extended_protein)
    relevant_db = get_relevant_db(fasta_db, i['AC'])
    my_windows.append(
        get_windows(
            relevant_db,
            i['AC'],
            i['position']))

from Bio import motifs
from Bio.Seq import Seq


def forcemotif(sequences):
    if len(sequences) != 0:
        return(motifs.create(sequences))
    else:
        return([])


my_motifs = [forcemotif(window['window']) for
             window in my_windows]

len(my_motifs)

##################################################################

# TODO test if the query by kinase of by gene name give same results
# 1. other file that returns the list of positions

import pandas as pd
import copy as copy
from Bio import SeqIO
from Bio.Alphabet import generic_protein

fasta_db = SeqIO.parse("./ModelOrganisms/UP000005640_9606.fasta",
                      "fasta",
                      generic_protein)

testingsubstrates = list(masubstrates.values())[2]
testingsubstrates = copy.deepcopy(testingsubstrates)

frame = pd.DataFrame()

frame['identifiers'], frame['windows'] = identifiers, windows
frame['residue'] = testingsubstrates['code'].tolist()
frame['position'] = testingsubstrates['position'].tolist()

print(frame)

# some_list = [["AAAA"], ["BBBB"], ["CCCC"], ["abc"]]
# matching = [s[0][1:2] for s in some_list if any(xs in s[0] for xs in ["abc"])]

# 1. script that given a series of sequences returns a series of PWM
# 1.
# 1.
