#!/usr/bin/python

# This is simply a wrapper that returns an empty list if given an empty
# series of sequences and in the other case calls motif.create


def forcemotif(sequences):
    from Bio import motifs
    from Bio.Seq import Seq
    if len(sequences) != 0:
        return(motifs.create(sequences))
    else:
        return([])
