#!/usb/bin/python

# calculate alignment is based on 
# the one in Bio.motifs.matrix.PositionSpecificScorinMatrix.calculate
# It differs from the default implementation in that it dows not have a
# CPython implementation and it does not exclusively work for DNA sequences

# calculate_alignment_scores uses a position specific scoring matrix
# (pssm) and a sequence and return a series of scores for the align
# at each positions of this pssm.


from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pandas as pd

def _calculate_alignment_scores(pssm, sequence, m, n):
    sequence = sequence.upper()
    scores = []
    for i in range(n - m + 1):
        score = 0.0
        for position in range(m):
            letter = sequence[i + position]
            try:
                score += pssm[letter][position]
            except KeyError:
                score = float("nan")
                break
        scores.append(score)
    return scores


def calculate_alignment_scores(pssm, sequence):
    """Returns the PWM score for a given sequence for all positions.
    Notes:

     - the search is performed only on one strand
     - the result is a one-dimensional list or numpy array
    """
    sequence = str(sequence)
    m = pssm.length
    n = len(sequence)

    scores = _calculate_alignment_scores(pssm, sequence, m, n)
    return scores


def cross_score(pssms, fasta_database, start = None, end = None):
    # takes a fasta database locations and a list of pssm's and
    # returns list of dataFrames, a data frame of the maximum score
    # and ids database entries for each of the queried pssms

    fasta_db = SeqIO.to_dict(
        SeqIO.parse(fasta_database,
                    "fasta",
                    IUPAC.extended_protein))

    score_lists = []
    i = 1

    for pssm in pssms:
        pssm_scores = []
        for keys, values in list(fasta_db.items())[start:end]:
            score = []
            if len(pssm) == 0:
                score = []
            else:
                score = calculate_alignment_scores(pssm, values.seq)
                score = max(score)
            pssm_scores.append(score)
        DF = pd.DataFrame()
        DF['scores'] = pssm_scores
        DF['id'] = [i.id for i in list(fasta_db.values())[start:end]]
        score_lists.append(DF)
        print('motif', i, 'of', len(pssms))
        i += 1

    return(score_lists)
