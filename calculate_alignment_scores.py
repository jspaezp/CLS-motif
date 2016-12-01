#!/usb/bin/python

# Modified from Bio.motifs.matrix.PositionSpecificScorinMatrix.calculate

# This couple of functions get a position specific scoring matrix
# (pssm) and a sequence and return a series of scores for the align
# at each positions of this pssm.

# It differs from the default implementation in that it dows not have a
# CPython implementation and it does not exclusively work for DNA sequences

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
