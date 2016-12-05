def cross_score(pssms, fasta_database, start = None, end = None):
    # takes a fasta database locations and a list of pssm's and
    # returns list of dataFrames, a data frame of the maximum score
    # and ids database entries for each of the queried pssms
    from Bio import SeqIO
    from Bio.Alphabet import IUPAC
    from calculate_alignment_scores import calculate_alignment_scores
    import pandas as pd
    
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