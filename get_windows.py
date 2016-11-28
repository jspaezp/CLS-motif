#!/usr/bin/python

# takes as an input an iterable database (with .id and .seq attributes)
# and an iterable of identifiers, an iterable of positions
# and returns a DataFrame with the sequence windows


def fill_sequence(sequence, length, fill_right=True, filler="X"):
    offset = length - len(sequence)
    if (offset == 0):
        return(sequence)
    if (fill_right):
        return(sequence+offset*"X")
    else:
        return(offset*"X"+sequence)


def get_windows(database, identifiers, positions, fill=True, length=7):
    import pandas as pd
    import copy

    # Generate a list with only the entries
    # that match any of our IDs, INTO MEMMORY!
    # this section is possibly redundant

    relevant_entries = \
        copy.deepcopy(
            [entry for entry in database if
             any(substrate_id in
                 entry.id for substrate_id in identifiers)])

    aminoacid,  upstream_window, downstream_window = [], [], []

    for ident, pos in zip(identifiers, positions):
        for entry in relevant_entries:
            if ('|' + ident in entry.id):
                upstream_window.append(entry.seq[max(pos-1-length, 0):pos-1])
                downstream_window.append(entry.seq[pos:pos+length])
                aminoacid.append(entry.seq[pos-1:pos])

    if fill:
        upstream_window = [fill_sequence(i, length, fill_right=False) for
                           i in upstream_window]
        downstream_window = [fill_sequence(i, length) for
                             i in downstream_window]

    frame = pd.DataFrame()
    frame['aminoacid'] = aminoacid
    frame['upstream'] = upstream_window
    frame['downstream'] = downstream_window
    frame['window'] = frame['upstream']+frame['aminoacid']+frame['downstream']

    return(frame)
