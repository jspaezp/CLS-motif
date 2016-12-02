#!/usr/bin/python

# takes as an input an iterable database (with .id and .seq attributes)
# and an iterable of identifiers, an iterable of positions
# and returns a DataFrame with the sequence windows

import pandas as pd
import copy


def fill_sequence(sequence, length, fill_right=True, filler="X"):
    # fills a string to match a length
    # for instance 'A' for length 5 would be 'AXXXX'

    offset = length - len(sequence)
    if (offset == 0):
        return(sequence)
    if (fill_right):
        return(sequence+offset*"X")
    else:
        return(offset*"X"+sequence)


def get_window_strings(entry, position, length):
    # gets a string, a position and a length,
    # returns the character at given position,
    # and the window up and downstream of given length

    aminoacid = entry[position-1:position]
    upstream = entry[max(position-length-1, 0):position-1]
    downstream = entry[position:position+length]
    return aminoacid, upstream, downstream


def get_windows(database, identifiers, positions, fill=True, length=7):

    # Generate a list with only the entries
    # that match any of our IDs, INTO MEMMORY!
    # this section is possibly redundant

    relevant_entries = \
        copy.deepcopy(
            [entry for entry in database if
             any(substrate_id in
                 entry.id for substrate_id in identifiers)])

    aminoacid,  upstream_window, downstream_window = [], [], []

    regions = [aminoacid,  upstream_window, downstream_window]

    for ident, pos in zip(identifiers, positions):
        for entry in relevant_entries:
            [region.append(string_win) for
             region, string_win in
             zip(regions, get_window_strings(entry.seq, pos, length)) if
             ('|' + ident) in entry.id]

    if fill:
        upstream_window = \
            [fill_sequence(i, length, fill_right=False) for
             i in upstream_window]
        downstream_window = \
            [fill_sequence(i, length) for
             i in downstream_window]

    frame = pd.DataFrame()
    frame['aminoacid'] = aminoacid
    frame['upstream'] = upstream_window
    frame['downstream'] = downstream_window
    frame['window'] = frame['upstream'] + \
        frame['aminoacid'] + \
        frame['downstream']

    return(frame)
