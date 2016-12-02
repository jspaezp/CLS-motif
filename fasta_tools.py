#!/usr/bin/python

import copy

def get_relevant_db(parsed_db, identifiers):
    # takes as an input a parser (with .id and .seq attributes)
    # and an iterable set of identifiers
    # and returns a list of entries (in memmory) of the
    # entries that match any of the elements in the set

    relevant_entries = copy.deepcopy(
        [entry for entry in parsed_db if
         any(substrate_id in entry.id for
             substrate_id in identifiers)]
    )

    return(relevant_entries)
