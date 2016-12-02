#!/usr/bin/python


import pandas as pd

def get_kinase_group(source, group):
    # takes a database file and group name and returns
    # a list of the kinases
    # that match given group.

    kinase_db = pd.read_csv(source, sep="\t")
    kinase_db = kinase_db[kinase_db['Group'].notnull()]

    kinase_list = kinase_db[kinase_db['Group'] == str(group)]
    print(kinase_list.head())
    kinase_list = kinase_list["Kinase"].tolist()
    return(kinase_list)

def get_substrates(db_source, kinase_list):
    # takes a database file and a list of uniprot ID's and
    # returns a nested DataFrame with the kinases as a column
    # and a data frame of data frames.
    import pandas as pd

    substrate_db = pd.read_csv(db_source, sep="\t")
    substrate_db = substrate_db[
        substrate_db['catalytic kinase'].notnull()
    ]

    substrate_list = []

    for kinase in kinase_list:
        substrates = substrate_db[
            substrate_db['catalytic kinase'].str.contains('^'+kinase)
        ]
        substrate_list.append(substrates)

    assert(len(kinase_list) == len(substrate_list)), \
        "Differing length in lists"

    kin_subs = pd.DataFrame()
    kin_subs['kinase'], kin_subs['substrates'] = \
        kinase_list, substrate_list

    return(kin_subs)
