#!/usr/bin/python

# takes a database file and group name and returns
# a list of the kinases
# that match given group.

def get_kinase_group(source, group):
    import pandas as pd

    kinase_db = pd.read_csv(source, sep="\t")
    kinase_db = kinase_db[kinase_db['Group'].notnull()]

    kinase_list = kinase_db[kinase_db['Group'] == str(group)]
    print(kinase_list.head())
    kinase_list = kinase_list["Kinase"].tolist()
    return(kinase_list)
