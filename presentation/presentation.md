
# Phosphorylation motif finder

---

## Overview of the project

----

### 1. Given a list of phosphorylation sites

### 2. Get the known phosphorylation sites of a group of kinases

### 3. Know how similar is another set of proteins 

----

![](./regphos.png)

---

## Dependencies

----

#### Biopython

1. SeqIO Module - Fasta parsing
2. IUPAC Module - Implements "alphabets"
3. motifs Module - provides the "Motif" object

----

#### Other dependencies

1. Numpy - Numeric calculations
2. Pandas - provides object "DataFrame" 
3. ggplot2 - provides a nice plotting framework

---

## Project Structure

----

Here goes a super cool image of the dependency-tree

----

#### reg\_phos\_reader.py

Program built from two smaller pandas dependant programs
- get_kinase_group.py(source, group)
--Takes a database file and group name and returns a \nlist of the kinases that match given group.

- get_substrates(db_source, kinase_list)
--Takes a database file and a list of uniprot ID's and
returns a nested DataFrame with the kinases as a column
and a data frame of data frames.
----

#### get_windows.py
Three smaller programs built upon each other
- fill_sequence(sequence, length, fill_right, filler)
-- Fills a string to match a length for instance 'A' for length 5 would be 'AXXXX'

- get_window_strings(entry, position, length)
-- gets a string, a position and a length,
    returns the character at given position,
    and the window up and downstream of given length

- get_windows(database, identifiers, positions, fill, length)
-- Generate a list with only the entries that match any of our IDs, INTO MEMMORY!

----

#### calculate_allignment_scores.py
Three functions used to calculate the aligment scores and cross score of a sequence

- _calculate_alignment_scores(pssm, sequence, m, n)

- calculate_alignment_scores(pssm, sequence)

- cross_score(pssms, fasta_database, start = None, end = None)

----



----

