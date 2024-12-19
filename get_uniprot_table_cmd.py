#Name: Pierre Wensel
#Filename: get_uniprot_table_cmd.py
#Applications
#date: December 13, 2023
#Description
#This program includes a function uniprotquery2pandas(query) that sends a query to the UniProt site and returns the data in a panda's DataFrame. 
#As an example of a query, the url for the query "muscarinic+receptor" is this one.

#Hint
#df = pd.DataFrame(data=data_in_rows, columns=column_names)
#where data_in_rows is a list of list (with the actual data) and column_names a list with just the attribute names

'''Tests
    >>> df = uniprotquery2pandas("muscarinic+receptor")
    >>> df.columns[1]
    'Entry Name'
    >>> len(df) > 10
    True
    >>> df.loc['O00622']['Entry Name']
    'CCN1_HUMAN'''

#Entry	Reviewed	Entry Name	Protein names	Gene Names	Organism	Length
#D3Z752	reviewed	SPXN_MOUSE	Spexin (NPQ) (Neuropeptide Q) (Spexin hormone) [Cleaved into: Spexin-1; Spexin-2]	Spx	Mus musculus (Mouse)	116
G5ECB2	reviewed	GABR2_CAEEL	Gamma-aminobutyric acid type B receptor subunit 2	gbb-2 ZK180.1	Caenorhabditis elegans	842
#H2L0Q3	reviewed	GABR1_CAEEL	Gamma-aminobutyric acid type B receptor subunit 1	gbb-1 Y41G9A.4	Caenorhabditis elegans	899'''

# Get data from the Uniprot using requests

import requests
import pandas as pd
import sys

#Accepting argument at command line (Note if file, then argument is 0) 
# Works great!!!
query=sys.argv[1]
#str(query)

def uniprotquery2pandas(query_text):
    """Send a query to Uniprot and put data in a panda's DataFrame"""
    # query_text = "muscarinic+receptor"
    url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv&query=%28{query_text}%29"
    req = requests.get(url)
    lines = req.text.splitlines()
    column_names = lines[0].split("\t")
    data_in_rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(data=data_in_rows, columns=column_names)
    df.set_index("Entry", inplace=True)
    return df

df = uniprotquery2pandas(query)
print(df.head())
