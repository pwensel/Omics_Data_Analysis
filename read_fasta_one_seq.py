#Name: Pierre Wensel
#Filename: read_fasta_one_seq.py
#Applications
#Date: December 13, 2023
#Description

#Description
#This program features a function read_fasta(filename) that is able to read a fasta file with only one sequence that fits in one line (as in here). 
#The function should return a tuple (sequence_id, sequence). In the output, the sequence_id should not contain '>'.

#Tests
#    >>> read_fasta('thefastcat.fasta')
 #   ('fast_cat', 'THEFASTCAT')'''

def read_fasta(filename):
    with open(filename, "r") as f:
        header=f.readline().strip(">\n")
        content=f.readline()
    return (header, content)

#Calling function to test:
print(read_fasta("thefastcat.fasta"))
#This should output the tuple: ('fast_cat', 'THEFASTCAT') and IT DOES!!!