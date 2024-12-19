#Student: Pierre Wensel
#Date: January 13, 2024
#FileName: count_matches_fasta.py
#Assignment: Submission#7
#Instructor: Dr. A. Cordomi

#This program involves a function column_matches_fasta(fasta_file) that computes how many positions (columns) have a match in a .fasta file 
#containing only two sequences. Assume both sequences are aligned and have the same length (with gaps).


import Bio
from Bio import SeqIO
import sys
#from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
#I am also using AlignIO which has similar functionality of SeqIO

def count_matches_fasta(filename):
    for i, record in enumerate(SeqIO.parse(filename, 'fasta')):
        if i == 0:
             reference = record
    """Return Hamming distance between equal-length sequences"""
    if len(reference.seq) != len(record.seq):
        raise ValueError("Undefined for sequences of unequal length")
    count=str(sum(reference.seq == record.seq for reference.seq, record.seq in zip(reference.seq, record.seq)))
    #print("\t".join([reference.id, record.id, str(count)]))
    return count


print("NOW TESTING PROGRAM WITH HARD-CODED INPUT FILE cat_and_rat.fasta")
#Hard-coding and specifying filename here as test docstring
#filename1="cat_and_rat.fasta"
print(count_matches_fasta("cat_and_rat.fasta"))

#Additionally, specifying filename to be provided by user at command line
print("FOR ADDITIONAL TESTING, PLEASE ALSO PROVIDE A FILENAME ARGUMENT AT COMMAND LINE PROMPT WHEN EXECUTING PROGRAM")
filename2=sys.argv[1]
print(count_matches_fasta(filename2))

#The following PSSM-based approach to iterate position specific score matrices (PSSMs) to get count of matches did not work but was included here to reference:
#Position specific score matrices (PSSMs) summarize the alignment information in a different way than a consensus, and may be useful 
#for different tasks. Basically, a PSSM is a count matrix. For each column in the alignment, the number of each alphabet letters is counted
#and totaled. The totals are displayed relative to some representative sequence along the left axis. 
#This sequence may be the consensus sequence, but can also be any sequence in the alignment.

'''alignment = AlignIO.read("cat_and_rat.fasta", 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)
second_seq = alignment[1] 
my_pssm = summary_align.pos_specific_score_matrix(second_seq, chars_to_ignore=["N"])
print(my_pssm)

#Iterating through PSSM matrix did not work:
count=0
#Iterate over the rows of matrix:
for row in my_pssm:
    #iterate over the elements in the row:
    for element in row:
        #if(int(element)==2.0):
            count+=1

print("Total matches counted are " + str(count))

#Iterating through PSSM matrix did not work:
count=0
for i in range(len(my_pssm)):  #Loop for access the row
    for j in range(len(my_pssm[i])): #Loop for accessing the column
        if(my_pssm[i][j]==2.0):
            count+=1
print("Total matches counted are " + str(count))'''

#Also tried numpy and converting PSSM to array did not work:
#import numpy as np
#A = np.squeeze(np.asarray(my_pssm)) 
#shape(A)





