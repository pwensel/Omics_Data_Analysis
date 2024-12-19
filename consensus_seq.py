#Student: Pierre Wensel
#Date: January 13, 2024
#FileName: consensus_seq.py
#Assignment: Submission#7
#Instructor: Dr. A Cordomi

#This program involves a function consensus_seq that returns the consensus sequence of an alignment in a .fasta file.
#PLEASE PROVIDE A FILENAME ARGUMENT AT COMMAND LINE PROMPT ALSO WHEN EXECUTING

from __future__ import annotations
from pathlib import Path
from itertools import chain
import sys
import Bio
from Bio import SeqIO
#from Bio.Alphabet import IUPAC
from Bio import AlignIO
#I am also using AlignIO which has similar functionality of SeqIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

#This consensus_seq(filename) function generates a consensus sequence from a specified fasta file containing multiple sequence alignment
#Three arguments can be supplied at the command line execution whereby 
#Argument1=python filename
#Argument2=fasta filename (PLEASE PROVIDE A FILENAME ARGUMENT AT COMMAND LINE PROMPT ALSO WHEN EXECUTING PROGRAM)
#Argument3=the percentage of sequences to call a position in the consensus sequence; 
##This is the threshold specifying how common a particular residue has to be at a position before it is added. The default is usually 0.7 (meaning 70% 
# of the sequences to call that position. However, for this assignment, a value of 0 was hardcoded in to obtain the correct test results. 
#THUS, NO NEED TO SUPPLY THIS ARGUMENT AT COMMAND LINE PROMPT WHEN EXECUTING PROGRAM

def consensus_seq(filename):
    alignment = AlignIO.read(filename, 'fasta')
    #common_alignment = MultipleSeqAlignment(chain(*AlignIO.parse(filename, "fasta")))
    summary_align = AlignInfo.SummaryInfo(alignment)
    #summary_align.dumb_consensus(float(sys.argv[2]))
    #The ambiguous character default value is 'N'
    consensus=summary_align.dumb_consensus(0, 'N') #hardcoded values
    return consensus

print("NOW TESTING PROGRAM WITH HARD-CODED INPUT FILE tm3.fa")
#Initially hard-coding and specifying filename here as test docstring
#filename1="tm3.fa"
print(consensus_seq("tm3.fa"))

#Additionally, specifying filename to be provided by user at command line
print("FOR ADDITIONAL TESTING, PLEASE ALSO PROVIDE A FILENAME ARGUMENT AT COMMAND LINE PROMPT WHEN EXECUTING PROGRAM")
filename2=sys.argv[1]
print(consensus_seq(filename2))

