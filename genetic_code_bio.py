#Student: Pierre Wensel
#Date: November 19, 2023
#File: genetic_code_bio.py
#Course: Applications

#This Python program will accept a sequence of letters from a user, determine if it is dna. rna, or protein, and then employ a user-specified 
#frame number to determine the appropriate codons in dna and rna. It will also print the sequence in dna, rna, and protein alphabet.

#Importing from BioPython:

from Bio.Seq import Seq
from Bio.Data import IUPACData
from Bio.Data import CodonTable
from Bio.Seq import translate

#Assigning a value to back-table disctionary provided by Bio-Python
codons=CodonTable.unambiguous_dna_by_name['Standard'].back_table

'''#Reversing the dictionary value-key pairing:
reverseCodons=dict()
for key in codons:
    val = codons[key]
    reverseCodons[val] = key
print(reverseCodons)'''

#IUPACData.protein_letters_1to3
#IUPACData.ambiguous_dna_complement
#IUPACData.protein_letters
#codons.start_codons
#codons.stop_codons
#codons.forward_table

# Assign a string variable to store the value of the input from prompted user 
sequence=Seq(input("Enter sequence of letters:"))
#Convert lower case letters to upper case:
sequence=sequence.upper()

#Assign integer value for reading frame entered by user:
frame=int(input("Please enter your frame shift number (e.g. integer value of 0,1, or 2)"))

#This function determines if user-entered sequence is DNA, RNA, or protein:
def which_seqtype(sequence):
    rna_dna_sub = 'AGTCU'
    for letter in sequence:
        if letter not in rna_dna_sub:
            return 'PROTEIN'
    if 'U' in sequence:
        return 'RNA'
    else:
        return 'DNA'
    
#This function provides one possible DNA encoding user's PROTEIN alphabet input sequence
def prot2dna(sequence):    
# Iterate through the elements of the Seq sequence
    sequence2=Seq("")
    for amino in sequence:
    # For each match of value to the unique key in the codons dictionary, print element value of dictionary
        sequence2+=codons[amino]
    return sequence2

#This function generates a sequence containing appropriate codons based on use-specified frame number corresponding to start-positoin (0,1, or 2):
def frames(sequence, frame):
    sequence=str(sequence)
    codon_list=list()
    if(frame==0):
        length = len(sequence)
        remainder= length%3
        for i in range (0, (length-remainder),3):
            codon_list.append(sequence[i:i+3])
    if(frame==1):
        length = len(sequence)-1
        remainder= length%3
        for i in range (1, (length-remainder),3):
            codon_list.append(sequence[i:i+3])
    if(frame==2):
        length = len(sequence)-2
        remainder= length%3
        for i in range (2, (length-remainder),3):
            codon_list.append(sequence[i:i+3])
    #List accumulation to generate string from list of codons
    codon_string="".join(codon_list)
    #Convert string to Seq data type
    codon_seq=Seq(codon_string)
    return codon_seq

#The following code will output the input sequence, along with the sequencein dna, rna, and protein alphabets:
print("THE INPUT SEQUENCE YOU ENTERED IS: " + sequence ) 
print("THIS IS " + which_seqtype(sequence))
print()
if (which_seqtype(sequence)== 'PROTEIN'):
    print("YOUR SEQUENCE IN DNA ALPHABET IS:(i.e. One possible DNA encoding your PROTEIN)): " + prot2dna(sequence))
    print("YOUR SEQUENCE IN RNA ALPHABET IS: " + (prot2dna(sequence)).transcribe())
    print("YOUR SEQUENCE IN PROTEIN ALPHABET IS : " + sequence)

if (which_seqtype(sequence) == 'DNA'):
    print("YOUR SEQUENCE IN DNA ALPHABET IS : " + sequence)
    print("YOUR SEQUENCE IN RNA ALPHABET IS: " + sequence.transcribe())
    print("YOUR SEQUENCE IN RNA ALPHABET WITH APPROPRIATE CODONS READ IN SPECIFIED FRAME OF " + str(frame) + " is: " + str(frames(sequence.transcribe(),frame)))
    print("YOUR SEQUENCE IN PROTEIN ALPHABET IN SPECIFIED FRAME OF " + str(frame) + " is: " + str(frames(sequence.transcribe(),frame).translate()))
if (which_seqtype(sequence) == 'RNA'):
    print("YOUR SEQUENCE IN DNA ALPHABET IS : " + str(sequence).replace('U', 'T'))
    print("YOUR SEQUENCE IN RNA ALPHABET IS: " + sequence)
    print("YOUR SEQUENCE IN RNA ALPHABET WITH APPROPRIATE CODONS READ IN SPECIFIED FRAME OF " + str(frame) + " is: " + str(frames(sequence,frame)))
    print("YOUR SEQUENCE IN PROTEIN ALPHABET IN SPECIFIED FRAME OF " + str(frame) + " is: " + str(frames(sequence,frame).translate()))
