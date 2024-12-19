#Student: Pierre Wensel
#Date: November 19, 2023
#File: promoters.py
#Course: Applications

#This Python program reads in a fasta DNA file and uses regular expressions functionality to count the number of times a specified DNA motif is p
#is present on both the forward DNA strand and the reverse complimentary DNA strand:


#Importing regular expressions:
import re

#Defining all searched DNA motifs per assignment description:
motif1 = 'TTGATT'
motif2 = 'TA[AGCT]A[AGCT]T'
motif3 = '[AT][AT]TGCTT[TC]A'
motif4 = 'TA[AGTC][AGTC][AGTC]T'
motif5 = '[AT][AT]TGCTT[TC]A[ATGC]{8,10}TG[AGTC]TATAAT'

#The following function reads in a fasta file, skips the first header line, and returns a string representing the genome of T4
#FASTA FILE ept4.fasta corresponding to ACCESSION NUMBER AF158101 was downloaded from ncbi

def read_fasta(filename):
 # Open the file in read mode
    try:
        #with open(file_path, 'r') as dna_file: 
        with open(filename) as dna_file:  
            lines = dna_file.readlines()[1:]
            if len(lines) == 0:
                raise ValueError("The file is empty.")
            # Check if the file is in FASTA format
            #if not lines[0].startswith('>'):
                #raise ValueError("The file is not in FASTA format.")
            #Initializing list
            dna_list = list()
            for line in lines:
                dna_list.append(line.strip())
            # Alternatively, remove the first line and store it as a variable
            #first_line = lines.pop(0).strip()
            dna = "".join(dna_list)
            #Convert all characters of string to uppercase
            dna = str.upper(dna)
            # Return a string
            return dna
 
    except FileNotFoundError:
        raise FileNotFoundError("The specified file path does not exist.")
    
#Calling  read_fasta() t oread in fasta file
dna = read_fasta("epT4.fasta")

print("The following non-DNA nucleotide letters in the original fasta file will be replaced with an A:")
      
for i in range(0, len(dna)):
    if (dna[i] not in ("AGTC")):
        print(dna[i], end=" ")

#Replacing non-DNA letters in dna string:     
dna = re.sub('[^AGTC]', 'A', dna)
 
#print("The Length of forward dna strand=" + str(len(dna)))

#This function first converts the string sequence of DNA characters (obtained from reading fasta file) into its reverse complementary DNA string sequence
#It then finds all matches for motif on forward and rever complementary strands

def find_matches(motif, dna):
# First, the reverse complementary DNA strand is obtained:
    compl_dict = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
    #Using list accumulation 
    compl_dna = "".join(compl_dict[base] for base in dna)
    rev_compl_dna = compl_dna[::-1]

    #Generating matches on forward dna for specified motif
    matches_forw = re.finditer(motif, dna)
    #print("Printing matches for forward strand:")
    #for match in matches_forw:
        #print(match)
        #print(match.start(), match.end(), match.group())

    #Generating matches on reverse complimentary dna steand for specified motif    
    matches_revc = re.finditer(motif, rev_compl_dna)
    #print("Printing matches for reverse complimentary strand:")
    #for match in matches_revc:
        #print(match)
        #print(match.start(), match.end(), match.group())
       
    return matches_forw, matches_revc

#This function returns the total number of instances of motif in the forward and reverse complementary dna strands
def find_nummatches(motif, dna):
    matches_forw, matches_revc = find_matches(motif, dna)
    count_forw = 0
    count_revc= 0
    for match in matches_forw:
        count_forw+=1
    for match in matches_revc:
        count_revc+=1

    return count_forw,count_revc

print()
count_forw,count_revc = find_nummatches(motif1, dna)
print("The number of matches for motif " + motif1 + " in the forward and reverse DNA strands, respectively, are " + str(count_forw) + " and " +  str(count_revc))
print("The total matches for motif " + motif1 + " is " + str(count_forw+count_revc))

print()
count_forw,count_revc = find_nummatches(motif2, dna)
print("The number of matches for motif " + motif2 + " in the forward and reverse DNA strands, respectively, are " + str(count_forw) + " and " +  str(count_revc))
print("The total matches for motif " + motif2 + " is " + str(count_forw+count_revc))

print()
count_forw,count_revc = find_nummatches(motif3, dna)
print("The number of matches for motif " + motif3 + " in the forward and reverse DNA strands, respectively, are " + str(count_forw) + " and " +  str(count_revc))
print("The total matches for motif " + motif3 + " is " + str(count_forw+count_revc))

print()
count_forw,count_revc = find_nummatches(motif4, dna)
print("The number of matches for motif " + motif4 + " in the forward and reverse DNA strands, respectively, are " + str(count_forw) + " and " +  str(count_revc))
print("The total matches for motif " + motif4 + " is " + str(count_forw+count_revc))

print()
count_forw,count_revc = find_nummatches(motif5, dna)
print("The number of matches for motif " + motif5 + " in the forward and reverse DNA strands, respectively, are " + str(count_forw) + " and " +  str(count_revc))
print("The total matches for motif " + motif5 + " is " + str(count_forw+count_revc))






