#Student Name: Pierre Wensel
#Filename: complementary.py
#Date: November 15, 2023
#Course: Applications

#This code converts a string of DNA characters into its complementary DNA string sequence

# Assign a string variable to store the value of the input from prompted user 
sequence=input("Enter DNA sequence:")
#Convert lower case letters to upper case:
sequence=str.upper(sequence)

def comp(sequence):
    complement = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
    #Using list accumulation and returning 
    return "".join(complement[base] for base in sequence)

print(str(comp(sequence)))
