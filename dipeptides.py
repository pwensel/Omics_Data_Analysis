#Pierre Wensel
#Assignment#1
#Filename: dipeptides.py
#Course: Applications
#Date: 31-10-2023

#Write a code that, given a string with residues, generates a list (and prints it) with all possible dipeptides that can be constructed 
#combining the residues in the string. Repetition is possible.

# Define a string variable to store the input from prompted user 
residues = input("Please enter peptide residue sequence containing at least 2 letters corresponding to the 20 amino acids (Input made upper-case by default): ")

#Change any lower case letters to upper case letters
residues=str.upper(residues)
def dipeptide(residues):
    #Intializing an empty list to contain dipeptide elements:
    results=list()
    #Enforcing data validation by ensuring there are at least 2 amino acid residues for a peptide
    if(len(residues)>=2):
        #results=[]
        for i in range(0, len(residues)):
                for j in range(0, len(residues)):
                 #Adding dipeptide elements to list
                 results.append(residues[i]+residues[j])
    else:
        print("Residue sequence too short. Please re-run and re-enter")
    return results
                
#Call function and print list results: 
print(*dipeptide(residues), sep="\n")
