#Student Name: Pierre Wensel
#Filename: amino_3to1.py
#Date: November 6, 2023
#Course: Applications

#This code converts a string of type Ala-Phe-Gly-Ala-Gly-Phe from three letter to one letter amino acid code (AFGAGF) using a dictionary.

#Define dictionary object with unique keys and corresponding single-letter amino acid values:

#Originally, this python program was going to open the amino.txt file, read in the line, and execute the lines of code as follows:

#with open("amino.txt") as f:
    #for line in f:
    #str1=f.read()
    #exec(str1)

#However, calling the exec() function was avoided to minimize undesirable consequences, and the contents of the file were 
#instead copy-pasted here to define the dictionary object containing the unique keys and corresponding amino acid single-letter values:

three2one = dict(
    [
        ["ALA", "A"],
        ["ARG", "R"],
        ["ASP", "D"],
        ["THR", "T"],
        ["PRO", "P"],
        ["HIS", "H"],
        ["SER", "S"],
        ["TRP", "W"],
        ["GLY", "G"],
        ["PHE", "F"],
        ["GLU", "E"],
        ["CYS", "C"],
        ["TYR", "Y"],
        ["VAL", "V"],
        ["LYS", "K"],
        ["GLN", "Q"],
        ["ASN", "N"],
        ["LEU", "L"],
        ["MET", "M"],
        ["ILE", "I"],
    ]
)

# Assign a string variable to store the value of the input from prompted user 
peptide_sequence=input("Enter hyphen-delimited peptide sequence comprised of three-letter names for amino acids (e.g. Ala-Phe-Gly-Ala-Gly-Phe):")
#Convert lower case amino acid names to upper case:
peptide_sequence=str.upper(peptide_sequence)
#Call split() method to convert user-entered string into a comma-delimited and iterable list:
#amino_list=list()
amino_list=peptide_sequence.split("-")

# Iterate through the three-letter elements of the list
for amino in amino_list:
    # For each match of value to the unique key in the dictionary, print element value of dictionary
    print(three2one[amino], end="")

#Print space for new line
print("")

