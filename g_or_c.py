#Pierre Wensel
#Course: Bioinformatics
#Filename: g_or_c.py
#Date: 31-10-23
# This code prints, for each character of a sequence, a message telling if the character contains G, or another type of base.

# Define a string variable to store the input from prompted user
dna = input("Please enter your sequence containing letters(Program will by default convert input to upper-case letters):")
#Change lower case letters to upper case letters
dna=str.upper(dna)
#For loop iteration to evaluate and print that base (string character) is or is not a G or C with concatenation:
for i in range(0, len(dna)):
    if (dna[i]=="C"):
        print("base " + str(i) + " is a C")
    elif (dna[i]=="G"):
        print("base " + str(i) + " is a G")
    else:
        print("base " + str(i) + " is not a G or C")