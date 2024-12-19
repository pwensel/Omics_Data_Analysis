#Student Name: Pierre Wensel
#Filename: rna_to_prot.py
#Date: November 14, 2023
#Course: Applications

#This code converts a list of DNA or RNA codons into a protein sequence 

# Assign a string variable to store the value of the input from prompted user 
sequence=input("Enter DNA or RNA sequence:").strip()
#Convert lower case letters to upper case:
sequence=str.upper(sequence)

#Initiliazing list
codon_list = []
#Making a comma-delimited list of codon strings
if len(sequence)%3 == 0:
    for i in range(0, len(sequence), 3):
         codon_list.append(sequence[i : i + 3])
else:
    print("Sequence of characters must be divisible by 3. PLease re-enter")
        
#Reads the file and return a dictionary with the genetic code
def read_gencode(inputfile):
    with open(inputfile, "r") as f:
        # Read the contents of the file into a list
        lines = f.readlines()
        # Create an empty dictionary 
        #dict() would also work
        code = {} 
    # Loop through the list of lines 
        for line in lines: 
            key, value = line.split()
        # Store the key-value pairs in the dictionary 
            code[key] = value 

        # Return the dictionary 'code' that now contains the contents of the text file 

        return code

#Takes list of triplets and return a string with the protein sequence.
def rna2prot(codon_list):
    #intiliaze output string
    prot = ""
    
    # Iterate through the three-letter string elements of the list
    if codon_list:
        for codon in codon_list:
    # For each codon, convert to rna. Replace substring in list of strings to convert DNA codons to RNA codons
                codon=codon.replace('T', 'U')
                #Create concatenated string of amino acids using the visible dictionary
                prot = prot + code[codon]
                #Return protein sequence string
    
    return prot
     

#Call method to create globally-visible disctionary before calling rna2prot method
code=read_gencode("gencode.txt")
#Print results and protein sequence
print("The corresponding translated protein sequence is: " + rna2prot(codon_list))