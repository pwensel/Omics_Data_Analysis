#Student Name: Pierre Wensel
#Filename: count_matches.py
#Date: November 6, 2023
#Course: Applications

# Given two sequences with the same length, this code counts the number of matched bases at the same position. 
# and prints a string that has ‘*’ at the matched positions and spaces at the unmatched positions

# Assign a string variable the value of the first input from prompted user 
sequence1=input("Please enter first nucleic sequence (DNA or RNA):")
# Convert all characters in string to upper case:
sequence1=str.upper(sequence1)
# Assign a string variable the value of the second input from prompted user 
sequence2=input("Please enter second nucleic sequence (DNA or RNA) of equal length to the first:")
#Convert all characters in string to upper case:
sequence2=str.upper(sequence2)

# Check to make sure both sequences are of the same length:

if (len(sequence1)==len(sequence2)):
        #Intialize counter of matches
        count=0
        #Intialize the string that will include asterisk match characters
        matches = ""
        # Iterate through the sequences
        for i in range(len(sequence1)):
            # For each iteration, compare the characters
            if (sequence1[i]==sequence2[i]):
                   # Increment match counter if there is a match
                   count+=1
                   # Concatenate an asterisk to the output string at that string index if there is a match
                   matches += "*"
            else:
                   # Concatenate a space to the output string at that string index if there is a match
                   matches += " "
else:
        print("Your 2 sequences are not of equal length. Please try-again")

# Print output
print("")
print (str(count) + " Matches ")
print("")
print (sequence1)
print (sequence2)
print (matches)
print ("")
