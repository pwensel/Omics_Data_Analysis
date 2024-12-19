#Student Name: Pierre Wensel
#Assignment: Submission 5
#Filename: naive_match.py
#Date: November 29,2023

#This function naive_match(p, t) returns a list with the offsets of all ocurrences of pattern p in the string t. It doe snot use .find(). 
#This slides the pattern over text one by one and checks for a match (Naive Match Algorithm).

def naive_match(p,t):
    lt=len(t) #length of the string
    lp=len(p) #length of the substring(pattern)
    match_list=[] #Initializing empty list
    if (lt<lp):
        return match_list
    for i in range(lt-lp+1):
        j=0
        while(j<lp):
            if t[i+j]==p[j]:
                j+=1
            else:
                break
        else:
            match_list.append(i) #Appending list with offsets for matches found
    return match_list

print("The following is output to test string and pattern inputs from assignment statement:")
print(naive_match('FAST', 'THEFASTCATTHEFASTRATAFASTRAT'))
print(naive_match('FAST', 'FAST'))
print(naive_match('FASTA', 'FAST'))

t=input("Please enter new string to test: ")
p=input("Please enter new pattern to test: ")
print('The following is output to test your string and pattern inputs: ')
print(naive_match(p,t))