#Student Name: Pierre Wensel
#Assignment: Submission 5
#Filename: kmers.py
#Date: November 29,2023


# This kmers(query, k) function returns a list with all the words of length k contained in the query.

def kmers(query, k):
    word_list=list() #Initializing empty list
    seq_length = len(query)
    if (k > seq_length):
        return word_list
    #remainder= seq_length%k #Calculating remainder to omit last characters from iterative range (This did not work)
    #for i in range (0, (seq_length-(k-remainder))):(This did not work)
    for i in range(len(query)-k+1):
        word_list.append(query[i:i+k]) #Appending list of words
    
    return word_list

print("The following is output to test inputs from assignment statement:")
print(kmers("QLNFQLMSAGQLQ", 3))
#['QLN', 'LNF', 'NFQ', 'FQL', 'QLM', 'LMS', 'MSA', 'SAG', 'AGQ', 'GQL', 'QLQ']
print(kmers("QLNFQLMSAGQLQ", 4))
#['QLNF', 'LNFQ', 'NFQL', 'FQLM', 'QLMS', 'LMSA', 'MSAG', 'SAGQ', 'AGQL', 'GQLQ']

query=input("Please enter new string query to test: ")
k=int(input("Please enter new integer word length to test: "))
print('The following is output to test your string query and word length inputs: ')
print(kmers(query, k))