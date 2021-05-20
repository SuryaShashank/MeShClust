import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from bisect import bisect_left
from scipy.spatial import distance


z = input("Please enter file path:")


def reading_fasta_file(file_name):
    count = 0
    # Dictionary
    dict = {}
    ls=[]
    # Finding the initial center.
    for record in SeqIO.parse(file_name, "fasta"): # passing the file
        #print(record.seq)
        ls.append(record.seq)
        count = len(record.seq)
        dict[count] = record.seq
    return ls


# z=file_upload()


# # K-mer details

ls=reading_fasta_file(z)
#print(ls)
quant_values=[]
for i in range(0,len(ls)):
    quant_values.append(str(ls[i]).replace('A','0').replace('C','1').replace('G','2').replace('T','3'))

def cMers(Tsequence, k):
    kFreq = {} # Dictionary to store and count the frequency of each sequence.
    quant_values=[]
    # For loop to read the entire file.
    for i in range(0, len(Tsequence) - k + 1):
        # based on the K value it slices the target sequence.
        kmer = Tsequence[i:i + k]
        if kmer in kFreq:# If sequence already in the dcitionary, add 1 and increment the count.
            kFreq[kmer] += 1
        else:# Add the value to the dictionary and set it's count to 1.
            kFreq[kmer] = 1
    #quant_values
    quant_values=[]
    for i in kFreq.keys():
        quant_values.append(int(str(i).replace('A','0').replace('C','1').replace('G','2').replace('T','3'),base=4))
    #print(quant_values)
    df_KFreq=pd.DataFrame({"K-mer":list(kFreq.keys()),"Count":list(kFreq.values()), "Index":quant_values}) 
    #print(df_KFreq)
    return df_KFreq

Kmerseq=pd.DataFrame({"K-mer":[],"Count":[], "Index":[]}) 
average_SQ = [len(i) for i in ls]
avg_len = sum(average_SQ)/len(average_SQ)
k_count = round(math.log(avg_len,4)-1)
print("K-value is: ",k_count)
for i in range(0,len(ls)):
    Kmers=cMers(str(ls[i]),k_count)
    #print(Kmer)
    Kmers['Seq'] = [ls[i]] * len(Kmers)
    print(Kmers)
    Kmerseq=Kmerseq.append(Kmers,ignore_index = True)
print(Kmerseq)


# # Identity Score

df_seqcomp=pd.DataFrame({"Seq1":[],"Seq2":[], "Identity Score":[]}) 
dict_per= {} # percentage of the calculated score.
dict_sequence= {}
dict_comparison_seq= {}
seqs=0
for i in range(0,len(ls)):
    for j in range(0,len(ls)):
        if ((ls[i]!=ls[j]) and j>=i) and (len(df_seqcomp[df_seqcomp['Identity Score']>70])<50):
            #if j in dict_keys_1 :  
            alignments = pairwise2.align.globalxx(ls[i], ls[j], score_only=True)
            #         print("score :"+ str(alignments)+ "  length :" + str(len(dict[j])))
            seqs+=1
                #print(dict[j])
            dict_per[j] = (round((int(alignments)/int(max(len(ls[i]),len(ls[j])))*100),2))
            dict_sequence[j] = str(ls[j])
            dict_comparison_seq[j] = len(ls[i])
            df_seq=pd.DataFrame({"Seq1":[ls[i]],"Seq2":[ls[j]], "Identity Score":[dict_per[j]]})
            df_seqcomp=df_seqcomp.append(df_seq,ignore_index = True)
        j+=1
        
    i+=1

#df_seq =  pd.DataFrame(list(dict_sequence.items()),columns= ['Sequence_length','sequence'])
#print(list(df_seqcomp.items()))
#print(pd.DataFrame(list(df_seqcomp.items()),columns=['seq1','seq2','Identity Score']))
# print(len(df_seqcomp[df_seqcomp['Identity Score']>70]))
print(df_seqcomp)




# In[12]:
Euclidean_list=[]
for i in range(0,len(df_seqcomp)):
    p1=len(df_seqcomp.Seq1[i])
    p2= len(df_seqcomp.Seq2[i])
    d = distance.euclidean(p1, p2)
    Euclidean_list.append(d)
print("euclidean distance : ", Euclidean_list)



# In[10]:





# def cMers(Tsequence, k):
#     df_KFreq=pd.DataFrame({"K-mer":[],"Count":[], "Index":[]}) 
#     kFreq = {} # Dictionary to store and count the frequency of each sequence.
#     count=0

#     # For loop to read the entire file.
#     for i in range(0, len(Tsequence) - k + 1):
#         # based on the K value it slices the target sequence.
#         kmer = Tsequence[i:i + k]
        
#         if kmer in kFreq:# If sequence already in the dcitionary, add 1 and increment the count.
#             kFreq[kmer] += 1
#         else:# Add the value to the dictionary and set it's count to 1.
#             kFreq[kmer] = 1
        
#     # plt.bar(kFreq.keys(), kFreq.values(), 20, color='g')
#     # plt.show()
#     # plt.close()
#     return kFreq

# tk.Button(master, text='Quit', command=master.quit).grid(row=3, column=0, sticky=tk.W, pady=4)
# tk.Button(master, text='Open', command=file_upload).grid(row=3, column=4, sticky=tk.W, pady=4)


# count = 0
# for i in range(len(readtf)):
#     if '\n' in readtf:
#         count += 1
# targetFile.close()

# # Calculating average sequence length.
# sequence_count = (len(readtf.splitlines()) - len(readtf.split()))+1
# average_SQ = count/sequence_count

# # Calculating the K value based on the average sequence length.
# k_count = round(math.log(average_SQ, 4)-1)
# print(k_count)

# rf = cMers(Tsequence, k_count) # Passing the K value along with the Target Sequence.

# Listrf = [[rf[sequence], sequence] for sequence in rf]
# # Sort the list in descending order.
# Listrf.sort(reverse=True)

# print(Listrf)