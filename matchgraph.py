import edlib
import matplotlib.pyplot as plt
import numpy as np
import heapq
from Bio.Seq import Seq

def rc(seq):
    return str(Seq(seq).reverse_complement())
def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
    A = ref[ref_st: ref_en]
    a = query[query_st: query_en]
    _a = rc(query[query_st: query_en])
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

def calculate_value(points, ref, query):  
    editdistance = 0
    aligned = 0
    for onetuple in points:
        query_st, query_en, ref_st, ref_en = onetuple[0], onetuple[1], onetuple[2], onetuple[3]
        editdistance += calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)
        aligned += query_en - query_st
    return max(aligned - editdistance, 0)

def dna_match(reference,sample,seglen,topk):
    ans=[]
    for i in range(0,len(sample),seglen):
        ri=min(i+seglen,len(sample))
        scorelist=[]
        scorelist1=[]
        for j in range(0,len(reference),seglen):
            rj=min(j+seglen,len(sample))
            scorelist.append((ri-i-edlib.align(reference[j:rj],sample[i:ri])['editDistance'],j,rj))
            scorelist1.append((ri-i-edlib.align(rc(reference[j:rj]),sample[i:ri])['editDistance'],rj,j))
        scorelist=heapq.nlargest(topk,scorelist)
        scorelist1=heapq.nlargest(topk,scorelist1)

        for j in range(topk):
            color=1/topk*j
            if(scorelist[j][0]>0):
                plt.plot((i,ri),(scorelist[j][1],scorelist[j][2]),color=(1-color,0,0))
            if(scorelist1[j][0]>0):
                plt.plot((i,ri),(scorelist1[j][1],scorelist1[j][2]),color=(0,0,1-color))
        if(scorelist[0][0]>0):
            ans.append((i,ri,scorelist[0][1],scorelist[0][2]))
    print(ans)
    print(calculate_value(ans,reference,sample))
            

with open("reference.txt",'r') as file:
    reference=file.read()
with open("sample.txt",'r') as file:
    sample=file.read()

dna_match(reference,sample,100,20)
plt.show()