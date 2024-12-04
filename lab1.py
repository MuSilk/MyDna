from Bio import Align
import numpy as np

def get_score(ref,sam,split1,split2,aligner):
    score=0
    if split1[0]>0 and split2[0]>0:
        score+=aligner.score(ref[0:split1[0]],sam[0:split2[0]])
    if split1[1]>split1[0] and split2[1]>split2[0]:
        score+=aligner.score(ref[split1[0]:split1[1]],sam[split2[0]:split2[1]],"-")
    if split1[1]<len(ref) and split2[1]<len(sam):
        score+=aligner.score(ref[split1[1]:],sam[split2[1]:])
    return score

with open("reference.txt",'r') as file:
    reference=file.read()
with open("sample.txt",'r') as file:
    sample=file.read()

# reference="TATTATAGTCTTCATTCTGTGTATTAGATTACTAAAGCATATTACTTCTGTCTAAATGAAATTT"
# sample="TATTATAGTCTTCATTCTGTATGCTTTAGTAATCTAATACATTACTTCTGTCTAAATGAAATTT"A

aligner=Align.PairwiseAligner()

split1=np.array([0,len(reference)])
split2=np.array([0,len(sample)])
maxscore=aligner.score(reference,sample,"-")
t=max(len(reference),len(sample))

split1=np.array([6500,23200])
split2=np.array([6500,23200])
maxscore=get_score(reference,sample,split1,split2,aligner)
t=100

while(t>=0.01):
    for i in range(20):
        vt=max(int(t),1)
        dv=np.random.randint(-vt,vt+1,size=2)
        nxt_split1=np.clip(dv+split1,0,len(reference))
        while(nxt_split1[1]<nxt_split1[0]):
            dv=np.random.randint(int(-t),int(t+1),size=2)
            nxt_split1=np.clip(dv+split1,0,len(reference))

        dv=np.random.randint(int(-t),int(t+1),size=2)
        nxt_split2=np.clip(dv+split2,0,len(sample))
        while(nxt_split2[1]<nxt_split2[0]):
            dv=np.random.randint(int(-t),int(t+1),size=2)
            nxt_split2=np.clip(dv+split2,0,len(reference))

        new_score=get_score(reference,sample,nxt_split1,nxt_split2,aligner)
        
        delta=new_score-maxscore
        if(np.exp(delta/t)>np.random.rand()):
            split1,split2=nxt_split1,nxt_split2
            maxscore=new_score
            print(maxscore,split1,split2)
    t*=0.97
    print(t)

print(split1,split2,maxscore)