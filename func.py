import matplotlib.pyplot as plt
import edlib 
def show_pairing(points, ref, sam):  
    for onetuple in points:
        sam_st, sam_en, ref_st, ref_en = onetuple[0], onetuple[1], onetuple[2], onetuple[3]

        alignment = edlib.align(ref[ref_st:ref_en], sam[sam_st:sam_en], mode="HW", task="path")
        result=edlib.getNiceAlignment(alignment, ref, sam, gapSymbol='-')
        sam_align=result['query_aligned']
        path=result['matched_aligned']
        ref_align=result['target_aligned']

        # print(sam_align)
        # print(path)
        # print(ref_align)

        cur_sam,cur_ref=sam_st,ref_st
        cur=0
        while(cur<len(path)):
            if path[cur]=='|':
                curr=cur
                while(curr<len(path) and path[curr]=='|'):curr+=1
                plt.plot([cur_sam, cur_sam+curr-cur], [cur_ref, cur_ref+curr-cur], color='black')
                cur_sam+=curr-cur
                cur_ref+=curr-cur
                cur=curr
            elif path[cur]=='.':
                cur_sam+=1
                cur_ref+=1
                cur+=1
            else:
                if(sam_align[cur]=='-'):cur_ref+=1
                else: cur_sam+=1
                cur+=1
    plt.show()

with open("reference.txt",'r') as file:
    reference=file.read()
with open("sample.txt",'r') as file:
    sample=file.read()


# reference="TATTATAGTCTTCATTCTGTGTATTAGATTACTAAAGCATATTACTTCTGTCTAAATGAAATTT"
# sample="TATTATAGTCTTCATTCTGTATGCTTTAGTAATCTAATACATTACTTCTGTCTAAATGAAATTT"
ans=[(0, 11294, 870, 12198), (11294, 14767, 34391, 37884), (14767, 20634, 32342, 38225), (20634, 27263, 32684, 39307), (27263, 30427, 33763, 36981), (30427, 37415, 31436, 38427), (37415, 43487, 32885, 38997), (43487, 55825, 33451, 45878), (55825, 68470, 29234, 41995), (68470, 81586, 19810, 32992), (81586, 90548, 16362, 25418), (90548, 103637, 20259, 33389), (103637, 116064, 16752, 29201), (116064, 129078, 12563, 25645), (129078, 145063, 36740, 52815)]
show_pairing(ans, reference, sample)
