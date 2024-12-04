import edlib,time

with open("reference.txt",'r') as file:
    reference=file.read()
with open("sample.txt",'r') as file:
    sample=file.read()

st=time.time()
result=edlib.align(reference,sample)
ed=time.time()
print(result['editDistance'],ed-st)