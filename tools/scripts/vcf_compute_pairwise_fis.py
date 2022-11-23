#Computes fis for two individuals in a vcf

import sys, re
from collections import Counter

vcf=sys.argv[1]
ind1=sys.argv[2]
ind2=sys.argv[3]
out1=sys.argv[4] #complete output
out2=sys.argv[5] #per gene mean output
out3=sys.argv[6] #overall mean output

def fis2inds(i1a1, i1a2, i2a1, i2a2):
    if i1a1 == "." or i1a2 =="." or i2a1 == "." or i2a2 == ".":
        return "NA"
    tab=dict(Counter([i1a1, i1a2, i2a1, i2a2]).items())
    if len(tab) == 1 or len(tab) == 4:
        return "NA"
    h={a:0 for a in tab.keys()}
    if i1a1 != i1a2:
        h[i1a1]=h[i1a1]+1
        h[i1a2]=h[i1a2]+1
    if i2a1 != i2a2:
        h[i2a1]=h[i2a1]+1
        h[i2a2]=h[i2a2]+1
    pbar=[i/4 for i in tab.values()]
    b=[2*((pbar[i]*(1-pbar[i]))-3*hv/16) for i,hv in enumerate(h.values())]
    c=[i/4 for i in h.values()]
    fis=1-sum(c)/sum(b+c) 
    return fis

with open(out1, "w") as out:
    for n,l in enumerate(open(vcf, "r")):
        if re.search("^#", l):
            if re.search("^#CHROM", l):
                l2=l.rstrip().split("\t")
                posind1=l2.index(ind1)
                posind2=l2.index(ind2)
        else:
            if n % 100000 == 0: 
                print("pair: "+ind1+"/"+ind2+". "+str(n)+" sites done.")
            l2=l.rstrip().split("\t")
            gene=l2[0]
            pos=l2[1]
            if l2[4] == ".":
                fis="NA"
            else:
                alls1=l2[posind1].split(":")[0].split("/")
                alls2=l2[posind2].split(":")[0].split("/")
                fis=fis2inds(alls1[0], alls1[1], alls2[0], alls2[1])
            if fis != "NA":
                out.write(gene+" "+pos+" "+ind1+" "+ind2+" "+str(fis)+"\n")

#mean values
fisdic={} #indiviual values
fisdic2={} #per gene means
for l in open(out1, "r"):
    l2=l.rstrip().split(" ")
    gene=l2[0]
    if gene not in fisdic.keys():
        fisdic[gene]=[]
        fisdic2[gene]=""
    fisdic[gene].append(float(l2[4]))

for g in fisdic.keys():
    fisdic2[g]=sum(fisdic[g])/len(fisdic[g])

tot=[]
with open(out2, "w") as out:
    for g in fisdic2.keys():
        out.write(g+" "+ind1+" "+ind2+" "+str(fisdic2[g])+"\n")
        tot.append(fisdic2[g])

with open(out3, "w") as out:
    out.write(ind1+" "+ind2+" "+str(sum(tot)/len(tot))+"\n")









