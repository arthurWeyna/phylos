#Computes fis for two individuals in a vcf

import sys, re

vcf=sys.argv[1]
out=sys.argv[2] 
discard_multi=True #Only take into account mononucleotide sites?

cntdic={}
cnt=0
for l in open(vcf, "r"):
    cnt+=1
    if cnt % 1000 == 0:
        print(cnt)
    if re.search("^#", l):
        if re.search("^##contig", l):
            cntdic[re.sub("##contig=<ID=", "", re.sub(",length=.*$", "", l.rstrip()))]={}
        if re.search("^#CHROM", l):
            head=l.rstrip().split("\t")
            for i in range(9, len(head)):
                for j in cntdic.keys():
                    cntdic[j][head[i]]=[0,0]
    else:
        l2=l.rstrip().split("\t")
        if discard_multi:
            alls=[l2[3]]+l2[4].split(",")
            if any([len(a)>1 for a in alls]):
                continue
        for i in range(9, len(head)):
            gt=l2[i].split(":")[0].split("/")
            if gt[0] == "." or gt[1] == ".":
                continue
            if gt[0] != gt[1]:
                cntdic[l2[0]][head[i]][0]+=1
            cntdic[l2[0]][head[i]][1]+=1

with open(out, "w") as out:
    out.write("chrom ind het tot\n")
    for c in cntdic.keys():
        for i in range(9, len(head)):
            out.write(c+" "+head[i]+" "+str(cntdic[c][head[i]][0])+" "+str(cntdic[c][head[i]][1])+"\n")
