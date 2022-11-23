#Computes fis for two individuals in a vcf

import sys, re

bed=sys.argv[1]
id=sys.argv[2] #id to output alonside results
out1=sys.argv[3] #per contig mean cov
out2=sys.argv[4] #overall mean cov
min_cov=sys.argv[5] #Only take into account site with more than min_cov coverage.


cntdic={}
for l in open(bed, "r"):
    l2=l.strip().split("\t")
    if l2[0] not in cntdic.keys():
        cntdic[l2[0]]=[0,0]
    cov=int(float(l2[3]))
    if cov >= int(min_cov):
        rang=int(l2[2])-int(l2[1])
        cntdic[l2[0]][0]+=rang
        cntdic[l2[0]][1]+=cov*rang

with open(out1, "w") as out:
    out.write("id chrom len cov\n")
    cov=0
    siz=0
    for c in cntdic.keys():
        if cntdic[c][0]>0:
            ratio=cntdic[c][1]/cntdic[c][0]
        else:
            ratio=0
        out.write(id+" "+c+" "+str(cntdic[c][0])+" "+str(ratio)+"\n")
        cov+=cntdic[c][1]
        siz+=cntdic[c][0]

with open(out2, "w") as out:
    out.write("id len cov\n")
    out.write(id+" "+str(siz)+" "+str(cov/siz)+"\n")
