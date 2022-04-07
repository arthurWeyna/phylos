

import sys, re
import copy

vcf=sys.argv[1]

for l in open(vcf, "r"):
    if re.search("^##", l):
        print(l.rstrip())
    elif re.search("^#", l):  
        l2=l.rstrip().split("\t")
        l2[len(l2)-1]=l2[len(l2)-1]+"\t"+l2[len(l2)-1]+"_mat\t"+l2[len(l2)-1]+"_pat"
        print("\t".join(l2))
    else:
        l2=l.rstrip().split("\t")
        l2[4]=l2[4]+",N"
        nball=len(l2[4].split(","))
        mom_all=list(set(l2[len(l2)-2].split(":")[0].split("/")))
        dau_gt=l2[len(l2)-1].split(":")
        dau_all=list(set(dau_gt[0].split("/")))
        mat_gt=copy.copy(dau_gt)
        pat_gt=copy.copy(dau_gt)
        if len(dau_all) == 2: 
            if len(mom_all) == 1 and mom_all[0] in dau_all:
                mat_gt[0]="/".join([mom_all[0]]*2)
                pat_gt[0]="/".join([[a for a in dau_all if a != mom_all[0]][0]]*2)
            else:
                mat_gt[0]="/".join([str(nball)]*2)
                pat_gt[0]="/".join([str(nball)]*2)
        print("\t".join(l2)+"\t"+":".join(mat_gt)+"\t"+":".join(pat_gt))
        
        



