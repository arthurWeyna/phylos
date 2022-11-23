

import sys, re
import copy

vcf=sys.argv[1]

for l in open(vcf, "r"):
    if re.search("^##", l):
        print(l.rstrip())
    elif re.search("^#", l):  
        l2=l.rstrip().split("\t")
        l2[len(l2)-1]=l2[len(l2)-1]+"\t"+l2[len(l2)-1]+"_all1\t"+l2[len(l2)-1]+"_all2"
        print("\t".join(l2))
    else:
        l2=l.rstrip().split("\t")
        gt=l2[len(l2)-1].split(":")
        alls=gt[0].split("/")
        gt1=copy.copy(gt)
        gt2=copy.copy(gt)
        gt1[0]="/".join([alls[0]]*2)
        gt2[0]="/".join([alls[1]]*2)
        print("\t".join(l2)+"\t"+":".join(gt1)+"\t"+":".join(gt2))
        
        



