
import sys, re

vcf=sys.argv[1]

if vcf == "-":
    con=sys.stdin
else:
    con=open(vcf, "r")

for l in con:
    if re.search("^#", l):
        print(l.rstrip())
    else:
        l2=l.rstrip().split("\t")
        for ind in range(9, len(l2)):
            gt=l2[ind].split(":")
            dp=[int(n) for n in gt[2].split(",")]
            maj=str(dp.index(max(dp)))
            gt[0]=maj+"/"+maj 
            l2[ind]=":".join(gt)
        print("\t".join(l2))


con.close()
