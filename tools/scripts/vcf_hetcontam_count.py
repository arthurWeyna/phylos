
import sys, re
import scipy.stats

vcf=sys.argv[1]
bed=sys.argv[2]
mincov=int(sys.argv[3])
prob_alt_null=float(sys.argv[4])

firstread=True
contigs={}
covs={}
for l in open(vcf, "r"):
    l2=l.rstrip()
    if re.search("^#", l):
        if re.search("^##contig", l):
            nam=re.sub(",length.*$", "", re.sub("^.*ID=", "", l2))
            contigs[nam]=[re.sub(">$", "", re.sub("^.*,length=", "", l2)), 0]
            covs[nam]=[]
    else:
        if firstread:
            for c in open(bed, "r"):
                c2=c.rstrip().split("\t")
                if int(c2[3]) < mincov:
                    if len(covs[c2[0]]) > 0 and covs[c2[0]][len(covs[c2[0]])-1][1] == int(c2[1]): 
                        covs[c2[0]][len(covs[c2[0]])-1][1]=int(c2[2])
                    else:
                        covs[c2[0]].append([int(c2[1]), int(c2[2])])
            firstread=False
        l3=l2.split("\t") 
        gt=l3[9].split(":")
        if not re.search("\.", gt[0]):
            if sum([ 1 for p in covs[l3[0]] if (int(l3[1]) > p[0]) and (int(l3[1]) <= p[1]) ]) == 0:
                dps=sorted([int(n) for n in gt[2].split(",")], reverse=True)
                probhet=scipy.stats.binom.pmf(dps[0], dps[0]+dps[1], 0.5)
                probhomo=scipy.stats.binom.pmf(dps[1], dps[0]+dps[1], prob_alt_null)
                if probhet > probhomo:
                    contigs[l3[0]][1]=contigs[l3[0]][1]+1 


for i in contigs.keys():
    goodcovlen=int(contigs[i][0])-sum([p[1]-p[0] for p in covs[i]])
    if goodcovlen > 0:
        print(i+" "+str(goodcovlen)+" "+str(contigs[i][1]))
