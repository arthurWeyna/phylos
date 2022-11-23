

import sys, os, re
import shutil
import subprocess
import concurrent.futures
from itertools import combinations


fis_exec=sys.argv[1]
vcf=sys.argv[2]
threads=int(sys.argv[3])

working_dir=re.sub(".vcf$", ".fis", vcf)
out1=re.sub(".vcf$", ".fis.pergene.txt", vcf)
out2=re.sub(".vcf$", ".fis.total.txt", vcf)

def compute_fis_pair(pair):
    file1=working_dir+"/"+re.sub("^.*\/", "", re.sub(".vcf$", ".fis."+pair[0]+"_"+pair[1], vcf))+".all.txt"
    file2=working_dir+"/"+re.sub("^.*\/", "", re.sub(".vcf$", ".fis."+pair[0]+"_"+pair[1], vcf))+".pergene.txt"
    file3=working_dir+"/"+re.sub("^.*\/", "", re.sub(".vcf$", ".fis."+pair[0]+"_"+pair[1], vcf))+".total.txt"
    command="python "+fis_exec+" "+vcf+" "+" ".join(pair)+" "+file1+" "+file2+" "+file3
    if os.path.exists(file3):
        command="echo 'Pair already exists'"
    subprocess.run(command, shell=True)


##Directory prep
if not os.path.exists(working_dir):
    os.makedirs(working_dir)

##Get pairs of individuals
for l in open(vcf, "r"):
    if re.search("^#CHROM", l):
        l2=l.rstrip().split("\t")
        inds=l2[9:len(l2)]
        pairs=[i for i in combinations(inds, 2) if i[0]!=i[1]]
        break


##Run
executor = concurrent.futures.ProcessPoolExecutor(threads)
futures = [executor.submit(compute_fis_pair, p) for p in pairs]
concurrent.futures.wait(futures)

#Gather results
outfiles=[f for f in os.listdir(working_dir) if re.search("pergene.txt$", f)]
with open(out1, "w") as out:
    for f in outfiles:
        for l in open(working_dir+"/"+f, "r"):
            out.write(l)

outfiles=[f for f in os.listdir(working_dir) if re.search("total.txt$", f)]
with open(out2, "w") as out:
    for f in outfiles:
        for l in open(working_dir+"/"+f, "r"):
            out.write(l)

