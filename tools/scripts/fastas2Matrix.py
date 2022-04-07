import sys
import re 
import os

#1: list with all fasta files to concatenate
#2: output
#all fastas are >>>already aligned<<< fastas with the same >labels to concatenate in a matrix filling the blanks
output=sys.argv[2]
fastas=[f.rstrip() for f in open(sys.argv[1], "r")]



inds=[]
fastas_length=[]
alldata_occ={}


alldata={}

for f in fastas:
    alldata[f]={}
    for nlin, line in enumerate(open(f, "r")):
        if (re.search(">", line)):
            ind=line.rstrip()
            alldata[f][ind]=""
            if (ind not in inds):
                inds.append(ind)
        if (not re.search(">", line)):
            alldata[f][ind]=alldata[f][ind]+line.rstrip()

for f in fastas:
    alldata_occ[f]={}
    if (os.stat(f).st_size == 0): seqlen=0
    else: seqlen=len(alldata[f][list(alldata[f].keys())[0]])
    for i in inds:
        alldata_occ[f][i]="1"
        if (i not in alldata[f]):
            alldata_occ[f][i]="0"
            alldata[f][i]="".join(["-"]*seqlen)

partfile=open(output+".part", "w")
partfile.write(str(len(inds))+" lines, "+str(len(fastas))+" columns"+"\n")
maxnamlen=sorted([len(x) for x in inds])[-1]
cnt=1
for f in fastas:
    seqlen=len(alldata[f][list(alldata[f].keys())[1]])
    partfile.write(f+" : "+str(cnt)+"-"+str(cnt+seqlen-1)+"\n")
    cnt=cnt+seqlen
for i in inds:
    partfile.write(i+"".join([" "]*(maxnamlen-len(i)))+" : ")
    for f in fastas:
        partfile.write(alldata_occ[f][i])
    partfile.write("\n")
partfile.close()

outfile=open(output, "w")
for i in inds:
    line=""
    for f in fastas:
        line=line+alldata[f][i]
    outfile.write(i+"\n"+line+"\n")

outfile.close()
print(output)
