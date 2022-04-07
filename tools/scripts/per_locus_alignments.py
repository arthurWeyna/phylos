

import sys, os, re
import shutil
import subprocess
import concurrent.futures


fasta_list=sys.argv[1]
working_dir=sys.argv[2]
align_command=sys.argv[3]
align_threads=sys.argv[4]

def fasta_align(fasta):
    nind=sum([1 for l in open(fasta, "r")  if re.search("^>", l)])
    out=re.sub("\.fasta$", ".aligned.fasta", fasta)
    print(out)
    if nind > 1:
        command=align_command+" "+fasta+" | awk '/>/{head=$0}; !/>/{seq[head]=seq[head]''$0}; END{PROCINFO[\"sorted_in\"]=\"@ind_str_asc\"; for (i in seq){print i; print toupper(seq[i])}}' > "+out
        print(command)
        subprocess.run(command, shell=True)
    else:
        open(out, "x").close()

##Directory prep
if os.path.exists(working_dir):
    shutil.rmtree(working_dir)
os.makedirs(working_dir)

##Read fastas and output in per locus files
fastas=sorted([f.rstrip() for f in open(fasta_list, "r")])
outfiles=[]
for f in fastas:
    for l in open(f, "r"):
        if re.search("^>", l):
            outfile=working_dir+"/"+re.sub("^>", "", l).rstrip()+".fasta"
            if outfile not in outfiles: outfiles.append(outfile)
            with open(outfile, "a") as out:
                out.write(">"+re.sub("^.*\/", "", f)+"\n")
                out.close()
        else:
            with open(outfile, "a") as out:
                out.write(l.rstrip()+"\n")
                out.close()

##Align per locus files
executor = concurrent.futures.ProcessPoolExecutor(2)
futures = [executor.submit(fasta_align, f) for f in outfiles]
concurrent.futures.wait(futures)


