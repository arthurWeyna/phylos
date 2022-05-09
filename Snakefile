
import os 
import numpy as np
import re


##PARAMETERS (feel free to edit)
#Path to directory for data files
DATADIR="/media/bigvol/arthur/phylogenies/phylos/data"
#Path to tools directory
TOOLDIR="/media/bigvol/arthur/phylogenies/phylos/tools"
#Path to run file (group, id)
RUN_FILE="/media/bigvol/arthur/phylogenies/phylos/phylo_list"
MOM_FILE="/media/bigvol/arthur/phylogenies/phylos/mother_daughter_list"

FASTP_THREADS=5
FASTP_OPT="-q 20 -u 70 -n 40 -l 40 -w 1"
MEGAHIT_THREADS=5
MEGAHIT_OPT="--k-min 31 --k-max 101 --k-step 10"
BUSCO_THREADS=5
BUSCO_OPT="-m genome"
BWA_THREADS=10
BWA_OPT="-k 19"
STAN_THREADS=2
STAN_OPT="2 1000 10000"
FREEBAYES_OPT="--min-alternate-count 1 -z 0.05"
VCFTOOLS_FILTER_OPT="--remove-indels --minQ 30 --minDP 5"
BCFTOOLS_CONSENSUS_MINDP=4 #Min depth to report nucleotide in consensus. Positions with less will be "-"
BCFTOOLS_CONSENSUS_MAXGAPPROP=0.8 #remove consensus sequences with more than this proportion of gaps.
BCFTOOLS_MERGE_THREADS=10
HETCOUNTS_CONTAM_OPT="10 0.1" #1: min dp to consider snp; #2: alt allele prob for contam filter
GROUPALIGN_MAFFT_OPT="--auto --adjustdirection"
GROUPALIGN_THREADS=20
TRIMAL_OPT="-gt 0.9 -st 0.8 -w 5"
IQTREE_OPT="-m GTR+I+F+G4 -bb 1000"
IQTREE_THREADS=30

#Construct dic with ids for each group in RUN_FILE
GRPDIC={k:[] for k in [l.rstrip().split(" ")[0] for l in open(RUN_FILE, "r") if not re.search("^#", l)]}
for l in open(RUN_FILE, "r"):
    if not re.search("^#", l):
	    GRPDIC[l.rstrip().split(" ")[0]].append(l.rstrip().split(" ")[1])

#Construct dic with "mom" for each id
MOMDIC={l.rstrip().split(" ")[0]:l.rstrip().split(" ")[1] for l in open(MOM_FILE, "r") if not re.search("^#", l)}

#List of mitochondrial genes of interest
MITOGENES=[re.sub("^>", "", l.rstrip()) for l in open(DATADIR+"/Mscabrinodis_mitoref.fasta", "r") if re.search("^>", l)]
#List of nuclear genes of interest
NUCGENES=[re.sub("^>", "", l.rstrip()) for l in open(DATADIR+"/Y15260_1_buscoref.fasta", "r") if re.search("^>", l)]


##CONSTRAINTS
wildcard_constraints:
    ID="[^/\.]*",
    GRP="[^/\.-]*",
    REF="[^/\.-]*",
    EXT="[^/]*",
    EXT2="[^/-]*",
    SP="[^/\.]*",
    DB="[^/\.]*",
    FILE="[^/]*"


ruleorder:
	fastp_pese_run > fastp_pe_run > fastp_se_run
ruleorder:
	megahit_pese_run > megahit_pe_run > megahit_se_run
ruleorder:
	bwa_pese_run > bwa_pe_run > bwa_se_run
ruleorder:
	bcftools_consensus_run > mother_phase_run
ruleorder:
    get_backup_consensus_sequences > bcftools_consensus_run

#######################################################
#####TARGETS###########################################
#######################################################
#Require final output. Comment out unwanted output.
rule all: 
	input: 
		#[ID]-[REFID][EXTENSION FOR REF FILE]bwa[EXTENSION FOR FASTQ OF ID]freebayes[EXT]
        ###SCANS
		expand("{PATH}/{ID}-Y15260_1_buscoref.bwa.fastp.{MET}.divestim.txt", PATH=DATADIR, ID=GRPDIC["Messor"], MET="angsd contamhet".split(" ")),
		#expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.{MET}.divestim.gathered.txt", PATH=DATADIR, GRP="Messor".split(" "), MET="angsd contamhet".split(" ")),
        ###CONSENSUS SEQUENCES ONLY
		#expand("{PATH}/{ID}-Mscabrinodis_mitoref.bwa.fastp.freebayes.bcfconsensus.fasta", PATH=DATADIR, ID=GRPDIC["Messor"]),
		#expand("{PATH}/{ID}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.fasta", PATH=DATADIR, ID=GRPDIC["Messor"]),
		###PHYLOS
        ##Mito
        #Full
		#expand("{PATH}/{GRP}-Mscabrinodis_mitoref.bwa.fastp.freebayes.bcfconsensus.matrix.trimal.fasta.treefile", PATH=DATADIR, GRP="Messor".split(" ")),
        #Per-gene
		#expand("{PATH}/{GRP}-Mscabrinodis_mitoref.bwa.fastp.freebayes.bcfconsensus.{GENE}_singlegene.trimal.fasta.treefile", PATH=DATADIR, GRP="Messor".split(" "), GENE=MITOGENES),
        ##Nuclear
        #Full
		#expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.matrix.trimal.fasta.treefile", PATH=DATADIR, GRP="Messorphased".split(" ")),
        #Per-gene
		#expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.{GENE}_singlegene.trimal.fasta.treefile", PATH=DATADIR, GRP="Messor".split(" "), GENE=NUCGENES),
        ###SNPS
		#expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfmerged.vcf", PATH=DATADIR, GRP="allstructor".split(" ")),
        
		


#######################################################
#####RULES#############################################
#######################################################

#######################################################
#run fastp on single-end data
rule fastp_se_run:
    input:	
        single="{PATH}/{ID}{EXT}fastq.gz"
    output:
        single=temp("{PATH}/{ID}{EXT}fastp.fastq.gz"),
        json=temp("{PATH}/{ID}{EXT}fastp.json"),
        html=temp("{PATH}/{ID}{EXT}fastp.html")
    threads:
        int(FASTP_THREADS)
    shell:
        "fastp -j {output.json} -h {output.html} {FASTP_OPT} -i {input.single} -o {output.single} -w {threads};"

#######################################################
#run fastp on paired-end data
rule fastp_pe_run:
    input:	
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz"
    output:
        forward=temp("{PATH}/{ID}_1{EXT}fastp.fastq.gz"),
        reverse=temp("{PATH}/{ID}_2{EXT}fastp.fastq.gz"),
        single=temp("{PATH}/{ID}{EXT}fastp.fastq.gz"),
        json=temp("{PATH}/{ID}{EXT}fastp.json"),
        html=temp("{PATH}/{ID}{EXT}fastp.html")
    threads:
        int(FASTP_THREADS)
    shell:
        "fastp -j {output.json} -h {output.html} {FASTP_OPT} --detect_adapter_for_pe -c -i {input.forward} -I {input.reverse} -o {output.forward} -O {output.reverse} --unpaired1 {output.single} --unpaired2 {output.single} -w {threads};"

#######################################################
#run fastp when there is both se and pe data
rule fastp_pese_run:
    input:	
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz",
        single="{PATH}/{ID}{EXT}fastq.gz"
    output:
        forward=temp("{PATH}/{ID}_1{EXT}fastp.fastq.gz"),
        reverse=temp("{PATH}/{ID}_2{EXT}fastp.fastq.gz"),
        singlep=temp("{PATH}/{ID}{EXT}fastp.pe.fastq.gz"),
        singles=temp("{PATH}/{ID}{EXT}fastp.se.fastq.gz"),
        single=temp("{PATH}/{ID}{EXT}fastp.fastq.gz"),
        jsonp=temp("{PATH}/{ID}{EXT}fastp.pe.json"),
        jsons=temp("{PATH}/{ID}{EXT}fastp.se.json"),
        htmlp=temp("{PATH}/{ID}{EXT}fastp.pe.html"),
        htmls=temp("{PATH}/{ID}{EXT}fastp.se.html")
    threads:
        int(FASTP_THREADS)
    shell:
        "fastp -j {output.jsonp} -h {output.htmlp} {FASTP_OPT} --detect_adapter_for_pe -c -i {input.forward} -I {input.reverse} -o {output.forward} -O {output.reverse} --unpaired1 {output.singlep} --unpaired2 {output.singlep} -w {threads};"
        "fastp -j {output.jsons} -h {output.htmls} {FASTP_OPT} -i {input.single} -o {output.singles} -w {threads};"
        #Put se data and orphan reads from pe data in the same file
        "nonemptysingle=`ls -l {output.singles} {output.singlep} | awk '{{if ($5 > 0) print $9}}' | tr '\n' ' '`;"
        "echo $nonemptysingle | xargs zcat | pigz -p{threads} > {output.single}"
		
#######################################################
#run megahit (assembly) on single-end data
rule megahit_se_run:
    input:
        single="{PATH}/{ID}{EXT}fastq.gz"
    output: 
        cont="{PATH}/{ID}{EXT}megahit.fasta",
        log="{PATH}/{ID}{EXT}megahit.log",
        dir=temp(directory("{PATH}/{ID}{EXT}megahit")),
    threads: int(MEGAHIT_THREADS)
    shell:
        "megahit -t {threads} -r {input.single} -o {output.dir} {MEGAHIT_OPT};"
        "mv {output.dir}/log {output.log};"
        """ awk '!/>/{{print}}; />/{{cnt++; print ">contig"cnt}}' {output.dir}/final.contigs.fa > {output.cont}; """

#######################################################
#run megahit (assembly) on paired-end data
rule megahit_pe_run:
    input:
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz",
    output: 
        cont="{PATH}/{ID}{EXT}megahit.fasta",
        log="{PATH}/{ID}{EXT}megahit.log",
        dir=temp(directory("{PATH}/{ID}{EXT}megahit")),
    threads: int(MEGAHIT_THREADS)
    shell:
        "megahit -t {threads} -1 {input.forward} -2 {input.reverse} -o {output.dir} {MEGAHIT_OPT};"
        "mv {output.dir}/log {output.log};"
        """ awk '!/>/{{print}}; />/{{cnt++; print ">contig"cnt}}' {output.dir}/final.contigs.fa > {output.cont}; """


#######################################################
#run megahit (assembly) when there is both se and pe data 
rule megahit_pese_run:
    input:
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz",
        single="{PATH}/{ID}{EXT}fastq.gz"
    output: 
        cont="{PATH}/{ID}{EXT}megahit.fasta",
        log="{PATH}/{ID}{EXT}megahit.log",
        dir=temp(directory("{PATH}/{ID}{EXT}megahit")),
    threads: int(MEGAHIT_THREADS)
    shell:
        "megahit -t {threads} -1 {input.forward} -2 {input.reverse} -r {input.single} -o {output.dir} {MEGAHIT_OPT};"
        "mv {output.dir}/log {output.log};"
        """ awk '!/>/{{print}}; />/{{cnt++; print ">contig"cnt}}' {output.dir}/final.contigs.fa > {output.cont}; """

#######################################################
#download database for busco (should run only once per unique database)
rule busco_database_download:
	output:
		directory("{PATH}/lineages/{DB}_odb{NB}")
	shell:
		"busco --download_path {wildcards.PATH} --download {wildcards.DB}_odb{wildcards.NB};"

#######################################################
#run busco (gene identification and extraction into fasta) using specified database
rule busco_run:
    input:
        cont="{PATH}/{ID}{EXT}fasta",
        db=expand("{DATA}/busco_database/lineages/{{DB}}", DATA=DATADIR)
    output:
        dir=temp(directory("{PATH}/{ID}{EXT}busco-{DB}-{SP}")),
        fasta="{PATH}/{ID}{EXT}busco-{DB}-{SP}.fasta",
        log="{PATH}/{ID}{EXT}busco-{DB}-{SP}.log",
    threads:
        int(BUSCO_THREADS)
    shell:
        "busco -f --download_path {wildcards.PATH} -i {input.cont} --out_path {wildcards.PATH} -o {output.dir} -l {wildcards.DB} -c {threads} --augustus --augustus_species {wildcards.SP} {BUSCO_OPT};"
        "cp {output.dir}/logs/busco.log {output.log};"
        """ awk 'FNR==1{{file=FILENAME; sub(".fna", "", file); sub("^.*\\/", "", file)}}; FNR!=1{{seq[file]=seq[file]""$0}}; END{{PROCINFO["sorted_in"]="@ind_str_asc"; for(f in seq){{print ">busco_"f; print seq[f]}}}}' {output.dir}/run_{wildcards.DB}/busco_sequences/single_copy_busco_sequences/*fna > {output.fasta}; """

#######################################################
#index fasta file with bwa
rule bwa_indexing:
    input:
        cont="{PATH}.fasta"
    output:
        amb=temp("{PATH}.fasta.amb"),
        ann=temp("{PATH}.fasta.ann"),
        pac=temp("{PATH}.fasta.pac"),
        bwt=temp("{PATH}.fasta.bwt.2bit.64"),
        num=temp("{PATH}.fasta.0123"),
    shell:
        "bwa-mem2 index {input.cont}"	

#######################################################
#align single-end reads using bwa 
rule bwa_se_run:
    input:
        index=expand("{{PATH}}/{{REF}}{{EXT}}fasta.{EXT2}", EXT2="amb ann pac bwt.2bit.64 0123".split(" ")),
        single="{PATH}/{ID}{EXT2}fastq.gz",
        cont="{PATH}/{REF}{EXT}fasta",
    output:
        bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
    threads: 
        int(BWA_THREADS)
    shell:
        "bwa-mem2 mem {BWA_OPT} -R '@RG\tID:{wildcards.ID}\tSM:{wildcards.ID}' -t {threads} {input.cont} {input.single} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bam}; """

#######################################################
#align paired-end reads using bwa 
rule bwa_pe_run:
	input:
		index=expand("{{PATH}}/{{REF}}{{EXT}}fasta.{EXT2}", EXT2="amb ann pac bwt.2bit.64 0123".split(" ")),
		forward="{PATH}/{ID}_1{EXT2}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT2}fastq.gz",
		cont="{PATH}/{REF}{EXT}fasta",
	output:
		bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
	threads: 
		int(BWA_THREADS)
	shell:
		" bwa-mem2 mem {BWA_OPT} -R '@RG\tID:{wildcards.ID}\tSM:{wildcards.ID}' -t {threads} {input.cont} {input.forward} {input.reverse} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bam}; "

#######################################################
#align reads using bwa when there is both se and pe data
rule bwa_pese_run:
	input:
		index=expand("{{PATH}}/{{REF}}{{EXT}}fasta.{EXT2}", EXT2="amb ann pac bwt.2bit.64 0123".split(" ")),
		forward="{PATH}/{ID}_1{EXT2}fastq.gz",
		reverse="{PATH}/{ID}_2{EXT2}fastq.gz",
		single="{PATH}/{ID}{EXT2}fastq.gz",
		cont="{PATH}/{REF}{EXT}fasta",
	output:
		bamse=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}se.bam"),
		bampe=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}pe.bam"),
		bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
	threads: 
		int(BWA_THREADS)
	shell:
		"bwa-mem2 mem {BWA_OPT} -R '@RG\tID:{wildcards.ID}\tSM:{wildcards.ID}' -t {threads} {input.cont} {input.single} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bamse}; """
		" bwa-mem2 mem {BWA_OPT} -R '@RG\tID:{wildcards.ID}\tSM:{wildcards.ID}' -t {threads} {input.cont} {input.forward} {input.reverse} | samtools view -bu -F 260 -@ {threads} | samtools sort -m 1G -@ {threads} > {output.bampe}; """
		#merge pe and se alignments
		"samtools merge -@ {threads} -o {output.bam} {output.bamse} {output.bampe};"

#######################################################
#index bam file with samtools
rule samtools_bam_index:
	input:
		"{PATH}.bam"
	output:
		temp("{PATH}.bam.bai")
	shell:
		"samtools index {input};"

#######################################################
#index fasta file with samtools
rule samtools_fasta_index:
	input:
		"{PATH}.fasta"
	output:
		temp("{PATH}.fasta.fai")
	shell:
		"samtools faidx {input};"

#######################################################
#call snp using freebayes
rule freebayes_run:
	input:
		bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
		ref="{PATH}/{REF}{EXT}fasta",
		refindex="{PATH}/{REF}{EXT}fasta.fai",
	output:
		vcf="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf",
	shell:
		"freebayes -f {input.ref} {FREEBAYES_OPT} {input.bam} | bcftools norm -f {input.ref} - > {output}"

#######################################################
#run bedtools to compute coverage
rule bedtools_genomecov_run:
	input:
		bam="{PATH}.bam",
	output:
		bed="{PATH}.coverage.bed",
	shell:
		"bedtools genomecov -bga -max 100 -ibam {input.bam} > {output.bed};"
		
#######################################################
#Consensus sequence with bcftools
rule bcftools_consensus_run:
	input:
		covbed="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed",
		vcf="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf",
		ref="{PATH}/{REF}{EXT}fasta",			
	output:
		covbed=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed.gz"),
		vcfgz=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.gtmajor.vcf.gz"),
		index=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.gtmajor.vcf.gz.csi"),
		fasta="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.bcfconsensus.fasta",
	shell:
		"awk -v min={BCFTOOLS_CONSENSUS_MINDP} '$4<min' {input.covbed} | bgzip > {output.covbed};"
		"python {TOOLDIR}/scripts/vcf_setGT_to_major.py {input.vcf} | bgzip > {output.vcfgz};" 
		"bcftools index -f {output.vcfgz};"
		""" bcftools consensus -s {wildcards.ID} --mask {output.covbed} --mask-with "-" -f {input.ref} {output.vcfgz} | awk -v maxgap={BCFTOOLS_CONSENSUS_MAXGAPPROP} '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]""$0}}; END{{for(i in seq){{l=seq[i]; if((gsub("-", "", l)/length(seq[i]))<maxgap){{print i; print seq[i]}}}}}}' > {output.fasta}; """


#######################################################
#Make use of back up sequences in consensus_sequences directory instead or running the whole thing again
rule get_backup_consensus_sequences:
    input:
        "{PATH}/consensus_sequences/{FILE}.bcfconsensus.backup.fasta"
    output:
        "{PATH}/{FILE}.bcfconsensus.fasta"
    shell:
        "cp {input} {output};"
        
#######################################################
#Phase daughter's vcf according to mother's
rule mother_phase_run:
	input:
		vcfgz="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz",
		vcfindex="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz.csi",
		momvcfgz=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz", PATH=wildcards.PATH, ID=MOMDIC[wildcards.ID], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
		momvcfindex=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz.csi", PATH=wildcards.PATH, ID=MOMDIC[wildcards.ID], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
		covbed="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed",
		momcovbed=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed", PATH=wildcards.PATH, ID=MOMDIC[wildcards.ID], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
		ref="{PATH}/{REF}{EXT}fasta",			
	output:
		badcovbed=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}badcoverage.bed"),
		mombadcovbed=temp("{PATH}/{ID}-mom-{REF}{EXT}bwa{EXT2}badcoverage.bed"),
		totbadcovbed=temp("{PATH}/{ID}-tot-{REF}{EXT}bwa{EXT2}badcoverage.bed"),
		totvcf=temp("{PATH}/{ID}-tot-{REF}{EXT}bwa{EXT2}freebayes.vcf"),
		totvcfpha=temp("{PATH}/{ID}-totphased-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz"),
		totvcfphaindex=temp("{PATH}/{ID}-totphased-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz.csi"),
		mat="{PATH}/{ID}-mat-{REF}{EXT}bwa{EXT2}freebayes.bcfconsensus.fasta",
		pat="{PATH}/{ID}-pat-{REF}{EXT}bwa{EXT2}freebayes.bcfconsensus.fasta",
	shell:
		"awk -v min={BCFTOOLS_CONSENSUS_MINDP} '$4<min' {input.covbed}  > {output.badcovbed};"
		"awk -v min={BCFTOOLS_CONSENSUS_MINDP} '$4<min' {input.momcovbed}  > {output.mombadcovbed};"
		"cat {output.badcovbed} {output.mombadcovbed} | sort -k1,1 -k2,2n | bedtools merge > {output.totbadcovbed};"
		"bcftools merge -0 -m all {input.momvcfgz} {input.vcfgz} | bcftools norm -f {input.ref} > {output.totvcf};"
		"python {TOOLDIR}/scripts/vcf_mom_daughter_phase.py {output.totvcf} | bgzip > {output.totvcfpha};"
		"bcftools index -f {output.totvcfpha};"
		""" bcftools consensus -s {wildcards.ID}_mat --mask {output.totbadcovbed} --mask-with "-" -f {input.ref} {output.totvcfpha} | awk -v maxgap={BCFTOOLS_CONSENSUS_MAXGAPPROP} '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]""$0}}; END{{for(i in seq){{l=seq[i]; if((gsub("-", "", l)/length(seq[i]))<maxgap){{print i; print seq[i]}}}}}}' > {output.mat}; """
		""" bcftools consensus -s {wildcards.ID}_pat --mask {output.totbadcovbed} --mask-with "-" -f {input.ref} {output.totvcfpha} | awk -v maxgap={BCFTOOLS_CONSENSUS_MAXGAPPROP} '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]""$0}}; END{{for(i in seq){{l=seq[i]; if((gsub("-", "", l)/length(seq[i]))<maxgap){{print i; print seq[i]}}}}}}' > {output.pat}; """

#######################################################
#Compress vcf
rule vcf_compress:
    input:
        "{PATH}.vcf"
    output:
        vcf=temp("{PATH}.vcf.gz"),
        index=temp("{PATH}.vcf.gz.csi"),
    shell:
        "bgzip -c -f {input} > {output.vcf}; bcftools index -f {output.vcf};"

#######################################################
#Merge vcfs for a group
rule vcf_merge:
    input:
        ref="{PATH}/{REF}{EXT}fasta",	
        vcfs=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
        index=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz.csi", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2)
    output:
        "{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}freebayes.bcfmerged.vcf",
    threads:
        int(BCFTOOLS_MERGE_THREADS)
    shell:
        "bcftools merge --threads {threads} -m all {input.vcfs} | bcftools norm -f {input.ref} > {output};"

#######################################################
#group-wise per-locus alignments + supermatrix
rule group_supermatrix:
	input:
		fastas=lambda wildcards: expand("{PATH}/{ID}{EXT}.fasta", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], EXT=wildcards.EXT)
	output:
		list=temp("{PATH}/{GRP}{EXT}.fasta.list"),
		dir=directory("{PATH}/{GRP}{EXT}.matrix"),
		alignedlist=temp("{PATH}/{GRP}{EXT}.fasta.alignedlist"),
		matrix="{PATH}/{GRP}{EXT}.matrix.fasta",
		part="{PATH}/{GRP}{EXT}.matrix.fasta.part"
	threads:
		int(GROUPALIGN_THREADS)
	shell:
		""" echo {input.fastas} | tr " " "\n" > {output.list}; """
		""" python {TOOLDIR}/scripts/per_locus_alignments.py {output.list} {output.dir} "mafft {GROUPALIGN_MAFFT_OPT}" {threads}; """	
		""" ls {output.dir} | grep aligned.fasta$ | awk -v dir={output.dir} '{{print dir"/"$0}}' > {output.alignedlist}; """
		""" python {TOOLDIR}/scripts/fastas2Matrix.py {output.alignedlist} {output.matrix}; """	

#######################################################
rule perlocus_alignments_extract:
    input:
        dir=directory("{PATH}/{GRP}{EXT}.matrix"),
    output:
        "{PATH}/{GRP}{EXT}.{GENE}_singlegene.fasta",
    shell: 
        "cp {input.dir}/{wildcards.GENE}.aligned.fasta {output};"
        


#######################################################
#run trimal on fasta
rule trimal_run:
    input:
        "{PATH}.fasta"
    output:
        "{PATH}.trimal.fasta"
    shell:
        """
        if [ -s {input} ]; then
            trimal -in {input} {TRIMAL_OPT} | awk '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]""$0}}; END{{PROCINFO["sorted_in"]="@ind_str_asc"; for(f in seq){{print f; print seq[f]}}}}' > {output}
        else
            touch {output}
        fi;
        """

#######################################################
#run iqtree on fasta
rule iqtree_run:
    input:
        "{PATH}.fasta"
    output:
        tree="{PATH}.fasta.treefile",
        repor="{PATH}.fasta.iqtree",
        log=temp("{PATH}.fasta.log"),
        bionj=temp("{PATH}.fasta.bionj"),
        ckp=temp("{PATH}.fasta.ckp.gz"),
        cont=temp("{PATH}.fasta.contree"),
        mld=temp("{PATH}.fasta.mldist"),
        nex=temp("{PATH}.fasta.splits.nex")
    threads: int(IQTREE_THREADS)
    shell:
        "iqtree -s {input} {IQTREE_OPT} -nt {threads} -st DNA;"

#######################################################
#estimate number of heterozygous sites from a bam file per contig using ANGSD and bam file
rule angsd_run:
    input:
        bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
        bamindex="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam.bai",
        cont="{PATH}/{REF}{EXT}fasta",
        refindex="{PATH}/{REF}{EXT}fasta.fai",
    output:
        bamlist=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.bamlist"),
        safidx=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.saf.idx"),
        safpos=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.saf.pos.gz"),
        saf=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.saf.gz"),
        sfs=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.sfs"),
        thetaidx=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.thetas.idx"),
        theta=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.thetas.gz"),
        pestpg=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.thetas.idx.pestPG"),
        arg=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.arg"),
        stats="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.stats.txt",
        counts="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd.counts.txt",
    params:
        prefix="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}angsd",
        env="angsd_env"
    threads: 1
    shell:
##        "samtools faidx {input.cont};"
        "touch {input.refindex};"
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "echo {input.bam} > {output.bamlist};"
        "angsd -P {threads} -anc {input.cont} -doSaf 1 -GL 2 -bam {output.bamlist} -out {params.prefix};"
        "realSFS {output.safidx} > {output.sfs};"
        "angsd -P {threads} -anc {input.cont} -doThetas 1 -doSaf 1 -GL 2 -pest {output.sfs} -bam {output.bamlist} -out {params.prefix};"
        "thetaStat do_stat {output.thetaidx};"
        "cp {output.pestpg} {output.stats};"
        """ awk 'NR!=1{{print $2" "$14+1" "$4}}' {output.stats} | sort -k 1 > {output.counts}; """

#######################################################
#count number of heterozygous snps per contigs from vcf. Applies contamination filter and takes bed into account to compute contig length.
rule het_sites_count_run:
    input:
        covbed="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed",
        vcf="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf",
    output:
        counts="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}contamhet.counts.txt",
    shell:
        "python {TOOLDIR}/scripts/vcf_hetcontam_count.py {input.vcf} {input.covbed} {HETCOUNTS_CONTAM_OPT} > {output.counts};"
		
#######################################################
#compile divergence model for use by stan (should run only once)
rule stan_model_train:
    input:
        "{PATH}.stan"
    output:
        "{PATH}.rds"
    params:
        env="stan_env",
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "Rscript {TOOLDIR}/stan/stan_train_model.R {input};"

#######################################################
#estimate divergence using stan
rule divergence_estimations_run:
    input:
        rds=expand("{TOOLS}/stan/divergence_model.rds", TOOLS=TOOLDIR),
        counts="{PATH}.counts.txt",
    output:
        estimates="{PATH}.divestim.txt",
        posterior="{PATH}.divposterior.txt",
    threads:
        int(STAN_THREADS)
    params:
        env="stan_env",
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "Rscript {TOOLDIR}/scripts/stan_divergence_estimation.R {input.counts} {output.estimates} {output.posterior} {input.rds} {STAN_OPT} {threads};"

#######################################################
#gather results for all ids in a group
rule divergence_estimations_gather:
    input:
        lambda wildcards: expand("{PATH}/{ID}{EXT}txt", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], EXT=wildcards.EXT)
    output:
        "{PATH}/{GRP}{EXT}gathered.txt",
    shell:
        """ awk 'NR==1{{print}}; FNR!=1{{print}}' {input} > {output}; """



