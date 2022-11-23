
import os 
import numpy as np
import re


##PARAMETERS (feel free to edit)
#Path to directory for data files
DATADIR="/media/bigvol/arthur/phylogenies/phylos/data"
#Path to tools directory
TOOLDIR="/media/bigvol/arthur/phylogenies/phylos/tools"
#Path to run file (group, id)
RUN_FILE="/media/bigvol/arthur/phylogenies/phylos/phylo_list.my"
MOM_FILE="/media/bigvol/arthur/phylogenies/phylos/mother_daughter_list.my2"

FASTP_THREADS=5
FASTP_OPT="-q 20 -u 70 -n 40 -l 40 -w 1"
MEGAHIT_THREADS=5
MEGAHIT_OPT="--k-min 31 --k-max 101 --k-step 10"
BUSCO_THREADS=5
BUSCO_OPT="-m genome"
BWA_THREADS=3
BWA_OPT="-k 19"
STAN_THREADS=2
STAN_OPT="2 1000 10000"
FREEBAYES_OPT="--min-alternate-count 1 -z 0.05"
FREEBAYES_GROUP_THREADS=30
VCF_COMPRESS_THREADS=3
BCFTOOLS_CONSENSUS_MINDP=4 #Min depth to report nucleotide in consensus. Positions with less will be "-"
BCFTOOLS_CONSENSUS_MAXGAPPROP=0.8 #remove consensus sequences with more than this proportion of gaps.
BCFTOOLS_MERGE_MINDP=4 #Min depth for every individual to consider region in vcf merging
BCFTOOLS_MERGE_MINDP=2 #Min depth for every individual to consider region in vcf merging
BCFTOOLS_MERGE_THREADS=40
PLINK_GDIST_OPT="--distance square0 1-ibs"
PLINK_GDIST_THREADS=40
VCFTOOLS_FILTER_OPT="--remove-indels --min-meanDP 2.5 --min-alleles 2 --maf 0.02  --max-missing 0.9"
HETCOUNTS_CONTAM_OPT="10 0.1" #1: min dp to consider snp; #2: alt allele prob for contam filter
GROUPALIGN_MAFFT_OPT="--auto --adjustdirection"
GROUPALIGN_THREADS=40
TRIMAL_OPT="-gt 0.93 -st 0.8 -w 5"
IQTREE_OPT="-bb 1000"
IQTREE_OPT="-m HKY+I+F+G4 -bb 1000"
IQTREE_THREADS=40
GENOMESCOPE_THREADS=4
FIS_MATRIX_THREADS=40
LDHELMET_MAFFT_OPT="--auto --adjustdirection"
LDHELMET_THETA=0.002
LDHELMET_THREADS=5
PHYLOFLASH_THREADS=5


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
    FILE="[^/]*",
    KMER="[0-9]*"


ruleorder:
	fastp_pese_run > fastp_pe_run > fastp_se_run
ruleorder:
	megahit_pese_run > megahit_pe_run > megahit_se_run
ruleorder:
	bwa_pese_run > bwa_pe_run > bwa_se_run
ruleorder:
	bcftools_consensus_run > mother_phase_run
#ruleorder:
#    get_backup_consensus_sequences > bcftools_consensus_run
ruleorder:
	genomescope_pese_run > genomescope_pe_run > genomescope_se_run

#######################################################
#####TARGETS###########################################
#######################################################
#Require final output. Comment out unwanted output.
rule all: 
    input: 
		#[ID]-[REFID][EXTENSION FOR REF FILE]bwa[EXTENSION FOR FASTQ OF ID]freebayes[EXT]
        ###SCANS
        #expand("{PATH}/{ID}-Y15260_1_buscoref.bwa.fastp.{MET}.divestim.txt", PATH=DATADIR, ID=GRPDIC["switchNovogene"], MET="angsd contamhet".split(" ")),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.{MET}.divestim.gathered.txt", PATH=DATADIR, GRP="Messor".split(" "), MET="angsd contamhet".split(" ")),
        ###CONSENSUS SEQUENCES ONLY
        #expand("{PATH}/{ID}-Mscabrinodis_mitoref.bwa.fastp.freebayes.bcfconsensus.fasta", PATH=DATADIR, ID=GRPDIC["switchNovogene"]),
        #expand("{PATH}/{ID}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.fasta", PATH=DATADIR, ID=GRPDIC["switchNovogene"]),
        ###"phased" fasta for DILS
        #expand("{PATH}/{ID}-allele1-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.dils.fasta", PATH=DATADIR, ID=GRPDIC["Messor"]),
		###PHYLOS
        ##Mito
        #Full
        expand("{PATH}/{GRP}-Mscabrinodis_mitoref.bwa.fastp.freebayes.bcfconsensus.matrix.fasta.treefile", PATH=DATADIR, GRP="structoribericusMales".split(" ")),
        #Per-gene
		#expand("{PATH}/{GRP}-Mscabrinodis_mitoref.bwa.fastp.freebayes.bcfconsensus.{GENE}_singlegene.trimal.fasta.treefile", PATH=DATADIR, GRP="Messor".split(" "), GENE=MITOGENES),
        ##Nuclear
        #Full
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.matrix.fasta", PATH=DATADIR, GRP="MessorphasedNobad".split(" ")),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.matrix.trimal.fasta.treefile", PATH=DATADIR, GRP="MessorphasedNobad".split(" ")),
        #Per-gene
		#expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.{GENE}_singlegene.trimal.fasta.treefile", PATH=DATADIR, GRP="Messor".split(" "), GENE=NUCGENES),
        ###SNPS
        #expand("{PATH}/{ID}-Y15260_1_buscoref.bwa.fastp.rmonofreebayes.vcf.gz", PATH=DATADIR, ID=GRPDIC["switchNovogene"]),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayesgrp.filtered.vcf", PATH=DATADIR, GRP="barbarusref structorref".split(" ")),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayesgrp.vcf", PATH=DATADIR, GRP="barbarusref structorref".split(" ")),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.rmonofreebayes.bcfmerged.filtered.vcf", PATH=DATADIR, GRP="allstructor".split(" ")),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.rmonofreebayes.bcfmerged.keep-{GRP2}.filtered.vcf", PATH=DATADIR, GRP="allstructor".split(" "), GRP2="ibericusMales structorMMales structor3Males".split(" ")),
        ###GENETIC DISTANCES
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.rmonofreebayes.bcfmerged.nocontam0.1.filtered.plinkgdist.mdist.txt", PATH=DATADIR, GRP="allstructor".split(" ")),
        #expand("{PATH}/steiner.plinkgdist.mdist.txt", PATH=DATADIR, GRP="allstructor".split(" ")),
        ###NQUIRE
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.nquire.gathered.txt", PATH=DATADIR, GRP="Messor".split(" "))
        ###SMUDGEPLOT
        #expand("{PATH}/{ID}.fastp.gscope_k{KMER}.histo", PATH=DATADIR, ID="SH02-23 SH02-24".split(" "), KMER="11 15 21 26 31 36".split(" "))
        ###FIS MATRIX
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.rmonofreebayes.bcfmerged.keep-allstructorQueens.fis.total.txt", PATH=DATADIR, GRP="allstructor".split(" ")),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.rmonofreebayes.bcfmerged.fis.total.txt", PATH=DATADIR, GRP="allstructor".split(" ")),
        ###LDHELMET
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.freebayes.bcfconsensus.contig.{CONT}.mafft.ldhelmet.conf", PATH=DATADIR, GRP="WorkersFantomPat WorkersFantomMat WorkersStructorPat WorkersStructorMat".split(" "), CONT="11649at7399 123at7399 1602at7399 1895at7399 3237at7399 415at7399 439at7399 60at7399 61at7399 78at7399 88at7399".split(" ")),
        ###COVERAGE
        #expand("{PATH}/{ID}-Y15260_1_buscoref.bwa.fastp.coverage.overall.txt", PATH=DATADIR, ID=GRPDIC["switchNovogene"]),
        #expand("{PATH}/{GRP}-Y15260_1_buscoref.bwa.fastp.coverage.overall.gathered.txt", PATH=DATADIR, GRP="Messor".split(" ")),
        ###PHYLOFLASH
        #expand("{PATH}/RINSAM_phyloflash.phyloFlash.tar.gz", PATH=DATADIR, ID=GRPDIC["Messor"]),
        ###INDICES
        #expand("{PATH}/{ID}.indexcnt.txt", PATH=DATADIR, ID=GRPDIC["Messor2022"]),
        #expand("{PATH}/{ID}.lanecnt.txt", PATH=DATADIR, ID=GRPDIC["Messor2022"]),

#######################################################
#####RULES#############################################
#######################################################

#######################################################
#check index
rule index_pe_run:
    input:	
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz"
    output:
        out="{PATH}/{ID}{EXT}indexcnt.txt",
    shell:
        """ zcat {input.forward} | awk -v id={wildcards.ID} -v FS=":" 'NR%4==1{{!x[$10]++}}; END{{for(i in x){{print id,"forward",i,x[i]}}}}' > {output.out}; """
        """ zcat {input.reverse} | awk -v id={wildcards.ID} -v FS=":" 'NR%4==1{{!x[$10]++}}; END{{for(i in x){{print id,"reverse",i,x[i]}}}}' >> {output.out}; """
 
#######################################################
#check index
rule lane_pe_run:
    input:	
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz"
    output:
        out="{PATH}/{ID}{EXT}lanecnt.txt",
    shell:
        """ zcat {input.forward} | awk -v id={wildcards.ID} -v FS=":" 'NR%4==1{{!x[$2":"$3":"$4]++}}; END{{for(i in x){{print id,"forward",i,x[i]}}}}' > {output.out}; """
        """ zcat {input.reverse} | awk -v id={wildcards.ID} -v FS=":" 'NR%4==1{{!x[$2":"$3":"$4]++}}; END{{for(i in x){{print id,"reverse",i,x[i]}}}}' >> {output.out}; """

 
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
#call snp using freebayes and report monomorphic sites as well
rule freebayes_rmono_run:
	input:
		bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
		ref="{PATH}/{REF}{EXT}fasta",
		refindex="{PATH}/{REF}{EXT}fasta.fai",
	output:
		vcf="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}rmonofreebayes.vcf",
	shell:
		"freebayes -f {input.ref} --report-monomorphic {FREEBAYES_OPT} {input.bam} | bcftools norm -f {input.ref} - > {output}"


#######################################################
#call snp using freebayes
rule freebayes_group_run:
    input:
        bams=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
        bais=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam.bai", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
        ref="{PATH}/{REF}{EXT}fasta",
        refindex="{PATH}/{REF}{EXT}fasta.fai",
    output:
        list="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}freebayesgrp.list",
        vcf="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}freebayesgrp.vcf", 
    threads: int(FREEBAYES_GROUP_THREADS)
    shell:
        "echo {input.bams} | tr ' ' '\n' > {output.list};"
        "cp {output.list} {output.list}plop;"
        "freebayes-parallel <(fasta_generate_regions.py {input.refindex} 100000) {threads} -f {input.ref} --bam-list {output.list} {FREEBAYES_OPT} | bcftools norm -f {input.ref} - > {output.vcf};"
##        "freebayes -f {input.ref} --bam-list {output.list} {FREEBAYES_OPT} | bcftools norm -f {input.ref} - > {output.vcf};"


#######################################################
#Filter snps with samtools:
rule vcf_filter:
    input:
        vcf="{PATH}.vcf.gz"
    output:
        vcf="{PATH}.filtered.vcf"
    shell:
        "vcftools --gzvcf {input.vcf} {VCFTOOLS_FILTER_OPT} --recode -c > {output.vcf};"

#######################################################
#Filter snps with samtools:
rule vcf_keep:
    input:
        vcf="{PATH}.vcf.gz"
    output:
        vcf="{PATH}.keep-{GRP}.vcf",
        ids=temp("{PATH}.keep-{GRP}.ids")
    params:
        ids=lambda wildcards: expand("{ID}", ID=GRPDIC[wildcards.GRP])
    shell:
        """ echo {params.ids} | tr " " "\n" > {output.ids}; """
        "vcftools --gzvcf {input.vcf} --keep {output.ids} --recode -c > {output.vcf};"


#######################################################
#compute genetic distance matrix with plink
rule plink_gdist_matrix:
    input:
        vcf="{PATH}.vcf.gz"
    output:
        log="{PATH}.plinkgdist.log",
        mdist=temp("{PATH}.plinkgdist.mdist"),
        ids="{PATH}.plinkgdist.mdist.id",
        sex=temp("{PATH}.plinkgdist.nosex"),
        mdist2="{PATH}.plinkgdist.mdist.txt",
    params:
        pre="{PATH}.plinkgdist"
    threads:
        int(PLINK_GDIST_THREADS)
    shell:
        "plink --allow-extra-chr {PLINK_GDIST_OPT} --threads {threads} --vcf {input.vcf} --out {params.pre};"
        """ awk -v FS="\t" -v OFS="\t" -v ids="`cut -f1 {output.ids}`" 'BEGIN{{split(ids, ids2, "\\n"); head="ids"; for(i in ids2){{head=head""OFS""ids2[i]}}; print head}}; {{$0=ids2[NR]""OFS""$0; print $0}}' {output.mdist} > {output.mdist2}; """


#######################################################
#run bedtools to compute coverage
rule bedtools_genomecov_run:
	input:
		bam="{PATH}.bam",
	output:
		bed="{PATH}.coverage.bed",
	shell:
		"bedtools genomecov -bga -ibam {input.bam} > {output.bed};"
	
#######################################################
#compute mean coverage
rule mean_genomecov_run:
	input:
		bed="{PATH}/{ID}-{EXT}coverage.bed",
	output:
		pcont="{PATH}/{ID}-{EXT}coverage.percontig.txt",
		all="{PATH}/{ID}-{EXT}coverage.overall.txt",
	shell:
	    "python {TOOLDIR}/scripts/mean_coverage.py {input.bed} {wildcards.ID} {output.pcont} {output.all} {BCFTOOLS_CONSENSUS_MINDP};"

#######################################################
#gather mean coverage
rule mean_genomecov_gather:
    input:
        all=lambda wildcards: expand("{PATH}/{ID}-{EXT}coverage.overall.txt", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], EXT=wildcards.EXT),
        pcont=lambda wildcards: expand("{PATH}/{ID}-{EXT}coverage.percontig.txt", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], EXT=wildcards.EXT),
    output:
        all="{PATH}/{GRP}-{EXT}coverage.overall.gathered.txt",
        pcont="{PATH}/{GRP}-{EXT}coverage.percontig.gathered.txt",
    shell:
        """ awk 'FNR==2 || NR == 1' {input.all} > {output.all}; """
        """ awk 'FNR>1 || NR == 1' {input.pcont} > {output.pcont}; """
	
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
#rule get_backup_consensus_sequences:
#    input:
#        "{PATH}/consensus_sequences/{FILE}.bcfconsensus.backup.fasta"
#    output:
#        "{PATH}/{FILE}.bcfconsensus.fasta"
#    shell:
#        "cp {input} {output};"
        
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
#Output alelles in two different fasta (produce DILS input)
rule dils_phase_run:
	input:
		vcf="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf",
		covbed="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed",
		ref="{PATH}/{REF}{EXT}fasta",			
	output:
		badcovbed=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}badcoverage.bed"),
		vcfpha=temp("{PATH}/{ID}-totphased-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz"),
		vcfphaindex=temp("{PATH}/{ID}-totphased-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz.csi"),
		all1="{PATH}/{ID}-allele1-{REF}{EXT}bwa{EXT2}freebayes.bcfconsensus.dils.fasta",
		all2="{PATH}/{ID}-allele2-{REF}{EXT}bwa{EXT2}freebayes.bcfconsensus.dils.fasta",
	shell:
		"awk -v min={BCFTOOLS_CONSENSUS_MINDP} '$4<min' {input.covbed}  > {output.badcovbed};"
		"python {TOOLDIR}/scripts/vcf_dils_phase.py {input.vcf} | bgzip > {output.vcfpha};"
		"bcftools index -f {output.vcfpha};"
		""" bcftools consensus -s {wildcards.ID}_all1 --mask {output.badcovbed} --mask-with "-" -f {input.ref} {output.vcfpha} | awk -v maxgap={BCFTOOLS_CONSENSUS_MAXGAPPROP} '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]""$0}}; END{{for(i in seq){{l=seq[i]; if((gsub("-", "", l)/length(seq[i]))<maxgap){{print i; print seq[i]}}}}}}' > {output.all1}; """
		""" bcftools consensus -s {wildcards.ID}_all2 --mask {output.badcovbed} --mask-with "-" -f {input.ref} {output.vcfpha} | awk -v maxgap={BCFTOOLS_CONSENSUS_MAXGAPPROP} '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]""$0}}; END{{for(i in seq){{l=seq[i]; if((gsub("-", "", l)/length(seq[i]))<maxgap){{print i; print seq[i]}}}}}}' > {output.all2}; """


#######################################################
#Compress vcf
rule vcf_compress:
    input:
        "{PATH}.vcf"
    output:
        vcf="{PATH}.vcf.gz",
        index="{PATH}.vcf.gz.csi",
    threads:
        int(VCF_COMPRESS_THREADS)
    shell:
        "bgzip -c -@ {threads} -f {input} > {output.vcf}; bcftools index -f {output.vcf};"

#######################################################
#Merge vcfs for a group
rule vcf_merge:
    input:
        ref="{PATH}/{REF}{EXT}fasta",	
        vcfs=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
        index=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}freebayes.vcf.gz.csi", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
        beds=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}coverage.bed", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
    output:
        badbed="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}badcoverage.bed",
        goodbed="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}goodcoverage.bed",
        gen="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}refgenome.bed",
        vcf="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}freebayes.bcfmerged.vcf",
    threads:
        int(BCFTOOLS_MERGE_THREADS)
    shell:
        """ awk -v OFS="\t" '/>/{{head=$0; sub("^>", "", head)}}; !/>/{{print head,length($0)}}' {input.ref} | sort -k1,1 > {output.gen}; """
        "echo 'Computing coverage!';"
        "awk -v min={BCFTOOLS_MERGE_MINDP} '$4<min' {input.beds} | sort -k1,1 -k2,2n | bedtools merge > {output.badbed};"
        "bedtools complement -L -i {output.badbed} -g {output.gen} > {output.goodbed};"
        "echo 'Merging vcf files!';"
        "bcftools merge --threads {threads} -R {output.goodbed} -0 -m all {input.vcfs} | bcftools norm -f {input.ref} > {output.vcf};"

#######################################################
#Merge vcfs for a group
rule vcf_rmono_merge:
    input:
        ref="{PATH}/{REF}{EXT}fasta",	
        vcfs=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}rmonofreebayes.vcf.gz", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
        index=lambda wildcards: expand("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}rmonofreebayes.vcf.gz.csi", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], REF=wildcards.REF, EXT=wildcards.EXT, EXT2=wildcards.EXT2),
    output:
        vcf="{PATH}/{GRP}-{REF}{EXT}bwa{EXT2}rmonofreebayes.bcfmerged.vcf",
    threads:
        int(BCFTOOLS_MERGE_THREADS)
    shell:
        "bcftools merge --threads {threads} -m all {input.vcfs} | bcftools norm -f {input.ref} > {output.vcf};"

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
#Ensure all individuals in alignment are required from RUN_LIST. Useful to remove individuals post alignment.
rule alignment_trim:
    input:
        "{PATH}/{GRP}{EXT}.fasta"
    output:
        out="{PATH}/{GRP}{EXT}.trimmed.fasta",
        inds=temp("{PATH}/{GRP}{EXT}.trimmed.inds.txt")
    params:
        inds=lambda wildcards: expand("{ID}", ID=GRPDIC[wildcards.GRP])
    shell:
        "echo {params.inds} | tr ' ' '\n' | sed 's/^/>/g' > {output.inds};"
        """ grep -A1 -f {output.inds} {input} | grep -v "^\-\-$" > {output.out}; """


#######################################################
#run iqtree on fasta
rule iqtree_run:
    input:
        "{PATH}.fasta"
    output:
        tree="{PATH}.fasta.treefile",
        repor="{PATH}.fasta.iqtree",
        log="{PATH}.fasta.log",
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
rule estimations_gather:
    input:
        lambda wildcards: expand("{PATH}/{ID}{EXT}txt", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], EXT=wildcards.EXT)
    output:
        "{PATH}/{GRP}{EXT}gathered.txt",
    shell:
        """ awk 'NR==1{{print}}; FNR!=1{{print}}' {input} > {output}; """

#######################################################
#Estimate ploidy with nQuire
rule nquire_run:
    input:
        bam="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam",
        bamindex="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}bam.bai",
    output:
        bin=temp("{PATH}/{ID}-{REF}{EXT}bwa{EXT2}nquire.bin"),
        est="{PATH}/{ID}-{REF}{EXT}bwa{EXT2}nquire.txt"
    shell:
        "{TOOLDIR}/nQuire/nQuire create -f 0.2 -b {input.bam} -o {wildcards.PATH}/{wildcards.ID}-{wildcards.REF}{wildcards.EXT}bwa{wildcards.EXT2}nquire;"
        "cp {output.bin} {output.bin}plop;"
        """ {TOOLDIR}/nQuire/nQuire estmodel {output.bin} | awk -v head="file" -v val={output.est} -v FS=" = " '/loglik/{{sub("^.* ", "", $0); head=head" nQ_logLik"; val=val" "$0}}; /a =/{{cnt2++}}; /=/{{gsub(" ", "", $1); head=head" nQ_"cnt2"_"$1; val=val" "$2}}; END{{print head; print val}}' > {output.est}; """


#######################################################
#Align with MACSE
rule macse_run:
    input:
        "{PATH}.fasta"
    output:
        nt="{PATH}_NT.fasta",
        aa="{PATH}_AA.fasta"
    shell:
        "java -jar {TOOLDIR}/macse/macse_v2.06.jar -prog alignSequences -seq {input};"
    
#######################################################
#Run genomescope
rule genomescope_pese_run:
    input:	
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz",
        single="{PATH}/{ID}{EXT}fastq.gz"
    output:
        jf=temp("{PATH}/{ID}{EXT}gscope_k{KMER}.jf"),
        histo="{PATH}/{ID}{EXT}gscope_k{KMER}.histo",
        out=directory("{PATH}/{ID}{EXT}gscope_k{KMER}")
    threads: int(GENOMESCOPE_THREADS)
    params:
        env="genomescope_env",
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "zcat {input} | jellyfish count /dev/fd/0 -C -m {wildcards.KMER} -s 1000000000 -t {threads} -o {output.jf};"
        "jellyfish histo -t {threads} {output.jf} > {output.histo};"
        "genomescope2 -i {output.histo} -o {output.out} -k {wildcards.KMER};"


#######################################################
#Run genomescope
rule genomescope_pe_run:
    input:	
        forward="{PATH}/{ID}_1{EXT}fastq.gz",
        reverse="{PATH}/{ID}_2{EXT}fastq.gz",
    output:
        jf=temp("{PATH}/{ID}{EXT}gscope_k{KMER}.jf"),
        histo="{PATH}/{ID}{EXT}gscope_k{KMER}.histo",
        out=directory("{PATH}/{ID}{EXT}gscope_k{KMER}")
    threads: int(GENOMESCOPE_THREADS)
    params:
        env="genomescope_env",
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "zcat {input} | jellyfish count /dev/fd/0 -C -m {wildcards.KMER} -s 1000000000 -t {threads} -o {output.jf};"
        "jellyfish histo -t {threads} {output.jf} > {output.histo};"
        "genomescope2 -i {output.histo} -o {output.out} -k {wildcards.KMER};"

#######################################################
#Run genomescope
rule genomescope_se_run:
    input:	
        single="{PATH}/{ID}{EXT}fastq.gz"
    output:
        jf=temp("{PATH}/{ID}{EXT}gscope_k{KMER}.jf"),
        histo="{PATH}/{ID}{EXT}gscope_k{KMER}.histo",
        out=directory("{PATH}/{ID}{EXT}gscope_k{KMER}")
    threads: int(GENOMESCOPE_THREADS)
    params:
        env="genomescope_env",
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "zcat {input} | jellyfish count /dev/fd/0 -C -m {wildcards.KMER} -s 1000000000 -t {threads} -o {output.jf};"
        "jellyfish histo -t {threads} {output.jf} > {output.histo};"
        "genomescope2 -i {output.histo} -o {output.out} -k {wildcards.KMER};"



#######################################################
#Run genomescope 
rule smudgeplot_run:
    input:
        jf="{PATH}/{ID}{EXT}gscope_k{KMER}.jf",
        histo="{PATH}/{ID}{EXT}gscope_k{KMER}.histo",
    output:
        cut="{PATH}/{ID}{EXT}gscope_k{KMER}_cutoffs.tsv",
        cov="{PATH}/{ID}{EXT}gscope_k{KMER}_coverages.tsv",
        seq="{PATH}/{ID}{EXT}gscope_k{KMER}_sequences.tsv",
        sumt="{PATH}/{ID}{EXT}gscope_k{KMER}_summary_table.tsv",
        sumv="{PATH}/{ID}{EXT}gscope_k{KMER}_verbose_summary.txt",
        warn="{PATH}/{ID}{EXT}gscope_k{KMER}_warnings.txt",
        plot="{PATH}/{ID}{EXT}gscope_k{KMER}_smudgeplot.png",
        plotl10="{PATH}/{ID}{EXT}gscope_k{KMER}_smudgeplot_log10.png"
    params:
        env="genomescope_env",
        id="{PATH}/{ID}{EXT}gscope_k{KMER}"
    threads: int(GENOMESCOPE_THREADS) #necessary to limit memory usage...
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "L=$(smudgeplot.py cutoff {input.histo} L); U=$(smudgeplot.py cutoff {input.histo} U); echo $L $U > {output.cut};"
        
        "jellyfish dump -c -L $L -U $U {input.jf} | smudgeplot.py hetkmers -o {params.id};"
        "smudgeplot.py plot {output.cov} -o {params.id};"
    
      

#######################################################
#This needs a db constructed in miniconda3/path/envs/phyloflash_env/lib/phyloFlash
#Run phyloflash
rule phyloflash_se_run:
    input:
        forward="{PATH}/{ID}_1.fastq.gz",
        reverse="{PATH}/{ID}_2.fastq.gz"
    output:
        tar="{PATH}/{ID}_phyloflash.phyloFlash.tar.gz",
        log="{PATH}/{ID}_phyloflash.phyloFlash.log",
        html="{PATH}/{ID}_phyloflash.phyloFlash.html",
    params:
        env="phyloflash_env",
        tmpdir="{PATH}/{ID}_phyloflash.phyloFlash"
    threads:
        int(PHYLOFLASH_THREADS)
    shell:        
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};"
        "phyloFlash.pl -lib {wildcards.ID}_phyloflash -read1 {input.forward} -read2 {input.reverse} -CPUs {threads} -skip_emirge -almosteverything;"
        "mv {wildcards.ID}_phyloflash.phyloFlash* {wildcards.PATH}/;"
        "rm -rf {params.tmpdir};"
        "mkdir -p {params.tmpdir};"
        "cp {output.tar} {params.tmpdir};"
        "cd {params.tmpdir};"
        "tar -xzf {wildcards.ID}_phyloflash.phyloFlash.tar.gz;"
        "rm -rf {params.tmpdir};"

 
#######################################################
#Run fis matrix
rule fis_matrix:
    input:
        vcf="{PATH}.vcf"
    output:
        pge="{PATH}.fis.pergene.txt",
        tot="{PATH}.fis.total.txt",
    threads: int(FIS_MATRIX_THREADS)
    shell:
        "python {TOOLDIR}/scripts/vcf_compute_all_pairwise_fis.py {TOOLDIR}/scripts/vcf_compute_pairwise_fis.py {input.vcf} {threads};"

#######################################################
#Run ldhelmet on a given contig for a given groups
rule ldhelmet_prepare:
    input:
        fastas=lambda wildcards: expand("{PATH}/{ID}{EXT}.fasta", PATH=wildcards.PATH, ID=GRPDIC[wildcards.GRP], EXT=wildcards.EXT)
    output:
        fasta=temp("{PATH}/{GRP}{EXT}.contig.{CONT}.fasta"),
        fasta2="{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.fasta"
    shell:
        """ grep -A1 ">{wildcards.CONT}" {input.fastas} | grep -v "^\-\-$" | awk '/>/{{sub(":>.*$", "", $0); sub("^.*\\\/", "", $0); $0=">"$0}}; !/>/{{sub("^.*fasta-", "", $0)}}; {{print}}' > {output.fasta}; """
        """ mafft {LDHELMET_MAFFT_OPT} {output.fasta} | awk '/>/{{head=$0}}; !/>/{{seq[head]=seq[head]''$0}}; END{{PROCINFO[\"sorted_in\"]=\"@ind_str_asc\"; for (i in seq){{print i; print toupper(seq[i])}}}}' | awk '/>/{{print}}; !/>/{{gsub("-", "N", $0); print}}' > {output.fasta2}; """
        
#######################################################
#Run ldhelmet on a given contig for a given groups
rule ldhelmet_run:
    input:
        fasta="{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.fasta"
    output:
        conf=temp("{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.ldhelmet.conf"),
        lk=temp("{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.ldhelmet.lk"),
        pade=temp("{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.ldhelmet.pade"),
        post=temp("{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.ldhelmet.post"),
        resu="{PATH}/{GRP}{EXT}.contig.{CONT}.mafft.ldhelmet.results.txt"
    params:
        env="ldhelmet_env",
    threads: int(LDHELMET_THREADS)
    shell:
        "set +eu;"
        ". $(conda info --base)/etc/profile.d/conda.sh;"
        "conda activate {params.env};" 
        "echo {threads};"
        "ldhelmet find_confs --num_threads {threads} -w 50 -o {output.conf} {input.fasta};" 
        "ldhelmet table_gen --num_threads {threads} -c {output.conf} -t {LDHELMET_THETA} -r 0.0 0.1 10.0 1.0 100.0 -o {output.lk};"
        "ldhelmet pade --num_threads {threads} -c {output.conf} -t {LDHELMET_THETA} -x 11 -o {output.pade};"
        "ldhelmet rjmcmc --num_threads {threads} -w 50 -l {output.lk} -p {output.pade} -b 50.0 -s {input.fasta} --burn_in 100000 -n 1000000 -o {output.post};"
        "ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o {output.resu} {output.post};"





