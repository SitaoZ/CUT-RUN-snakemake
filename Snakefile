configfile: "config.yaml"

# software path 

WORKDIR = config["workdir"]
FILES = config["samples"]
SAMPLE = FILES.keys()
USED_GTF = config["GTF"]
USED_GFF = config["GFF"]


ADAPTER_FORWARD = config['adapter']['forward']
ADAPTER_BACKWARD = config['adapter']['backward']

CUTADAPT = config['software']['cutadapt']
BOWTIE2 = config['software']['bowtie2']
HISAT2 = config['software']['hisat2']
SAMTOOLS = config['software']['samtools']
RSCRIPT = config['software']['Rscript']
MACS2 = config['software']['macs2']


def all_input(wildcards):
    all_result = []
    # qc
    for sample in SAMPLE:
        # all_result.extend(expand(["data/samples/{sample}.fastq.gz"], sample=sample))
        # all_result.extend(expand(['01.fastqc/{sample}_fastqc.html', 
	#			 '01.fastqc/{sample}_fastqc.zip'], sample=sample))
        # all_result.extend(expand(["02.cutadapt/{sample}.clean.fq.gz"], sample=sample))
        # all_result.extend(expand(["03.sortmerna/{sample}/{sample}_non_rRNA.fq.gz"], sample=sample))
        # all_result.extend(expand(["04.bowtie2/{sample}/{sample}.bam"], sample=sample))
        # all_result.extend(expand(["05.sortbam/{sample}/{sample}.sort.bam"], sample=sample))
        # all_result.extend(expand(["06.markdup/{sample}/{sample}.rmdup.sort.bam"], sample=sample))
        all_result.extend(expand(["07.macs2/{sample}/{sample}.rmdup.sort.bam_peaks.xls"], sample=sample))
    return all_result

rule all:
    input: all_input

def get_input_fastqs(wildcards):
    sample_name = wildcards.sample
    return FILES[sample_name]

rule prepare_fq:
    input: 
         fq = lambda wildcards: FILES[wildcards.sample]
    output:
        "data/samples/{sample}.fastq.gz"
    log:
        stdout="logs/prepare_fq.{sample}.stdout",
        stderr="logs/prepare_fq.{sample}.stderr"
    shell:
        "if [ ! -d data/samples ]; then mkdir -p data/samples;fi; ln -s {input.fq} data/samples"

rule fastqc:
    input:
        "data/samples/{sample}.fastq.gz"
    output:
        out1="01.fastqc/{sample}_fastqc.html",
        out2="01.fastqc/{sample}_fastqc.zip"
    threads: 16
    resources:
        tmpdir="tmp"
    message: "fastqc {input}: {threads} threads"
    benchmark:
        "benchmarks/fastqc.{sample}.csv"
    log:
        stdout="logs/fastqc.{sample}.stdout",
        stderr="logs/fastqc.{sample}.stderr"
    shell:
        "fastqc -t {threads} -o 01.fastqc {input} > {log.stdout} 2>{log.stderr}"

rule cutadapt:
    input:
        qc_html="01.fastqc/{sample}_fastqc.html",
        qc_zip="01.fastqc/{sample}_fastqc.zip",
        read="data/samples/{sample}.fastq.gz"
    output:
        "02.cutadapt/{sample}.clean.fq.gz",
    threads: 16
    resources:
        tmpdir="tmp"
    message: "cutadapt {input.read}: {threads} threads"
    benchmark:
        "benchmarks/cutadapt.{sample}.csv"
    log:
        stdout="logs/cutadapt.{sample}.stdout",
        stderr="logs/cutadapt.{sample}.stderr"
    params:
        "-m 20 -q 20 --max-n 0.2 -e 0.08 -Z --trim-n"
    shell:
        "{CUTADAPT} {params} -j {threads} -a {ADAPTER_FORWARD} -o {output} {input.read}> {log.stdout} 2>{log.stderr}"

rule sortmerna:
    input:
        read="02.cutadapt/{sample}.clean.fq.gz"
    output:
        out="03.sortmerna/{sample}/{sample}_non_rRNA.fq.gz",
    threads: 16
    resources:
        tmpdir="tmp"
    message: "sortmerna {input.read} : {threads} threads"
    benchmark:
        "benchmarks/sortmerna.{sample}.csv"
    log:
        stdout="logs/sortmerna.{sample}.stdout",
        stderr="logs/sortmerna.{sample}.stderr"
    params:
        fixed = "--sam --num_alignments 1 --fastx",
        aligned = "03.sortmerna/{sample}/{sample}_rRNA",
        other = "03.sortmerna/{sample}/{sample}_non_rRNA",
        workdir = "03.sortmerna/{sample}",
        idx = "03.sortmerna/{sample}/idx",
        kvdb = "03.sortmerna/{sample}/kvdb",
        readb = "03.sortmerna/{sample}/readb"
    shell:
        "sortmerna {params.fixed} --threads {threads} --workdir {params.workdir} \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \
         --ref /home/zhusitao/software/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \
         --reads {input.read} --aligned {params.aligned} --other {params.other} -v > {log.stdout} 2>{log.stderr};"
         "rm -rf {params.idx} {params.kvdb} {params.readb}"

rule bowtie2:
    input:
        "03.sortmerna/{sample}/{sample}_non_rRNA.fq.gz"
    output:
        "04.bowtie2/{sample}/{sample}.bam"
    threads: 16
    resources:
        tmpdir = "tmp"
    message: "bowtie2 {input} : {threads} threads"
    benchmark:
        "benchmarks/bowtie2.{sample}.csv"
    log:
        stdout="logs/bowtie2.{sample}.stdout",
        stderr="logs/bowtie2.{sample}.stderr"
    params:
        "-q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
         --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200"
    shell:
        "{BOWTIE2} {params} -p {threads} -x /home/zhusitao/database/plant/ath/tair10/RSEM/tair10_rsem \
         -U {input} 2>{log.stderr} | {SAMTOOLS} view -b -o {output} -"

rule sort_bam:
    input:
        "04.bowtie2/{sample}/{sample}.bam"
    output:
        "05.sortbam/{sample}/{sample}.sort.bam"
    threads: 16
    resources:
        tmpdir="tmp"
    message: "samtools sort {input} : {threads} threads"
    benchmark:
        "benchmarks/sort_bam.{sample}.csv"
    log:
        stdout="logs/sort_bam.{sample}.stdout",
        stderr="logs/sort_bam.{sample}.stderr"
    threads: 16
    shell:
        """
        export TMPDIR=05.sortbam/{wildcards.sample}/tmp
        mkdir -p ${{TMPDIR}}
        {SAMTOOLS} sort -@ {threads} -T ${{TMPDIR}} -o {output} {input}
        """

rule mark_dup:
    input:
        "05.sortbam/{sample}/{sample}.sort.bam"
    output:
        "06.markdup/{sample}/{sample}.rmdup.sort.bam"
    threads: 1
    resources:
        tmpdir="tmp"
    message: "picards {input} : {threads} threads"
    benchmark:
        "benchmarks/mark_dup.{sample}.csv"
    log:
        stdout="logs/markdup.{sample}.stdout",
        stderr="logs/markdup.{sample}.stderr"
    shell:
        "java -jar /home/zhusitao/software/picard/picard.jar MarkDuplicates  --REMOVE_DUPLICATES true  --INPUT {input}  --OUTPUT {output} --METRICS_FILE 06.markdup/{wildcards.sample}/{wildcards.sample}.sort.matrix > {log.stdout} 2> {log.stderr}"


rule macs2:
    input:
        "06.markdup/{sample}/{sample}.rmdup.sort.bam"
    output:
        "07.macs2/{sample}/{sample}.rmdup.sort.bam_peaks.xls"
    threads: 1
    resources:
        tmpdir="tmp"
    message: "picards {input} : {threads} threads"
    benchmark:
        "benchmarks/macs2.{sample}.csv"
    log:
        stdout="logs/macs2.{sample}.stdout",
        stderr="logs/macs2.{sample}.stderr"
    params:
        base_file="{sample}.rmdup.sort.bam",
        outdir="07.macs2/{sample}",
        genome_size=100000000
    shell:
        """
        {MACS2} callpeak -t {input} -g {params.genome_size} -f BAM -n {params.base_file} --outdir {WORKDIR}/{params.outdir} -q 0.01 -B --SPMR --keep-dup all --nomodel 2>{log.stderr}
        """

