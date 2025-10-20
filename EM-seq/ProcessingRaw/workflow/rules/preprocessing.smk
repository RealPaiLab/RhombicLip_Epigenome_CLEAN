
import string

# Get config file
configfile: f"{workflow.basedir}/../config/config.yaml"

SAMPLES = []

path: str = config["directory"]["input"]

for fastq in os.scandir(path):
    if fastq.name[0] != '.':
        SAMPLES.append(fastq.name)

# Define input and output file paths
out_path: str = config["directory"]["output"]
fastqc_dir = f"{out_path}/prealignment/qc"
fastp_dir = f"{out_path}/prealignment/trim"

# Define targets for all samples
fastqc_targets = []
fastp_targets = []
fastp_target_display = [] # json and html
count_targets = []

for sample in SAMPLES:
    tmp_sample = re.split(r'\.fastq', sample)[0]
    sample_name = re.split(r'_R[12]\.', sample)[0]

    fastqc_files = expand("{fastqc_dir}/fastqc/{tmp_sample}_fastqc.{ext}", fastqc_dir=fastqc_dir, tmp_sample=tmp_sample, ext=["html", "zip"])
    fastp_files = expand("{fastp_dir}/{sample_name}_{read}.trimmed.fastq.gz", fastp_dir=fastp_dir, sample_name=sample_name, read=["R1", "R2"])
    fastp_reports = expand("{fastp_dir}/{sample_name}.report.{ext}", fastp_dir=fastp_dir, sample_name=sample_name, ext=["json", "html"])
    count_files = expand("{out_path}/stats/{sample_name}.qc.counts.tsv", out_path=out_path, sample_name=sample_name)

    fastqc_targets.extend(fastqc_files)
    fastp_targets.extend(fastp_files)
    fastp_target_display.extend(fastp_reports)
    count_targets.extend(count_files)


rule all:
    input:
        expand(fastqc_targets),
        expand(fastp_targets),
        expand(fastp_target_display),
        expand(count_targets)


# Define rule for quality control w/ fastqc
rule fastqc:
    input:
        fastq = os.path.join(path, "{tmp_sample}.fastq.gz")
    output:
        "{fastqc_dir}/fastqc/{tmp_sample}_fastqc.html",
        "{fastqc_dir}/fastqc/{tmp_sample}_fastqc.zip"
    threads: 32
    shell:
        '''
        #!/bin/bash

        module load fastqc/0.11.9
        fastqc --threads {threads} --outdir={fastqc_dir}/fastqc {input.fastq}
        module unload fastqc/0.11.9
        '''

# Rule to generate trimmed fastq files using fastp
rule fastp:
    input:
        read_one = os.path.join(path, "{sample_name}_R1.fastq.gz"),
        read_two = os.path.join(path, "{sample_name}_R2.fastq.gz")
    output:
        trim_one = "{fastp_dir}/{sample_name}_R1.trimmed.fastq.gz",
        trim_two = "{fastp_dir}/{sample_name}_R2.trimmed.fastq.gz",
        json = "{fastp_dir}/{sample_name}.report.json",
        html = "{fastp_dir}/{sample_name}.report.html"
    threads: 4
    shell:
        '''
        #!/bin/bash
        module load fastp/0.23.2
        fastp --in1 {input.read_one} --in2 {input.read_two} \
            --out1 {output.trim_one} --out2 {output.trim_two} \
            --thread {threads} --dedup --html {output.html} --json {output.json}
        module unload fastp/0.23.2
        '''

rule count:
    input:
        r1 = os.path.join(path, "{sample_name}_R1.fastq.gz"),
        r2 = os.path.join(path, "{sample_name}_R2.fastq.gz"),
        t1 = os.path.join(fastp_dir, "{sample_name}_R1.trimmed.fastq.gz"),
        t2 = os.path.join(fastp_dir, "{sample_name}_R2.trimmed.fastq.gz")
    output:
        txt = "{out_path}/stats/{sample_name}.qc.counts.tsv"
    shell:
        """
        #!/bin/bash

        raw1=`zcat {input.r1} | wc -l | awk '{{rd=$1 / 4; print rd}}'`
        raw2=`zcat {input.r2} | wc -l | awk '{{rd=$1 / 4; print rd}}'`
        trimmedR1=`zcat {input.t1} | wc -l | awk '{{rd=$1 / 4; print rd}}'`
        trimmedR2=`zcat {input.t2} | wc -l | awk '{{rd=$1 / 4; print rd}}'`

        sum_raw=$(expr $raw1 + $raw2)
        sum_trimmed=$(expr $trimmedR1 + $trimmedR2)

        echo -e "{wildcards.sample_name}\t${{raw1}}\t${{raw2}}\t${{trimmedR1}}\t${{trimmedR2}}\t${{sum_raw}}\t${{sum_trimmed}}" >> {output}

        """
