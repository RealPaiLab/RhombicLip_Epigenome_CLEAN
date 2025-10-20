

# Get config file
configfile: f"{workflow.basedir}/../config/config.yaml"

SAMPLES = []

path: str = config["directory"]["input"]

for fastq in os.scandir(path):
    if fastq.name[0] != '.':
        SAMPLES.append(fastq.name)

# Define input and output file paths
out_path: str = config["directory"]["output"]
fastp_dir = f"{out_path}/prealignment/trim"
align_dir = f"{out_path}/alignment"
bwa_dir = f"{out_path}/alignment/bwa"

bam_targets = []
sam_targets = []
picard_targets = []
count_targets = []


for sample in SAMPLES:
    sample_name = re.split(r'_R[12]\.', sample)[0]

    sam_files = expand("{bwa_dir}/unsorted/{sample_name}.sam", sample_name=sample_name, bwa_dir=bwa_dir)
    bam_files = expand("{bwa_dir}/sorted/{sample_name}.aligned.sorted.bam", sample_name=sample_name, bwa_dir=bwa_dir)
    picard_files = expand("{bwa_dir}/sorted/{sample_name}.{ext}", sample_name=sample_name, bwa_dir=bwa_dir, ext=['filtered.bam', 'filtered.metrics.txt'])
    count_files = expand("{out_path}/stats/alignment/{sample_name}.alignment.counts.tsv", out_path=out_path, sample_name=sample_name)

    bam_targets.extend(bam_files)
    sam_targets.extend(sam_files)
    picard_targets.extend(picard_files)
    count_targets.extend(count_files)


rule all:
    input:
        expand(sam_targets),
        expand(bam_targets),
        expand(picard_targets),
        expand(count_targets)        



# Define rule to align reads to reference genome using BWA-Meth
rule align_reads:
    input:
        read_one = os.path.join(fastp_dir, "{sample_name}_R1.trimmed.fastq.gz"),
        read_two = os.path.join(fastp_dir, "{sample_name}_R2.trimmed.fastq.gz"),
    params:
       "/.mounts/labs/pailab/src/aligner_indexes/hg38/bwameth/genome_lambda_puc19.fa"
    output:
        pe = temp("{bwa_dir}/unsorted/{sample_name}.sam")
    log:
        "{bwa_dir}/unsorted/{sample_name}.bwameth.log"
    threads: 30
    shell:
        '''
        #!/bin/bash

        module load bwa-meth/0.2.7

        bwameth.py --threads {threads} --reference {params} {input.read_one} {input.read_two} > {output} 2> {log}

        module unload bwa-meth/0.2.7
        '''

# Define rule to convert SAM files to BAM files using samtools
rule sam_to_bam:
    input:
        sam = os.path.join(bwa_dir, "unsorted/{sample_name}.sam")
    output:
        bam = temp("{bwa_dir}/sorted/{sample_name}.aligned.sorted.bam")
    threads: 15
    shell:
        '''
        #!/bin/bash

        module load samtools/1.15

        samtools view -t 2 -Sb {input.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}

        module unload samtools/1.15

        '''


rule mark_duplicates:
    input:
        bamfiles = os.path.join(bwa_dir, "sorted/{sample_name}.aligned.sorted.bam")
    output:
        filtered = "{bwa_dir}/sorted/{sample_name}.filtered.bam",
        stats = "{bwa_dir}/sorted/{sample_name}.filtered.metrics.txt"
    shell:
        '''
        #!/bin/bash

        module load picard-tools/1.89 picard/2.27.5

        picard=/.mounts/labs/pailab/modulator/sw/Ubuntu18.04/picard-2.27.5/picard.jar

        java -Xmx4g -Xms4g -jar $picard MarkDuplicates \
            I={input.bamfiles} O={output.filtered} M={output.stats} \
            ASSUME_SORTED=false REMOVE_DUPLICATES=true CREATE_INDEX=true

        module unload picard-tools/1.89 picard/2.27.5
        '''


rule count:
    input:
        bamfile = os.path.join(bwa_dir, "sorted/{sample_name}.aligned.sorted.bam"),
        filtered_bam = os.path.join(bwa_dir, "sorted/{sample_name}.filtered.bam")
    output:
        txt = "{out_path}/stats/alignment/{sample_name}.alignment.counts.tsv",
    params:
        value = "{sample_name}",
    shell:
        """
        #!/bin/bash

        module load samtools/1.15

        num_aligned=$(samtools flagstat {input.bamfile} | head -n 12 | tail -1 | awk '{{print$1;}}')
        dedup_count=$(samtools flagstat {input.filtered_bam}| grep "properly paired" | awk '{{print $1;}}')

        echo -e "{params.value}\t${{num_aligned}}\t${{dedup_count}}" >> {output}

        module unload samtools/1.15

        """
