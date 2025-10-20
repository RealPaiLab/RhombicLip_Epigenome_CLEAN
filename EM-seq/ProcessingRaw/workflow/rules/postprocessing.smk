# Get config file
configfile: f"{workflow.basedir}/../config/config.yaml"

SAMPLES = []

path: str = config["directory"]["input"]

for fastq in os.scandir(path):
    if fastq.name[0] != '.' and 'R1' in fastq.name:
        SAMPLES.append(fastq.name)

# Define input and output file paths
out_path: str = config["directory"]["output"]
bwa_dir = f"{out_path}/alignment/bwa"
dackel_dir = f"{out_path}/alignment/methyldackel"

picard_targets = []
dackel_targets = []

for sample in SAMPLES:
    sample_name = re.split(r'_R[12]\.', sample)[0]

    dackel = expand("{dackel_dir}/plot/{sample_name}_OB.svg", dackel_dir=dackel_dir, sample_name=sample_name)
    mbias_files = expand("{dackel_dir}/plot/{sample_name}.combined_mbias.tsv", dackel_dir=dackel_dir, sample_name=sample_name)
    dackel_files = expand("{dackel_dir}/report/{sample_name}.cytosine_report.{ext}", dackel_dir=dackel_dir, sample_name=sample_name, ext=["txt", "txt.gz", "txt.gz.tbi"])
    reports = expand("{out_path}/stats/report/{sample_name}.cytosine.controls.tsv", out_path=out_path, sample_name=sample_name)

    dackel_targets.extend(mbias_files),
    dackel_targets.extend(dackel_files)
    dackel_targets.extend(dackel)
    dackel_targets.extend(reports)


rule all:
    input:
        expand(dackel_targets)


rule mbias:
    input:
        bam = os.path.join(bwa_dir, "sorted/{sample_name}.filtered.bam"),
    output:
        ob  = "{dackel_dir}/plot/{sample_name}_OB.svg",
        tsv = "{dackel_dir}/plot/{sample_name}.combined_mbias.tsv"
    threads: 16,
    params:
        chrs=("chr20 chr21 chr22"),
        GENOME_FA = "/.mounts/labs/pailab/src/aligner_indexes/hg38/bwameth/genome_lambda_puc19.fa",
        plot_prefix = "{dackel_dir}/plot/{sample_name}"
    log:
        "{dackel_dir}/plot/{sample_name}.mbias.log"
    shell:
        '''
        #!/bin/bash
        module load methyldackel/0.6.1

        for chr in {params.chrs}; do

            join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <(MethylDackel mbias --noSVG -@ {threads} -r $chr \
                {params.GENOME_FA} {input.bam} \
                | tail -n +2 | awk '{{print $1"-"$2"-"$3"\t"$0}}' | sort -k 1b,1) \
                <(MethylDackel mbias --noSVG --keepDupes -F 2816 -@ {threads} -r $chr \
                {params.GENOME_FA} {output.tsv} |  tail -n +2 | awk '{{print $1"-"$2"-"$3"\t"$0}}' | sort -k 1b,1) \
                        | sed "s/^/$chr\tCG\t/" >> {output.tsv}
            
            MethylDackel mbias -@ {threads} -r $chr {params.GENOME_FA} {input.bam} {params.plot_prefix} 2>> {log}
            
        done

        module unload methyldackel/0.6.1
        '''

rule cys_report:
    input:
        bam = os.path.join(bwa_dir, "sorted/{sample_name}.filtered.bam")
    output:
        txt = temp("{dackel_dir}/report/{sample_name}.cytosine_report.txt"),
        bgzip = "{dackel_dir}/report/{sample_name}.cytosine_report.txt.gz",
        tbi = "{dackel_dir}/report/{sample_name}.cytosine_report.txt.gz.tbi"
    threads: 16,
    params:
        GENOME_FA = "/.mounts/labs/pailab/src/aligner_indexes/hg38/bwameth/genome_lambda_puc19.fa",
        sample_name_prefix = "{dackel_dir}/report/{sample_name}",
        OT = config['boundaries']['OT'],
        OB = config['boundaries']['OB']
    shell:
        '''
        #!/bin/bash
        module load methyldackel/0.6.1

        MethylDackel extract --cytosine_report --OT {params.OT} --OB {params.OB} -@ {threads} {params.GENOME_FA} {input.bam} -o {params.sample_name_prefix}
        bgzip -c {output.txt} > {output.bgzip} --threads {threads}
        tabix -s 1 -b 2 -e 2 -S 1 {output.bgzip}

        module unload methyldackel/0.6.1
        '''

rule controls:
    input:
        report = os.path.join(dackel_dir, "report/{sample_name}.cytosine_report.txt")
    output:
        txt = "{out_path}/stats/report/{sample_name}.cytosine.controls.tsv"
    params:
        value = "{sample_name}"
    shell:
        '''
        #!/bin/bash

        puc=$(cat {input.report} | grep "pUC19" | awk -v file={params.value} '
            BEGIN {{
                c_count=0;
                t_count=0;
            }}
            {{
                c_count+=$4;
                t_count+=$5;
            }}
            END {{
                ncrate=c_count/(c_count+t_count);
                crate=1-ncrate;
                print crate;
            }}
        ')

        lambda=$(cat {input.report} | grep "J02459.1" | awk -v file={params.value} '
            BEGIN {{
                c_count=0;
                t_count=0;
            }}
            {{
                c_count+=$4;
                t_count+=$5;
            }}
            END {{
                ncrate=c_count/(c_count+t_count);
                crate=1-ncrate;
                print crate;
            }}
        ')

        echo -e "{params.value}\t{{$puc}}\t{{$lambda}}" > {output}

        '''
