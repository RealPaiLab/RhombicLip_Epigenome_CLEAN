#!/bin/bash

module load snakemake/6.4.0
username=$1
configfile=$2

### Run preprocessing Snakemake ###
echo "### Running pre-processing step ###"
snakemake \
    --jobs 150 \
    --cluster "qsub -P pailab -V -cwd -b y -N emseq_preprocessing -M ${username}@oicr.on.ca -m a -l h_rt=5000,h_vmem=4G -pe smp 32" \
    --snakefile ./workflow/Snakefile_preprocessing all \
    --configfile ${configfile} \
    --stats emseq_preprocessing_stats.json \
    --latency-wait 900

### Run processing Snakemake ###
echo "### Running processing step ###"
snakemake \
    --jobs 100 \
    --cluster "qsub -P pailab -q all.q -V -cwd -b y -N emseq_processing -M ${username}@oicr.on.ca -m ea -l h_rt=2:0:0:0,h_vmem=10G -pe smp 32" \
    --snakefile ./workflow/Snakefile_processing all --rerun-incomplete \
    --configfile ${configfile} \
    --stats emseq_processing_statistics.json \
    --latency-wait 1800


### Run postprocessing Snakemake ###
echo "### Running post-processing step ###"
snakemake \
    --jobs 100 \
    --cluster "qsub -P pailab -V -cwd -b y -N emseq_postprocessing -M ${username}@oicr.on.ca -m ea -l h_rt=2:0:0:0,h_vmem=10G -pe smp 32" \
    --snakefile ./workflow/Snakefile_postprocessing all --rerun-incomplete \
    --configfile ${configfile} \
    --stats emseq_pipeline_statistics.json \
    --latency-wait 1800

### Clean up ###
# rm -rf .snakemake
module unload snakemake/6.4.0

echo "### Finished ###"

