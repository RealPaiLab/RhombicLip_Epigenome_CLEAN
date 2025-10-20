#!/bin/bash

inDir=/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/oligoDesign/stats/
inFile=${inDir}/dmr_fetalCbEnh_240711.bed
jasparFile=${inDir}/../large_regions_jaspar_240731.out

# use bedtools to extract the rows from jasparFile that overlap with inFile
bedtools intersect -a ${jasparFile} -b ${inFile} -wa -wb > ${inDir}/jaspar_overlaps.bed 