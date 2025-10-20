### Config
baseDir=/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240522

# Input config
dmrFile=${baseDir}/DMRs.csv
keoughHARs=/.mounts/labs/pailab/private/xsun/output/HARs/Keough_2023_HARs/HARs_hg38.bed
mergedHARs=/.mounts/labs/pailab/private/xsun/output/HARs/multipleLabs_merged_HARs/nchaes_merged_hg38.bed

# Output config
KeoughIntersect=${baseDir}/DMR_PollardLab-Keough-HARs_intersect.bed
MergedIntersect=${baseDir}/DMR_PollardLab-Merged-HARs_intersect.bed

### Intersect RL-VZ/SVZ DMRs and HARs
# Intersect with Keough_2023 HARs
awk 'BEGIN {OFS="\t"} NR>1 {print}' ${dmrFile} | bedtools intersect -a - -b ${keoughHARs} -wo > ${KeoughIntersect}

# Intersect with Merged HARs from Pollard lab and others
awk 'BEGIN {OFS="\t"} NR>1 {print}' ${dmrFile} | bedtools intersect -a - -b ${mergedHARs} -wo > ${MergedIntersect}

