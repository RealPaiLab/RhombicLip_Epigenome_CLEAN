#!/bin/bash

# Renaming Sarropoulous ATAC fragment files so the cell names match those in the Seurat object.



fragDir=/.mounts/labs/pailab/src/neurodev-genomics/scMultiome/Sarropoulos_2026

###echo "Renaming sa038..."
###zcat ${fragDir}/sa038_fragments.tsv.gz | \
###    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) {$4="sa038_Human_Cerebellum_11wpc_M#"$4; print $0}}' \
###| gzip > ${fragDir}/sa038_fragments_renamed.tsv.gz

###echo "Renaming sa086..."
###zcat ${fragDir}/sa086_fragments.tsv.gz | \
###    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) {$4="sa086_Human_Cerebellum_11wpc_M#"$4; print $0}}' \
###| bgzip > ${fragDir}/sa086_fragments_renamed.tsv.gz
###
###echo "Renaming sa206..."
###zcat ${fragDir}/sa206_fragments.tsv.gz | \
###    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) {$4="sa206_Human_Cerebellum_16wpc_F#"$4; print $0}}' \
###| bgzip > ${fragDir}/sa206_fragments_renamed.tsv.gz

echo "Renaming sa207..."
zcat ${fragDir}/sa207_fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) {$4="sa207_Human_Cerebellum_17wpc_M#"$4; print $0}}' \
| bgzip > ${fragDir}/sa207_fragments_renamed.tsv.gz


