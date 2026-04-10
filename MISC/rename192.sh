#!/bin/bash

# Renaming Sarropoulous ATAC fragment files so the cell names match those in the Seurat object.

fragDir=/.mounts/labs/pailab/src/neurodev-genomics/scMultiome/Sarropoulos_2026

echo "Renaming sa192..."
zcat ${fragDir}/sa192_fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) {$4="sa192_Human_Cerebellum_15wpc_F#"$4; print $0}}' \
| bgzip > ${fragDir}/sa192_fragments_renamed.tsv.gz

echo "Renaming sa191..."
zcat ${fragDir}/sa191_fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) {$4="sa191_Human_Cerebellum_16wpc_F#"$4; print $0}}' \
| bgzip > ${fragDir}/sa191_fragments_renamed.tsv.gz