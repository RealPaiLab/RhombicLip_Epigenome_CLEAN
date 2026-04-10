#!/bin/bash

tabix=/u/spai/software/htslib-1.15/tabix
fragDir=/.mounts/labs/pailab/src/neurodev-genomics/scMultiome/Sarropoulos_2026

for f in $fragDir/sa207*renamed.tsv.gz; do
    echo "Indexing $f"
    $tabix -p bed $f
done