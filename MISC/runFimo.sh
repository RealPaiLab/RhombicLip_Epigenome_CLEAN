#!/bin/bash

memePath="/home/rstudio/software/meme/bin"
tf_db=/home/rstudio/isilon/src/transcription-factors/HOCOMOCOv12/H12CORE_meme_format.meme
inFile="/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260408_DMR_overlap_EPs/inferredEP_DMRs.fa" 

$memePath/fimo --oc $inFile"_fimo" --thresh 1e-4 --verbosity 2 $tf_db  $inFile