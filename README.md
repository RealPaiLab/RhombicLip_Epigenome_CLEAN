# MB_scDNAm

## EM-seq 
This folder contains all scripts for EM-seq raw data processing and downstream analyses. 

### Analysis

#### env
The environment settings for the analyses.

#### FETHB2-FETHB3
The EM-seq of the human fetal rhombic lip (RL) was performed in two batches (batch 1: FETHB2; batch 2: FETHB3). Scripts under this folder pool two batches together for downstream analyses. 

- **CpG** (Analyses focus on CpG sites.)
  - ***DMR***
    - 01b.identify_DMR_CTsnvExcluded_withoutBatchCorrection.R: `Rscript --no-save 01b.identify_DMR_CTsnvExcluded_withoutBatchCorrection.R` identifies DMRs (Differentially methylated regions) using DSS package, by first test for each CpG sites (hierarchical Bayesian model) and then merge significant sites into regions.
    - 01c.identify_DMR_CTsnvExcluded_withoutBatchCorrection_MethylKit.R: `Rscript --no-save 01c.identify_DMR_CTsnvExcluded_withoutBatchCorrection_MethylKit.R` identifies DMRs using MethylKit package. This method tests DMRs based on provided regions. I spliced the genome into 50 bp non-overlapping windows as it is the smallest size limit applied in DSS. Only windows covered in all samples were kept for DMR detection. This was used as an alternative validation of the DSS identified DMRs; therefore, I intersected the MethylKit DMRs and DSS DMRs to check overlap and to prioritize DMRs.
    - 02b.annotate_DMR_CTsnvExcluded_withoutBatchCorrection.R: `Rscript --no-save 02b.annotate_DMR_CTsnvExcluded_withoutBatchCorrection.R` annotates DMRs identified using DSS method. E.g. ENCODE cCREs, human fetal hindbrain ChIP-seq peaks, and transposable elements based on repeatMasker.
    - 03b.AME_CTsnvExcluded_withoutBatchCorrection.R: `Rscript --no-save 03b.AME_CTsnvExcluded_withoutBatchCorrection.R` runs MEME::AME on the DSS DMRs to find core TFs (transcription factors) potentially interacting with the DMRs based on HOCOMOCO v11 database. To reduce the searching space of motifs, only motifs of TFs highly expressed in RL-VZ/SVZ single cells were used (Details of highly expressed genes identification code can be found in MB_scDNAm/MISC/RL_genes section).
    - 04a.check_DMR_CTsnvExcluded_withoutBatchCorrection.R: `Rscript --no-save 04a.check_DMR_CTsnvExcluded_withoutBatchCorrection.R` plots violin plots of DNA methyaltion between RL-VZ/SVZ DMRs that are overlap promoter region of RL-VZ/SVZ DEGs.
    - 05a.compare_DMR_RNA_CTsnvExcluded_withoutBatchCorrection.R: `Rscript --no-save 05a.compare_DMR_RNA_CTsnvExcluded_withoutBatchCorrection.R` plots scatter plots comparing DMRs overlapping promoters/enhancers of RL-VZ/SVZ DEGs, and calculates Pearson correlation. 
    - archived: all outdated/not used scripts
  - ***DMV*** (Differentially methylated valleys)
    - 01a.methylKit_segmentation.R: `Rscript --no-save 01a.methylKit_segmentation.R` uses MethylKit package to separate the genome into sections based on CpG methylation status cut points for each sample.
    - 02a.DNAmValley.R: `Rscript --no-save 02a.DNAmValley.R` takes section information generated from the 01a.methylKit_segmentation.R and filter for DNA methylation valleys. Then, the valleys were interesected with the DSS DMRs to robustly check overall methylation status difference between RL-VZ and RL-SVZ. To gain a consensus valley map, valleys from all samples were merged when any overlap occurs.
  - ***overlapEnrichment***:
    to test for the enrichment of functional genomic regions in the DMRs using a permutation method (modified based on code from Shraddha Pai).
    - `Rscript --no-save DMROverlap_Stats.R` will test the enrichment of different functionally annotated regions in the given DMRs and generate related plots.
  - ***VariantRemoval***:
    code used for re-extracting cytosine report with potential C>T variants filtration. (Could be incorporated into the formal PrcessingRaw pipeline); **TODO: this removes sample-level heterozygous C>T variants but not necessarily all variants in other samples or the larger population. Should still consider gnomAD (together with the current method) for a complete removal?
- **CpH** (Analyses focus on non-CpG sites (e.g. CpA))
  - UNFINISHED. **TODO: resolve CpH loading error.

### Figures
Gather and generate figures used for publication purpose. `Rscript --no-save plot.R` **TODO: permutation test figures

### ProcessingRaw
Snakemake pipeline that uses OICR HPC to process EM-seq raw data to generate cytosine reports (based on Shraddha and Paul Wambo code). Use **config.yaml** under config to specify basic paramters and input/output. **env/hpc_ubuntu20_modulator** contains all modulator yaml files for OICR hpc. Use **run.sh** to adjust HPC resource and `bash run.sh` to initiate the EM-seq pipeline. Detailed sakemake rules can be found under the workflow folder. **TODO: improve pipeline stat summary and cytosine report extraction rules (1. spike-in removal; 2. CpG/CpH/CpA report; 3. C>T removal)

## MISC
This folder contains scripts for different tasks not directly for DNA methylation analyses.

### JASPAR_getTFBS
To get TFBS of given regions from JASPAR db. Adopted from Shraddha.

### RL_genes
  - outDated: previous scripts for TF db preparation
  - HOCOMOCOv12: Filter out TFs not in the provided gene list and generate new .meme format database based on HOCOMOCO v12 CORE motif database.
    - Hendrikse_2022.R: `Rscript --no-save Hendrikse_2022.R` identifies genes expressed above certain level (>1% cells in the specified cell type/cluster) using snRNA-seq of RL VZ/SVZ (Hendrikse et al. 2022).
    - filter_tf.py: `python3 filter_tf.py` reads provided gene lists and keep only motifs of TFs that are within both the gene list and HOCOMOCO v12 core TF motifs. The script accounts for discrepancy between gene name and TF name by referencing official H12CORE_annotation.jsonl. 
    - main.sh: `bash main.sh` runs Hendrikse_2022.R and filter_tf.py to generate RL active TF motif database used in /EM-seq/FETHB2-FETHB3/CpG/DMR/03b.AME_CTsnvExcluded_withoutBatchCorrection.R
  
### PCAWG_oncoplot
Use PCAWG hg19 SNV/INDEL/CNA information and RL-VZ/SVZ DMRs to generate oncoplot. Please see README.md files attached for each analysis. GISTIC2.0 for only Grp3 or Grp4 tumours were run using GenePattern online platform, output and description of params used can be found under /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/GenePattern.

## MPRA
For design and analyses of DMR validation using massively parallel reporter assay (MPRA).

### oligo_design
  - 01b.target_trim.R: trim DMRs into fragments of target size for MPRA (180bp in this case).
  - 02a.target_ranking.R: score and rank trimmed fragments by defined score matrix (currently for pilot; please check and modify for each case).
  - 03a.additional_targets.R: identifying addtional customized interested targets.
  - 04a.generate_targetProbes.R: generate target sequence for probe design. Note that there is one specific pair of oligos included in the script, which are USH2A DMRs (ol HAR) and its mutated sequence with a SNV (found in Grp3 MB).
  - 05a.positiveControl_ranking.R: generate and rank a list of postive controls for MPRA
  - 06a.negativeControl_ranking.R: generate and rank negative controls based on targets, using random shuffling.
  - 07a.summariseOligos.R: annotate the oligos with fetal developmental enhancers, ChIP-seq peaks (fCB, N2/N3), hars, neurodev/mb gene, snv/indel/sv ... 
  - 08a.generateFinalOligoPool.R: generate the final oligo pool and format the file for Agilent (3-col format). 
  - tracks/createTracks.R: create UCSC track hub tracks using the MPRA related information.
    
## Smart-RRBS
This folder contains code for Smart-RRBS data analysis of one MB tumour (194 cells). Stopped after EDA because of data sparsity even after merging CpGs by large windows (10-50kb).
