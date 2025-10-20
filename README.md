This repo contains analysis to reproduce the results and figures from the following preprint:

**Epigenomic landscape of the developing human rhombic lip reveals gene regulatory network and non-coding loci of developmental, evolutionary, and disease relevance**
Sun X., Menon S., Wambo P., Lungu I., Birth Defects Research Laboratory, Aldinger K.A., and S. Pai

* EM-seq preprocessing: `ProcessingRaw/` folder
* Exploratory data analysis: `EDA/` folder
* DMR analysis: `FETHB2-FETHB3/CpG/DMR/`
  * DMR calls: `01.identify_DMR_CTsnvExcluded_withoutBatchCorrection.R`
  * Plot sample DMRs & correlate with transcription: `01b.plotSampleDMRs_correlateWithRNA.R`
  * Annotate DMR by ENCODE cCREs `02.annotate_DMR_CTsnvExcluded_withoutBatchCorrection.R`
  * Annotate DMR by CpG island status: `02b.DMRs_CpGIslandShores.R`
  * Run TFBS enrichment: `03.AME_CTsnvExcluded_withoutBatchCorrection.R`
  * Plot TFBS enrichment results: `03b.plotTFBSHits.R`
  * Annotate DMRs with putative target genes and their neurodevelopmental or cancer relevance: `04.DMR_link2Genes_ABC.R`
  * Annotate DMRs with putative target genes using Mannens et al. fetal brain single-cell multiome dataset: `04c.linkDMRsToGenes.R`
  * Create GRN by inferring TF-enhancer-target gene triplets: `04b.DMR_link2Genes_AnnotateDMRs.R`
  * Annotate superenhancers with nearest target genes: `05.annotateSuperEnhancers.R`
