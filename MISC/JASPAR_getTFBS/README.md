
This folder contains code to extract TFBS from JASPAR for user-provided coordinates.
Note that it only returns TFBS with score of 400.
For details see: https://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=290259872_tndyJlqyvi4iWtlWaXIXbDZqvILC&db=hg19&c=chr6&g=jaspar

To run this code in R: 

setwd("../JASPAR_getTFBS") # cd to the current directory in R
source("getTFBSMotifMatrix.R")
mat <- getTFBSMotifMatrix(gr) # where gr is your interval set of interest. 


Source of code:
* bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html 
* bigBedToBed: https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
* extract_TFBSs_JASPAR: https://bitbucket.org/CBGR/jaspar_tfbs_extraction/src/main/