### extract CHN for non-CpG analysis (should variants also be excluded? This is not only C>T; maybe add variant identification params)

## CONFIG
spike_in=/.mounts/labs/pailab/src/aligner_indexes/hg38/bwameth/genome_lambda_puc19.fa

out_dir=/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/CpH

fethb2_bwa_dir=/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB2/output/alignment/bwa

fethb3_bwa_dir=/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/bwa/sorted

## GET FILES
files=$(find ${fethb2_bwa_dir} ${fethb3_bwa_dir} -type f -name "*filtered.bam")

## MAIN 
for file in ${files}; do 
	echo "--- Processing: "$file 
	filename=$(basename ${file})

	# Extract
	echo "--- Extracting CHN cytosime report."
	MethylDackel extract \
	 --cytosine_report \
	 --OT 0,0,0,147 \
	 --OB 3,0,5,0 \
	 -@ 16 \
	 ${spike_in} \
	 ${file} \
	 -o ${out_dir}/${filename%.filtered.bam}.CHN \
	 --CHH \
	 --CHG

	 # Compress
	 echo "--- Compressing."
	 bgzip -c ${out_dir}/${filename%.filtered.bam}.CHN.cytosine_report.txt > ${out_dir}/${filename%.filtered.bam}.CHN.cytosine_report.txt.gz

	 # Tabix
	 echo "--- Creating tabix."
	 tabix -s 1 -b 2 -e 2 -S 1 ${out_dir}/${filename%.filtered.bam}.CHN.cytosine_report.txt.gz

	 echo "--- Finished."
done







