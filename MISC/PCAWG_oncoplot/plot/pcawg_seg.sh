subgroups=(Group3 Group4)
output=/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/pcawg_bed

for subgroup in "${subgroups[@]}"; do
	# extract g3/4 sample id
	cut -f 4 /.mounts/labs/pailab/private/xsun/output/ncMutMB/20240314/mutation/PCAWG/merged/${subgroup}_PCAWG_snv-indel_hg38.bed | sort -u > ${output}/${subgroup}_pcawg.sample
	
	# extract g3/4 mb copy number seg
	zless /.mounts/labs/pailab/src/MB_genomics/WGS/PCAWG/consensus_cnv/GISTIC_analysis/focal_input.rmcnv.pt_170207.seg.txt.gz | grep -wFf ${output}/${subgroup}_pcawg.sample - > ${output}/${subgroup}_pcawg.hg19.seg
	
	# reformat
	awk 'BEGIN {OFS="\t"} {print "chr"$2,$3,$4,$1,$6"&"$5}' ${output}/${subgroup}_pcawg.hg19.seg > ${output}/${subgroup}_pcawg.hg19.seg.bed
	
	# liftOver to hg38
	liftOver ${output}/${subgroup}_pcawg.hg19.seg.bed /Users/xsun/db/liftOver/hg19ToHg38.over.chain.gz ${output}/${subgroup}_pcawg.hg38.seg.tmp.bed ${output}/${subgroup}_pcawg.hg38.seg.tmp.unlifted
	sed "s/&/\t/g" ${output}/${subgroup}_pcawg.hg38.seg.tmp.bed | cut -f1-5 > ${output}/${subgroup}_pcawg.hg38.seg.bed
	sed "s/&/\t/g" ${output}/${subgroup}_pcawg.hg38.seg.tmp.bed | sed "s/chr//g" | awk '{OFS="\t"} {print $4,$1,$2,$3,$6,$5}' > ${output}/${subgroup}_pcawg.hg38.seg

	# extract copy number change segs
	awk '{if ($5 >= 0.5 || $5 <= -0.5) print $0}' ${output}/${subgroup}_pcawg.hg38.seg.bed > ${output}/${subgroup}_pcawg.hg38.seg.amp_del.bed
	
	awk '{if ($5 >= 0.5) print $0}' ${output}/${subgroup}_pcawg.hg38.seg.bed > ${output}/${subgroup}_pcawg.hg38.seg.amp.bed

	awk '{if ($5 <= -0.5) print $0}' ${output}/${subgroup}_pcawg.hg38.seg.bed > ${output}/${subgroup}_pcawg.hg38.seg.del.bed

	# count dmr cna
	bedtools intersect -a /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/dmrs.bed -b ${output}/${subgroup}_pcawg.hg38.seg.amp.bed -wo | cut -f1-3,7 | sort -u | cut -f1-3 | sort | uniq -c | sort -k1n  > ${output}/${subgroup}_dmr_cna.hg38.amp.count
	bedtools intersect -a /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/dmrs.bed -b ${output}/${subgroup}_pcawg.hg38.seg.del.bed -wo | cut -f1-3,7 | sort -u | cut -f1-3 | sort | uniq -c | sort -k1n  > ${output}/${subgroup}_dmr_cna.hg38.del.count
done

# combine g3/4
cat ${output}/Group3_pcawg.hg38.seg.amp.bed ${output}/Group4_pcawg.hg38.seg.amp.bed | bedtools intersect -a /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/dmrs.bed -b - -wo | cut -f1-3,7 | sort -u | cut -f1-3 | sort | uniq -c | sort -k1n  > ${output}/Group34_dmr_cna.hg38.amp.count
cat ${output}/Group3_pcawg.hg38.seg.del.bed ${output}/Group4_pcawg.hg38.seg.del.bed | bedtools intersect -a /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/dmrs.bed -b - -wo | cut -f1-3,7 | sort -u | cut -f1-3 | sort | uniq -c | sort -k1n  > ${output}/Group34_dmr_cna.hg38.del.count



