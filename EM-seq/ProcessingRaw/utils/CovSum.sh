### Calculate cytosine report's 5x, 10x C percentage; median and mean coverage. Preferred to be used on C>T removed cytosine report 
ls /.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/CpG_snpFiltered/*.gz | while read line; 
	do echo $line; 
	# 5x, 10x percentage
	less $line | awk '{if ($1!="pUC19" && $1!="J02459.1") {count++; tmp=($4+$5); if (tmp>=5) five++; if (tmp>=10) ten++}} END {print five/count, ten/count}'

	# median and mean coverage
	less $line | awk '{if ($1!="pUC19" && $1!="J02459.1") print $4+$5}' | datamash median 1 mean 1
done
