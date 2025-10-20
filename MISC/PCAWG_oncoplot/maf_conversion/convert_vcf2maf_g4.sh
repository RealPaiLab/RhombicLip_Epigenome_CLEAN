# conda activate vcf2maf
cat ./sample_list/g4_pcawg.sample | while read sample; do
	vcf2maf.pl --input-vcf /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/vcf/${sample}.consensus.20160830.somatic.snv_mnv.hg38.vcf --output-maf /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/g4_pcawg_maf/${sample}.consensus.20160830.somatic.snv_mnv.hg38.maf --ref-fasta /.mounts/labs/pailab/private/xsun/Database/genomeRef/hg38/hg38.fa --vep-path /.mounts/labs/pailab/private/xsun/Software/conda_env/vcf2maf/bin/ --vep-data /.mounts/labs/pailab/private/xsun/Database/.vep/ --ncbi-build GRCh38 --tumor-id ${sample} --tmp-dir /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/tmp
done
