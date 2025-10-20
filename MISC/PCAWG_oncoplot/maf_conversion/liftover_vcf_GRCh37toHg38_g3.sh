# conda activate picard
cat ./sample_list/g3_pcawg.sample | while read sample; do
    picard -Xmx8g LiftoverVcf -I /.mounts/labs/pailab/src/MB_genomics/WGS/PCAWG/final_consensus_snv_indel_passonly_icgc.public/snv_mnv/${sample}.consensus.20160830.somatic.snv_mnv.vcf.gz -O /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/vcf/${sample}.consensus.20160830.somatic.snv_mnv.hg38.vcf -CHAIN /.mounts/labs/pailab/private/xsun/Database/liftover/GRCh37toHg38.over.chain -R /.mounts/labs/pailab/private/xsun/Database/genomeRef/hg38/hg38.fa -REJECT /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/vcf/${sample}.consensus.20160830.somatic.snv_mnv.hg38.unlifted.vcf
done

