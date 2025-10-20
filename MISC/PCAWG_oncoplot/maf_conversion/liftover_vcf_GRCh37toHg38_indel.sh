# conda activate picard
cat ./sample_list/g4_pcawg.sample ./sample_list/g3_pcawg.sample | while read sample; do
    picard -Xmx8g LiftoverVcf -I /.mounts/labs/pailab/src/MB_genomics/WGS/PCAWG/final_consensus_snv_indel_passonly_icgc.public/indel/${sample}.consensus.20161006.somatic.indel.vcf.gz -O /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/vcf/${sample}.consensus.20160830.somatic.indel.hg38.vcf -CHAIN /.mounts/labs/pailab/private/xsun/Database/liftover/GRCh37toHg38.over.chain -R /.mounts/labs/pailab/private/xsun/Database/genomeRef/hg38/hg38.fa -REJECT /.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/vcf/${sample}.consensus.20160830.somatic.indel.hg38.unlifted.vcf
done

echo "--- Finished liftover"

