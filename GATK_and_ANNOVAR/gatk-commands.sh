docker run -v /Volumes/ncatssctl/NGS_related/:/gatk/my_data -it broadinstitute/gatk

docker pull broadinstitute/gatk

docker run -it broadinstitute/gatk


java -jar gatk-package-4.0.11.0-local.jar \
    -T HaplotypeCaller \
    -R reference.fa \
    -I preprocessed_reads.bam \
    -L 20 \
    --genotyping_mode DISCOVERY \
    -stand_emit_conf 10 \
    -stand_call_conf 30 \
    -o raw_variants.vcf


GATK -m 7g RealignerTargetCreator \
  -R ref.fasta \
  -I input1.bam \
  -o output1.intervals \
  --known /fdb/GATK_resource_bundle/hg19/1000G_phase1.indels.hg19.vcf.gz
GATK -m 7g RealignerTargetCreator \
  -R ref.fasta \
  -I input2.bam \
  -o output2.intervals \
  --known /fdb/GATK_resource_bundle/hg19/1000G_phase1.indels.hg19.vcf.gz


GATK -m 7g HaplotypeCaller \
  -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
  -I /data/NCATS_ifx/data/WES/CT01/CT01.bam \
  -o /data/NCATS_ifx/data/WES/vcf/CT01.vcf \
  -L 20 \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30
