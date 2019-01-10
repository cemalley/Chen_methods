#!/bin/bash

module load GATK
module load bwa
module load picard
module load R
module load samtools
module load annovar

REF='/fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta'

ID='CT01'

cd /data/NCATS_ifx/data/WES/

#bwa mem -t 8 -R "@RG\tID:CT01\tPL:Illumina\tSM:${ID}" /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa <(zcat fastq/${ID}_S1_R1_001.fastq.gz) <(zcat fastq/${ID}_S1_R2_001.fastq.gz) | samtools view -bS - | samtools sort - > ${ID}/${ID}.bam; samtools index ${ID}/${ID}.bam

# java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
#   INPUT=${ID}/${ID}.bam \
#   OUTPUT=${ID}/${ID}_dedup.bam \
#   METRICS_FILE=metrics.txt

# java -jar picard.jar BuildBamIndex \
#     INPUT=${ID}/${ID}_dedup.bam
#

GATK -m 20g HaplotypeCaller \
  -R $REF \
  -I ${ID}/${ID}_dedup.bam \
  -o vcf/${ID}_raw_variants.vcf \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30
  -L /data/NCATS_ifx/data/WES/hglft_genome_1845_3870.bed

GATK -m 20g SelectVariants \
-R $REF \
-V $V \
-selectType SNP \
-o ${ID}_SNP_raw.vcf
  -L /data/NCATS_ifx/data/WES/hglft_genome_1845_3870.bed

GATK -m 20g SelectVariants \
-R $REF \
-V $V \
-selectType INDEL \
-o ${ID}_INDEL_raw.vcf
  -L /data/NCATS_ifx/data/WES/hglft_genome_1845_3870.bed

# SNPs
java -jar $GATK_JAR \
-T VariantFiltration \
-R $REF \
-V $V \
-filterName \"QD_filter\" \
-filter \"QD < 2.0\" \
-filterName \"FS_filter\" \
-filter \"FS > 60.0\" \
-filterName \"MQ_filter\" \
-filter \"MQ < 40.0\" \
-filterName \"SOR_filter\" \
-filter \"SOR > 4.0\" \
-o $O

# indels
java -jar $GATK_JAR \
-T VariantFiltration \
-R $REF \
-V $V \
-filterName \"QD_filter\" \
-filter \"QD < 2.0\" \
-filterName \"FS_filter\" \
-filter \"FS > 200.0\" \
-filterName \"SOR_filter\" \
-filter \"SOR > 10.0\" \
-o $O


# sbatch --mem=30g --gres=lscratch:30 --time=24:00:00 scripts/GATK.sh
# 14173414

# then BQSR
# then analyze covariates
# then apply BQSR
# then parse metrics
# split SNPs and indels
# filter variants
# then run annovar
