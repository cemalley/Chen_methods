library(data.table)

# Whole exome sequencing analysis: GATK Best Practices pipeline + ANNOVAR annotation
# Generate swarm files for batch running on NIH Biowulf HPC
# Claire Malley 2018
# NCATS NIH

samples <- c('CT01','CT02','CT03','CT04', 'CT05','CT06')

# GATK reference genome hg38----
ref <- '/fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta'
dbsnp <- '/fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz'
indels <- '/fdb/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
limit <- '/data/NCATS_ifx/data/WES/hglft_genome_subset.bed'

# Align----
for (i in samples){
  cat('bwa mem -t 8 -R "@RG\\tID:CT05\\tPL:Illumina\\tSM:',i,'" /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa <(zcat fastq/',i,'_R1_001.fastq.gz) <(zcat fastq/',i,'_R2_001.fastq.gz) | samtools view -bS - | samtools sort - > ',i,'/',i,'.bam; samtools index ',i,'/',i,'.bam', sep='')
  cat('\n\n')
}

# Mark duplicates----
for (i in samples){
  cat('cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates INPUT=',i,'/',i,'.bam OUTPUT=',i,'/',i,'_dedup.bam METRICS_FILE=metrics_',i,'.txt', sep='')
  cat('\n\n')
}

# Base quality score recalibration----

for (i in samples){
  cat('cd /data/NCATS_ifx/data/WES/ && GATK -m 15g BaseRecalibrator -R ',ref,' -I ',i,'/',i,'_dedup.bam -knownSites ',dbsnp,' -knownSites ',indels,' -o ',i,'/',i,'_recal_data.table -L ',limit,' && GATK PrintReads -R ',ref,' -I ',i,'/',i,'_dedup.bam -BQSR ',i,'/',i,'_recal_data.table -o ',i,'/',i,'_recal.bam -L ',limit, sep='')
  cat('\n\n')
}


# Caller----
for (i in samples){
  cat('cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R ',ref,' -I ',i,'/',i,'_recal.bam -o vcf/',i,'_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp ',dbsnp,' -L ',limit, sep='')
  cat('\n\n')
}

# Split to SNPs or Indels----

for (i in samples){
  cat('cd /data/NCATS_ifx/data/WES/ && GATK -m 15g SelectVariants -R ',ref,' -V vcf/',i,'_raw_variants.vcf -selectType SNP -o vcf/',i,'_SNP_raw.vcf -L ',limit,' && GATK -m 15g SelectVariants -R ',ref,' -V vcf/',i,'_raw_variants.vcf -selectType INDEL -o vcf/',i,'_INDEL_raw.vcf -L ',limit,sep='')
  cat('\n\n')
}

# Filter variants----

# Strict filter thresholds:
# SNPS: QD < 2.0, FS > 60.0, MQ < 40.0, MQRankSum < -12.5, ReadPosRankSum < -8.0
# INDELS: QD < 2.0, FS > 200.0, ReadPosRankSum < -20.0

# Relaxed filter thresholds:
# SNPS: QD < 2.0, FS > 60.0, MQ < 40.0, SOR > 4.0
# INDELS: QD < 2.0, FS > 200.0, SOR > 10.0


for (i in samples){
  cat('cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R ',ref,' -V vcf/',i,'_SNP_raw.vcf -filterName \"QD\" --filterExpression \"(QD < 2.0)\" -filterName \"FS\" --filterExpression \"(FS > 60.0)\" -filterName \"MQ\" --filterExpression \"(MQ < 40.0)\" -o vcf/',i,'_SNP_filtered2.vcf ', sep='')
  cat('\n\n')
  cat('cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R ',ref,' -V vcf/',i,'_INDEL_raw.vcf -filterName \"QD\" --filterExpression \"(QD < 2.0)\" -filterName \"FS\" --filterExpression \"(FS > 200.0)\" -filterName \"SOR\" --filterExpression \"(SOR > 10.0)\" -o vcf/',i,'_INDEL_filtered2.vcf ', sep='')
  cat('\n\n')
}

# merge snp and indel files----

for (i in samples){
  cat('bgzip -c ',i,'_SNP_filtered2.vcf > ',i,'_SNP_filtered2.vcf.gz ; tabix -p vcf ',i,'_SNP_filtered2.vcf.gz', sep='')
  cat('\n\n')
  cat('bgzip -c ',i,'_INDEL_filtered2.vcf > ',i,'_INDEL_filtered2.vcf.gz ; tabix -p vcf ',i,'_INDEL_filtered2.vcf.gz', sep='')
  cat('\n\n')
}


cat('bcftools merge -m all -O v -0 ', sep='')
for (i in samples){
  cat(i, '_SNP_filtered2.vcf.gz', sep='')
  cat(' ')
}
cat('bcftools merge -m all -O v -0 ', sep='')
for (i in samples){
  cat(i, '_INDEL_filtered2.vcf.gz', sep='')
  cat(' ')
}


samples <- c('CT_ALL')

for (i in samples){
  cat('cd /data/NCATS_ifx/data/WES/annotation && convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 ../vcf/',i,'_SNP_filtered4.vcf > ',i,'_SNP4.avinput && table_annovar.pl ',i,'_SNP4.avinput /fdb/annovar/current/hg38 -buildver hg38 -out ',i,'_SNP4_anno -remove -protocol refGene,avsnp150,1000g2015aug_all,exac03,ljb26_all,cosmic70 -operation g,f,f,f,f,f -nastring . -csvout -otherinfo', sep='')
  cat('\n\n')
  cat('cd /data/NCATS_ifx/data/WES/annotation && convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 ../vcf/',i,'_INDEL_filtered4.vcf > ',i,'_INDEL4.avinput && table_annovar.pl ',i,'_INDEL4.avinput /fdb/annovar/current/hg38 -buildver hg38 -out ',i,'_INDEL4_anno -remove -protocol refGene,avsnp150,1000g2015aug_all,exac03,cosmic70 -operation g,f,f,f,f -nastring . -csvout -otherinfo', sep='')
  cat('\n\n')
}




# ${ID}/${ID} --> ',i,'/',i,'



