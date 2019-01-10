cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT01/CT01_dedup.bam -o vcf/CT01_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp /fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz -L /data/NCATS_ifx/data/WES/hglft_genome_subset.bed

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT02/CT02_dedup.bam -o vcf/CT02_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp /fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz -L /data/NCATS_ifx/data/WES/hglft_genome_subset.bed

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT03/CT03_dedup.bam -o vcf/CT03_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp /fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz -L /data/NCATS_ifx/data/WES/hglft_genome_subset.bed

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT04/CT04_dedup.bam -o vcf/CT04_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp /fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz -L /data/NCATS_ifx/data/WES/hglft_genome_subset.bed

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT05/CT05_dedup.bam -o vcf/CT05_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp /fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz -L /data/NCATS_ifx/data/WES/hglft_genome_subset.bed

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g HaplotypeCaller -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT06/CT06_dedup.bam -o vcf/CT06_raw_variants.vcf --genotyping_mode DISCOVERY -stand_call_conf 30 --dbsnp /fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz -L /data/NCATS_ifx/data/WES/hglft_genome_subset.bed
