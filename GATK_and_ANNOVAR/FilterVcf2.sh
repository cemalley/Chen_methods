cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT01_SNP_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 60.0) & (MQ < 40.0)" -o vcf/CT01_SNP_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT01_INDEL_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 200.0) & (SOR > 10.0)" -o vcf/CT01_INDEL_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT02_SNP_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 60.0) & (MQ < 40.0)" -o vcf/CT02_SNP_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT02_INDEL_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 200.0) & (SOR > 10.0)" -o vcf/CT02_INDEL_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT03_SNP_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 60.0) & (MQ < 40.0)" -o vcf/CT03_SNP_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT03_INDEL_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 200.0) & (SOR > 10.0)" -o vcf/CT03_INDEL_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT04_SNP_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 60.0) & (MQ < 40.0)" -o vcf/CT04_SNP_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT04_INDEL_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 200.0) & (SOR > 10.0)" -o vcf/CT04_INDEL_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT05_SNP_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 60.0) & (MQ < 40.0)" -o vcf/CT05_SNP_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT05_INDEL_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 200.0) & (SOR > 10.0)" -o vcf/CT05_INDEL_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT06_SNP_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 60.0) & (MQ < 40.0)" -o vcf/CT06_SNP_filtered2.vcf

cd /data/NCATS_ifx/data/WES/ && GATK -m 15g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V vcf/CT06_INDEL_raw.vcf -filterName "my_filter" --filterExpression "(QD < 2.0) & (FS > 200.0) & (SOR > 10.0)" -o vcf/CT06_INDEL_filtered2.vcf 
