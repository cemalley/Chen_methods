# snps-----
GATK -m 9g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V CT_ALL_SNP_filtered2.vcf --genotypeFilterName "lowDP" --genotypeFilterExpression "(DP < 7)" --genotypeFilterName "lowGQ" --genotypeFilterExpression "(GQ < 30)" -o CT_ALL_SNP_filtered3.vcf

GATK -m 9g SelectVariants -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V CT_ALL_SNP_filtered3.vcf --setFilteredGtToNocall -o CT_ALL_SNP_filtered4.vcf


# indels-----
GATK -m 9g VariantFiltration -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V CT_ALL_INDEL_filtered2.vcf --genotypeFilterName "lowDP" --genotypeFilterExpression "(DP < 7)" --genotypeFilterName "lowGQ" --genotypeFilterExpression "(GQ < 30)" -o CT_ALL_INDEL_filtered3.vcf

GATK -m 9g SelectVariants -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -V CT_ALL_INDEL_filtered3.vcf --setFilteredGtToNocall -o CT_ALL_INDEL_filtered4.vcf
