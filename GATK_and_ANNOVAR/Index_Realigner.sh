cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar BuildBamIndex INPUT=CT01/CT01_dedup.bam && GATK -m 20g RealignerTargetCreator -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT01/CT01_dedup.bam -o CT01/CT01_realignment_targets.list

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar BuildBamIndex INPUT=CT02/CT02_dedup.bam && GATK -m 20g RealignerTargetCreator -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT02/CT02_dedup.bam -o CT02/CT02_realignment_targets.list

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar BuildBamIndex INPUT=CT03/CT03_dedup.bam && GATK -m 20g RealignerTargetCreator -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT03/CT03_dedup.bam -o CT03/CT03_realignment_targets.list

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar BuildBamIndex INPUT=CT04/CT04_dedup.bam && GATK -m 20g RealignerTargetCreator -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT04/CT04_dedup.bam -o CT04/CT04_realignment_targets.list

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar BuildBamIndex INPUT=CT05/CT05_dedup.bam && GATK -m 20g RealignerTargetCreator -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT05/CT05_dedup.bam -o CT05/CT05_realignment_targets.list

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar BuildBamIndex INPUT=CT06/CT06_dedup.bam && GATK -m 20g RealignerTargetCreator -R /fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.fasta -I CT06/CT06_dedup.bam -o CT06/CT06_realignment_targets.list
