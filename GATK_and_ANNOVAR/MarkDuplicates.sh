cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
  INPUT=CT01/CT01.bam \
  OUTPUT=CT01/CT01_dedup.bam \
  METRICS_FILE=metrics_CT01.txt

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
  INPUT=CT02/CT02.bam \
  OUTPUT=CT02/CT02_dedup.bam \
  METRICS_FILE=metrics_CT02.txt

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
  INPUT=CT03/CT03.bam \
  OUTPUT=CT03/CT03_dedup.bam \
  METRICS_FILE=metrics_CT03.txt

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
  INPUT=CT04/CT04.bam \
  OUTPUT=CT04/CT04_dedup.bam \
  METRICS_FILE=metrics_CT04.txt

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
  INPUT=CT05/CT05.bam \
  OUTPUT=CT05/CT05_dedup.bam \
  METRICS_FILE=metrics_CT05.txt

cd /data/NCATS_ifx/data/WES/ && java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates \
  INPUT=CT06/CT06.bam \
  OUTPUT=CT06/CT06_dedup.bam \
  METRICS_FILE=metrics_CT06.txt
