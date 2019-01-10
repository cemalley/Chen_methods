# Subset bedfile from WES targets list to only main chromosomes
# Claire Malley
# NCATS NIH

bed <- fread('/Volumes/ncatssctl/NGS_related/WES/hglft_genome_1845_3870.bed', sep="\t", header=F) #367779 rows
bed <- subset(bed, V1 %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')) #367520 rows
fwrite(bed, '/Volumes/ncatssctl/NGS_related/WES/hglft_genome_subset.bed', sep='\t', col.names=F)
