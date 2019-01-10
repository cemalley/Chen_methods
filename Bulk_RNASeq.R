# Bulk RNASeq analysis pipeline: FASTQ combining from different lanes, TRIMMOMATIC, STAR, HTSeq-count, DESeq2, pathway enrichment analysis with Enrichr (see enrichr.R) and REVIGO, linear regression models, and various plots of specific genes
# Claire Malley 2018
# NCATS NIH

library(data.table)
library(DESeq2)
library(stringr)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/')
files <- fread('./File_organization/file-paths.txt', header=F)
files[,(c('V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11')):= tstrsplit(V1,'/')]
files[,Folder:= paste('',V3, V4, V5, V6, V7, V8, V9, V10, sep='/')]


for (i in 1:nrow(files)){
  if (grepl('ISB008', files$V11[i])){
    sampleID <- str_split_fixed(files$V10[i], '_L', Inf)[1]
    files$Sample[i] <- sampleID
  }
  if (grepl('ISB008', files$V10[i])){
    sampleID <- str_split_fixed(files$V11[i], '_S', Inf)[1]
    sampleID <- gsub('-', '_', sampleID)
    files$Sample[i] <- sampleID
  }
}

files[,FullPath := Folder]
files[,Folder := V10]
files[,File := V11]
files <- files[,c('FullPath', 'Folder', 'File', 'Sample')]

grep('FASTQ_Generation_2018-10-11_16_48_59Z-129381954', files$FullPath[i])

files <- files[-c(1:32),]


setwd('./File_organization/')
fwrite(files, 'files.csv')

samples <- unique(files$Sample)

for (sample in samples){
  rows <- grep(paste0("^",sample,"$"), files$Sample)
  cat(sample)
  cat('\t')
  cat(rows)
  cat('\n')
}

rows <- fread('rows.csv')

# fastq cat -----
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp[c(1,3,5,7)]){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0(files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], ' '))
  }
  cat(paste0('> ', sample, '_R1_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
  
}

for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp[c(2,4,6,8)]){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0(files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], ' '))
  }
  cat(paste0('> ', sample, '_R2_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
  
}
#swarm -f fastq_cat_102318.sh -t 1 -g 1 --time=8:00:00 #12051249

# trimmomatic ----
dir <- '/data/NCATS_ifx/data/mRNASeq/ISB008'

for (sample in samples){
  cat(paste0("java -jar $TRIMMOJAR PE -phred33 ", dir ,'/', sample, '_R1_001.fastq.gz ', dir, '/',
             sample, '_R2_001.fastq.gz ', '-baseout ', sample, '.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36'), sep="\n")
  cat("\n")
}

# swarm -f trim_swarm_102318.sh -g 12 -t 8 --module trimmomatic --time=8:00:00  #12051880


# STAR -----
for (sample in samples){
  cat(paste0('cd ', dir, ' && mkdir -p bam/' , 'Sample_' , sample , ' && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within –-outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 -- sjdbScore 1 --readFilesIn ', dir, '/' , sample , '_1P.fastq.gz ', dir, '/', sample , '_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/Sample_' , sample ,'/' , sample, '_hg38'))
  cat("\n\n")
}

# swarm -f star_swarm_102318.sh -t 12 -g 40 --module STAR --time=8:00:00 #12054308

# htseq-----

indir <- '/data/NCATS_ifx/data/mRNASeq/ISB008'
outdir <- '/data/NCATS_ifx/data/mRNASeq/ISB008/htseq'
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/')
#samples <- unlist(fread('samples.txt'), use.names=F)

for (sample in samples){
  cat(paste0('htseq-count -f bam -r pos -s no -t exon -m union ', indir , '/bam/Sample_' , sample  , '/' , sample , '_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', outdir, '/', sample , '_htseq_counts.txt' , "\n\n"))
}
#swarm -f htseq_swarm_102318.sh -t 4 -g 4 --module htseq --time=10:00:00 #12061435

# start DESeq2 work-----
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("biomaRt")
library(biomaRt)
library(stringr)
library(dendextend)
library(ComplexHeatmap)
library(gtools)

# convert ENSG.# genes to ENSG, then gene symbols----

setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/htseq/Countfiles_ENSG/")
files <- Sys.glob("*counts.txt")

for (file in files){
  
  dt <- fread(file)
  dt <- dt[1:(nrow(dt)-5),]
  names(dt) <- c("ENSG.full", "Counts")
  
  ensg.genes <- data.table("ENSG.full" = dt$ENSG.full)
  
  ensg.genes[,ENSG.short := tstrsplit(ENSG.full, "\\.")[1]]
  
  genes <- ensg.genes$ENSG.short
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG.short", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG.full %in% ensg.genes$ENSG.full)
  dt <- merge(dt, ensg.genes, by="ENSG.full")
  dt <- dt[,c("external_gene_name", "Counts")]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "_htseq_counts.txt",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/htseq/Countfiles_gene_symbols/", sample_id, "_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
  
}
#start DESeq2 work----

setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008")
sampleFiles <- gtools::mixedsort(grep("counts.txt",list.files("./htseq/Countfiles_gene_symbols/"),value=TRUE))
sampleTable <- data.table('sampleFiles' = sampleFiles, 'fileName' = sampleFiles)
sampleTable[,condition:= tstrsplit(sampleFiles, '_Sample')[1]]
sampleTable
fwrite(sampleTable, 'sampletable.csv', sep=',', quote=F, col.names = T, row.names = F)

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/htseq/Countfiles_gene_symbols/"
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

save(dds, file='ISB008.DDS.RData')
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/')
load('ISB008.DDS.RData')

data <- as.data.frame(counts(dds, normalized=TRUE))
newnames <- row.names(dds@colData)

newnames <- vapply(newnames, function(x) str_split_fixed(x, '_gene_symbol_counts.txt', Inf)[1], character(1))


names(data) <- newnames
data
fwrite(data, 'ISB008.normalized.counts.csv', sep=',', quote=F, col.names = T, row.names = T)

data <- as.data.frame(counts(dds, normalized=FALSE))
names(data) <- newnames
data
fwrite(data, 'ISB008.raw.counts.csv', sep=',', quote=F, col.names = T, row.names = T)

cat(vapply(newnames, function(x) str_split_fixed(x, '_Sample', Inf)[1], character(1)), sep='\n')

# october 29, renaming a few samples -----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/')
load('ISB008.DDS-old-before-renaming.RData')
# A02 Lonza_CE_72hr_Sample_3 --> Lonza_CEPT_24hr_Sample_3
# B01 Lonza_CEPT_24hr_Sample_3 --> Lonza_CE_72hr_Sample_3
# H10 Lonza_CEPT_24hr_Sample_2a --> Lonza_CEPT_72hr_Sample_2

setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008")
sampleFiles <- gtools::mixedsort(grep("counts.txt",list.files("./htseq/Countfiles_gene_symbols/"),value=TRUE))
sampleTable <- data.table('sampleFiles' = sampleFiles, 'fileName' = sampleFiles)
sampleTable[,condition:= tstrsplit(sampleFiles, '_Sample')[1]]
sampleTable
fwrite(sampleTable, 'sampletable.csv', sep=',', quote=F, col.names = T, row.names = F)

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/htseq/Countfiles_gene_symbols/"
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

save(dds, file='ISB008.DDS.RData')
#setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/')
#load('ISB008.DDS.RData')

data <- as.data.frame(counts(dds, normalized=TRUE))
newnames <- row.names(dds@colData)

newnames <- vapply(newnames, function(x) str_split_fixed(x, '_gene_symbol_counts.txt', Inf)[1], character(1))

names(data) <- newnames
data
fwrite(data, 'ISB008.normalized.counts.csv', sep=',', quote=F, col.names = T, row.names = T)

data <- as.data.frame(counts(dds, normalized=FALSE))
names(data) <- newnames
data
fwrite(data, 'ISB008.raw.counts.csv', sep=',', quote=F, col.names = T, row.names = T)


# planning all DE tests----

celllines <- c('H9', 'Lonza')

times <- c('24hr', '72hr', 'Thaw')

conditions <- c('CE', 'CEPT', 'Y27', 'EDTA', 'CH')
expand.grid(celllines,conditions, times)


#i.e. H9_CE_24hr vs H9_CEPT_24hr


text <- capture.output( for (i in times){
  for (b in celllines){
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[2], '_', i, '\n'))
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[3], '_', i, '\n'))
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[4], '_', i, '\n'))
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  cat(paste0(b, '_', conditions[2], '_', i, ' vs ', b, '_', conditions[3], '_', i, '\n'))
  cat(paste0(b, '_', conditions[2], '_', i, ' vs ', b, '_', conditions[4], '_', i, '\n'))
  cat(paste0(b, '_', conditions[2], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  cat(paste0(b, '_', conditions[3], '_', i, ' vs ', b, '_', conditions[4], '_', i, '\n'))
  cat(paste0(b, '_', conditions[3], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  cat(paste0(b, '_', conditions[4], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  }
})



tests_table <- data.table('test'=text)
tests_table[,c('condition1', 'condition2'):=tstrsplit(test, ' vs ')]
tests_table <- tests_table[!grepl('EDTA_Thaw',tests_table$test),] # there is no EDTA_Thaw condition


# profiling doParallel----
# 
# library(parallel)
# library(doParallel)
# 
# ptime.parallel <- system.time({cl <- makeCluster(6)
# registerDoParallel(cl)
# foreach(i=1:60, .packages = "data.table")  %dopar% {print(tests_table[i])}
# stopCluster(cl)})[3]
# 
# 
# ptime.serial <- system.time({cl <- makeCluster(6)
# registerDoParallel(cl)
# foreach(i=1:60, .packages = "data.table")  %do% {print(tests_table[i])}
# stopCluster(cl)})[3]
# 

# parallelize deseq2-----

testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE/Volcano'

#syntax
#DESeq2_pipeline <- function(dds, condition1, condition2, testdir, volcanodir)

#DESeq2_pipeline(dds, 'H9_CE_24hr', 'H9_CEPT_24hr', testdir, volcanodir)

cl <- makeCluster(8)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]

  testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE'
  volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE/Volcano'

  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)
}

stopCluster(cl)


# find pathway commonalities across DE tests----

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE')
data <- fread('DE.H9_CEPT_Thaw.vs.H9_Y27_Thaw.csv')
fwrite(data[1:200,c('GeneId')], './Enrichr/H9_CEPT_Thaw.vs.H9_Y27_Thaw.enrichr.csv', col.names = F, row.names = F, quote=F, sep='\t')


# REVIGO ----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE/Enrichr/')
results <- fread('H9_CEPT_Thaw.vs.H9_Y27_Thaw.enrichr.csv,_SIG_EnrichR_Output.csv')
results <- results[order(-score)]

results.go <- results[grepl('GO', libName),]
results.go[,GOID := tstrsplit(term, "\\(GO")[2]]
results.go[,GOID := paste0('GO', GOID)]
results.go[,GOID := gsub("\\)", '', GOID)]
results.go
results.go[,term.short := tstrsplit(term, "\\(")[1]]
results.go
results.go[,term.short := trimws(term.short, which="right")]
results.go
results.go.clean <- results.go[,c('GOID','adjPval', 'term.short')]


results.go.clean
fwrite(results.go.clean, 'H9_CEPT_Thaw.vs.H9_Y27_Thaw.enrichr.csv,_SIG_EnrichR_Output_clean.csv')

library(ggplot2)
#ggplot(data=results.go.clean, aes(y=-log10(adjPval), x=reorder(term.short, -adjPval) )) + geom_bar(stat='identity') + coord_flip() # too many terms

results <- fread('H9_CEPT_Thaw.vs.H9_Y27_Thaw.enrichr.csv,_SIG_EnrichR_Output_REVIGO.csv')
results <- na.omit(results)
results <- results[dispensability < 0.7,]
ggplot(data=results, aes(y=-`log10 p-value`, x=reorder(description, -`log10 p-value`) )) + geom_bar(stat='identity') + coord_flip() +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=15), plot.title = element_text(size=20))+
  labs(y='-log10 p-value', x='', title='Enriched GO Biological Process terms in H9 CEPT Thaw vs. H9 Y27 Thaw')
#

# plot counts of particular genes-----
genelist <- c("NANOG", "POU5F1", "ESRG", 'CDKN1A', 'H1F0', 'GDF15', 'BTG2')

length(genelist)

dev.off()
par(mfrow=c(2,4))

for (i in 1:length(genelist)){
  gene <- genelist[i]
  main <- gene
  DESeq2::plotCounts(dds, gene=gene, intgroup="condition", main = main)
}
DESeq2::plotCounts(dds, gene='NANOG', intgroup="condition", main = main)

dds.subset <- as.data.frame(counts(dds, normalized=TRUE))
dds.subset <- subset(dds.subset, row.names(dds.subset) %in% c(genelist))
View(dds.subset)
dds.subset <- subset(dds.subset, select=c('H9_CEPT_Thaw_Sample_1_gene_symbol_counts.txt', 'H9_CEPT_Thaw_Sample_2_gene_symbol_counts.txt', 'H9_CEPT_Thaw_Sample_3_gene_symbol_counts.txt', 'H9_Y27_Thaw_Sample_1_gene_symbol_counts.txt', 'H9_Y27_Thaw_Sample_2_gene_symbol_counts.txt', 'H9_Y27_Thaw_Sample_3_gene_symbol_counts.txt'))
View(dds.subset)

data <- as.data.frame(counts(dds, normalized=TRUE))


# find similarities in top genes between H9 and Lonza, lowering thresholds for consideration----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE')
h9 <- fread('DE.H9_CEPT_Thaw.vs.H9_Y27_Thaw.csv')
lonza <- fread('DE.Lonza_CEPT_Thaw.vs.Lonza_Y27_Thaw.csv')

h9 <- h9[,c(1:4,8)]
lonza <- lonza[,c(1:4,8)]
names(h9)[4:5] <- paste0('H9_', names(h9)[4:5])
names(lonza)[4:5] <- paste0('Lonza_', names(lonza)[4:5])
merged <- merge(h9, lonza, by='GeneId', all=T)

fwrite(merged, 'Merged.H9.Lonza.CEPT.thaw.vs.Y27.csv')

merged[Lonza_padj < 0.001 & H9_padj < 0.001,]
cat(unlist(merged[Lonza_padj < 0.001 & H9_padj < 0.001,c('GeneId')], use.names=F), sep='\n')

edta <- fread('DE.H9_CEPT_24hr.vs.H9_EDTA_24hr.csv')

restoplot <- na.omit(edta)
anno <- subset(restoplot, ((abs(log2FoldChange) > 1) & (padj < 0.001)))

attach(edta)

ggplot(data=edta, aes(x=log(H9_CEPT_24hr), y=log(H9_EDTA_24hr))) + geom_point(color=ifelse((edta$padj < 0.001 & abs(edta$log2FoldChange) > 1), 'red', 'black'))+ labs(title='H9 CEPT 24hr vs. EDTA 24hr expression') +
  geom_text_repel(data=anno, aes(x=log(H9_CEPT_24hr), y=log(H9_EDTA_24hr), label=GeneId),
                  segment.size = 0.5)

linearMod <- lm(log(H9_EDTA_24hr+1) ~ log(H9_CEPT_24hr+1), data=edta)
print(linearMod)
summary(linearMod)

# y = 0.002467 + 0.999146x + 0.2163
# adjusted r-squared = 0.9866
#
fwrite(anno, 'H9_CEPT_24hr_vs_H9_EDTA_24hr_lreg_redpoints.csv')


## 72hr linear regression, CEPT vs EDTA 72hr ----
edta <- fread('./DE/DE.H9_CEPT_72hr.vs.H9_EDTA_72hr.csv')

restoplot <- na.omit(edta)
#version 1: anno <- subset(restoplot, ((abs(log2FoldChange) > 1) & (padj < 0.001)))
#new version:
anno <- subset(restoplot, restoplot$GeneId %in% c('CXCL5', 'DNMT3B', 'HESX1', 'IDO1', 'LCK', 'NANOG', 'POU5F1', 'SOX2', 'TRIM22', 'CDH1', 'CD30', 'SSEA4', 'SSEA3', 'GDF3', 'EPCAM', 'REX1', 'LIN28', 'PODXL', 'SALL4', 'EHMT2', 'APOE', 'CDH3', 'ERBB3', 'LEFTY1', 'GCTM2', 'CD24', 'DUSP6', 'ZFP42', 'ESRG'))
attach(edta)

# version 1 with red points: ggplot(data=edta, aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr))) + geom_point(color=ifelse((edta$padj < 0.001 & abs(edta$log2FoldChange) > 1), 'red', 'black'))+ labs(title='H9 CEPT 72hr vs. EDTA 72hr expression') +  geom_text_repel(data=anno, aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.5)

linearMod <- lm(log(H9_EDTA_72hr+1) ~ log(H9_CEPT_72hr+1), data=edta)
print(linearMod)
summary(linearMod)

anno[,model:=log(anno$H9_EDTA_72hr) - (coef(linearMod)[1] + coef(linearMod)[2]*log(anno$H9_CEPT_72hr))]

#new version: --plot all points--
# ggplot(data=edta, aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr))) + geom_point()+ labs(title='H9 CEPT 72hr vs. EDTA 72hr expression') +  
#   geom_text_repel(data=anno[model > 0,], aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.8, nudge_x = -1.5, size=4)+
#   geom_text_repel(data=anno[model < 0,], aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.8, nudge_x = 1.5, size=4)+
#   theme_light()+
#   theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=1, color='#e0e0e0'), panel.grid.major = element_line(size=1, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+
#   labs(x='log(H9 CEPT 72hr)', y='log(H9 EDTA 72hr)')+
#   annotate('text', x=0.5,y=11.5, label='Adjusted R-squared = 0.99\ny = -0.06 + 1.00x +/- 0.12', size=4)

# new version -- plot only pluripotency genes--
attach(anno)

fun.1 <- function(x) (coef(linearMod)[2])*x + coef(linearMod)[1]

ggplot(data=anno, aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr))) + geom_point()+ labs(title='H9 CEPT 72hr vs. EDTA 72hr expression') +  
  geom_text_repel(data=anno[model > 0,], aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.2, nudge_x = -1.5, size=4)+
  geom_text_repel(data=anno[model < 0,], aes(x=log(H9_CEPT_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.2, nudge_x = 1.5, size=4)+
  theme_light()+
  theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.ticks = element_line(size=1, color='#e0e0e0'))+
  labs(x='log(H9 CEPT 72hr)', y='log(H9 EDTA 72hr)')+
  annotate('text', x=1,y=11,hjust=0, label='Adjusted R-squared = 0.99\ny = 1.00x -0.06 +/- 0.12', size=4)+
  stat_function(fun = fun.1)+
  scale_x_continuous(expand=c(0,.1), limits=c(0,12)) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,12)) +
  coord_cartesian(xlim=c(0,12), ylim=c(0,12)) 




# y = -0.06499 + 1.00933x +/- 0.1156
# adjusted r-squared = 0.9937
#
fwrite(anno, 'H9_CEPT_72hr_vs_H9_EDTA_72hr_lreg_redpoints.csv')

## 72hr linear regression, H9 Y27 vs. EDTA ----
edta <- fread('./DE/DE.H9_Y27_72hr.vs.H9_EDTA_72hr.csv')
restoplot <- na.omit(edta)
anno <- subset(restoplot, restoplot$GeneId %in% c('CXCL5', 'DNMT3B', 'HESX1', 'IDO1', 'LCK', 'NANOG', 'POU5F1', 'SOX2', 'TRIM22', 'CDH1', 'CD30', 'SSEA4', 'SSEA3', 'GDF3', 'EPCAM', 'REX1', 'LIN28', 'PODXL', 'SALL4', 'EHMT2', 'APOE', 'CDH3', 'ERBB3', 'LEFTY1', 'GCTM2', 'CD24', 'DUSP6', 'ZFP42', 'ESRG'))
attach(anno)
linearMod <- lm(log(H9_EDTA_72hr+1) ~ log(H9_Y27_72hr+1), data=edta)
anno[,model:=log(anno$H9_EDTA_72hr) - (coef(linearMod)[1] + coef(linearMod)[2]*log(anno$H9_Y27_72hr))]

# ggplot(data=edta, aes(x=log(H9_Y27_72hr), y=log(H9_EDTA_72hr))) + geom_point()+ labs(title='H9 Y27 72hr vs. EDTA 72hr expression') +  
#   geom_text_repel(data=anno[model > 0,], aes(x=log(H9_Y27_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.8, nudge_x = -1.5, size=4)+
#   geom_text_repel(data=anno[model < 0,], aes(x=log(H9_Y27_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.8, nudge_x = 1.5, size=4)+
#   theme_light()+
#   theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=1, color='#e0e0e0'), panel.grid.major = element_line(size=1, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+
#   labs(x='log(H9 Y27 72hr)', y='log(H9 EDTA 72hr)')+
#   annotate('text', x=0.5,y=11.5, label='Adjusted R-squared = 0.99\ny = 0.01 + 0.99x +/- 0.16', size=4)

ggplot(data=anno, aes(x=log(H9_Y27_72hr), y=log(H9_EDTA_72hr))) + geom_point()+ labs(title='H9 Y27 72hr vs. EDTA 72hr expression') +  
  geom_text_repel(data=anno[model > 0,], aes(x=log(H9_Y27_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.2, nudge_x = -1.5, size=4)+
  geom_text_repel(data=anno[model < 0,], aes(x=log(H9_Y27_72hr), y=log(H9_EDTA_72hr), label=GeneId), segment.size = 0.2, nudge_x = 1.5, size=4)+
  theme_light()+
  theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.ticks = element_line(size=1, color='#e0e0e0'))+
  labs(x='log(H9 Y27 72hr)', y='log(H9 EDTA 72hr)')+
  annotate('text', x=1,y=11,hjust=0, label='Adjusted R-squared = 0.99\ny = 0.99x +0.01 +/- 0.16', size=4)+
  stat_function(fun = fun.1)+
  scale_x_continuous(expand=c(0,.1), limits=c(0,12)) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,12)) +
  coord_cartesian(xlim=c(0,12), ylim=c(0,12))
