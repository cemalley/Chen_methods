# Single cell RNASeq analysis pipeline with Seurat and gene set enrichment analysis
# Claire Malley 2019
# NCATS NIH

library(Seurat)
library(data.table)
library(biomaRt)

setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS010/')
files <- Sys.glob('*dense*csv')

for (file in files){
  dt <- as.data.frame(fread(file))
  names(dt)[1] <- 'ENSG'
  dt <- as.data.table(dt)
  ensg.genes <- data.table("ENSG" = dt$ENSG)
  genes <- ensg.genes$ENSG
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG %in% ensg.genes$ENSG)
  dt <- merge(dt, ensg.genes, by="ENSG")
  dt <- dt[,ENSG:=NULL]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "_dense_expression_matrix.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS010/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
  
}

files <- Sys.glob('*gene_symbol.csv')

for (file in files){
  dt <- fread(file, header=T)
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  sample_id <- str_split_fixed(file, "_gene_symbol.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS010/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
}


reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- x[external_gene_name != "NA",]
  x <- subset(x, select=c("external_gene_name", unlist(names(x))[2:(length(x)-1)] ))
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  barcodes <- names(x)[2:length(x)]
  barcodes <- paste0(barcodes, paste0(".", samplename))
  names(x)[2:length(x)] <- barcodes
  x <- CreateSeuratObject(raw.data = x, project = "IS010")
  x@meta.data$sample <- samplename
  #x <- NormalizeData(x)
 # x <- ScaleData(x, display.progress = F)
  return(x)
}

#setwd("/Volumes/ncatssctl/NGS_related/Chromium/IS010/")
setwd('/data/NCATS_ifx/iPSC/IS010')
files <- Sys.glob('*gene_symbol*.csv')

# on cluster, the order of files is as below. locally, it is different.
H9_CEPT <- reformat_for_seurat(fread(files[2]), "H9_CEPT")
H9_CE <- reformat_for_seurat(fread(files[1]), "H9_CE")
H9_Nociceptor_D28 <- reformat_for_seurat(fread(files[3]), "H9_Nocicptor_D28")
H9_Y27632 <- reformat_for_seurat(fread(files[4]), "H9_Y27632")
Lonza_CEPT <- reformat_for_seurat(fread(files[6]), "Lonza_CEPT")
Lonza_CE <- reformat_for_seurat(fread(files[5]), "Lonza_CE")
Lonza_Nociceptor_D28 <- reformat_for_seurat(fread(files[7]), "Lonza_Nociceptor_D28")
Lonza_Y27632 <- reformat_for_seurat(fread(files[8]), "Lonza_Y27632")

#
is010 <- MergeSeurat(H9_CE, H9_CEPT, project = "IS010", do.normalize = F)
is010 <- MergeSeurat(is010, H9_Y27632, project="IS010", do.normalize = F)
is010 <- MergeSeurat(is010, Lonza_CE, project="IS010", do.normalize = F)
is010 <- MergeSeurat(is010, Lonza_CEPT, project="IS010", do.normalize = F)
is010 <- MergeSeurat(is010, Lonza_Y27632, project="IS010", do.normalize = F)



#on cluster----
setwd('/data/NCATS_ifx/iPSC/IS010')

is010 <- NormalizeData(is010, normalization.method = "LogNormalize", scale.factor = 10000)
is010 <- ScaleData(is010)
is010 <- FindVariableGenes(object = is010, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is010 <- SetAllIdent(is010, id="sample")
is010 <- RunPCA(object = is010, pcs.print = 1:5, genes.print = 5)
PCAPlot(is010)

is010 <- RunTSNE(object = is010, dims.use = 1:8, do.fast = TRUE)
tsneplotdata <- TSNEPlot(object = is010)

pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS010.PCA.pdf"), width=10, height=10)
print(PCAPlot(is010))
dev.off()
pdf(file=paste0(pdfPath, "/IS010.TSNE.pdf"), width=10, height=10)
print(TSNEPlot(is010))
dev.off()
save(is010, file="IS010.Seurat.RData")

# h9-----

is010.h9 <- MergeSeurat(H9_CE, H9_CEPT, project = "IS010")
is010.h9 <- MergeSeurat(is010.h9, H9_Y27632, project="IS010")
is010.h9 <- NormalizeData(object = is010.h9, normalization.method = "LogNormalize", scale.factor = 10000)
is010.h9 <- ScaleData(is010.h9)

is010.h9 <- FindVariableGenes(object = is010.h9, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is010.h9 <- SetAllIdent(is010.h9, id="sample")
is010.h9 <- RunPCA(object = is010.h9, pcs.print = 1:5, genes.print = 5)
is010.h9 <- RunTSNE(object = is010.h9, dims.use = 1:8, do.fast = TRUE)
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS010.H9.PCA.pdf"), width=10, height=10)
print(PCAPlot(is010.h9))
dev.off()
pdf(file=paste0(pdfPath, "/IS010.H9.TSNE.pdf"), width=10, height=10)
print(TSNEPlot(is010.h9))
dev.off()
save(is010.h9, file="IS010.H9.Seurat.RData")
h9.CE.markers <- FindMarkers(object = is010.h9, ident.1 = "H9_CE", min.pct = 0.25)
h9.CEPT.markers <- FindMarkers(object = is010.h9, ident.1 = "H9_CEPT", min.pct = 0.25)
h9.Y27632.markers <- FindMarkers(object = is010.h9, ident.1 = "H9_Y27632", min.pct = 0.25)

h9.CE.markers$GeneId <- row.names(h9.CE.markers)
h9.CEPT.markers$GeneId <- row.names(h9.CEPT.markers)
h9.Y27632.markers$GeneId <- row.names(h9.Y27632.markers)

markers <- merge(h9.CE.markers, h9.CEPT.markers, all=T, by=c(names(h9.CE.markers)))
markers <- merge(markers, h9.Y27632.markers, all=T, by=c(names(h9.CE.markers)))

markers <- as.data.table(markers)
markers <- markers[order(markers$p_val_adj),]
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS010.H9.heatmap.pdf"), width=10, height=10)
print(DoHeatmap(is010.h9, genes.use=c(markers$GeneId), col.low = "green", col.high="red", remove.key = T, slim.col.label = TRUE))
dev.off()
fwrite(markers, 'H9.merged.markers.csv')
names(h9.CE.markers)[3:4] <- c('pct.H9_CE', 'pct.H9_rest')
fwrite(h9.CE.markers, 'H9.CE.markers.csv')

names(h9.CEPT.markers)[3:4] <- c('pct.H9_CEPT', 'pct.H9_rest')
fwrite(h9.CEPT.markers, 'H9.CEPT.markers.csv')

names(h9.Y27632.markers)[3:4] <- c('pct.H9_Y27632', 'pct.H9_rest')
fwrite(h9.Y27632.markers, 'H9.Y27632.markers.csv')

# plot enrichr results-----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS010/Enrichr/')
results <- fread('Lonza-Y27-input,_SIG_EnrichR_Output.csv')
results <- results[order(-score)]

results.go <- results[grepl('GO', libName),]
results.go[,GOID := tstrsplit(term, "\\(")[2]]
results.go[,GOID := gsub("\\(", '', GOID)]
results.go[,GOID := gsub("\\)", '', GOID)]
results.go
results.go[,term.short := tstrsplit(term, "\\(")[1]]
results.go
results.go[,term.short := trimws(term.short, which="right")]
results.go
results.go.clean <- results.go[,c('GOID','adjPval', 'term.short')]
fwrite(results.go.clean, 'Lonza-Y27.Enrichr.output.csv')

results <- fread('H9-CEPT-REVIGO.csv')
results <- na.omit(results)
results <- results[dispensability < 0.7 & eliminated==0,]
results <- results[1:30,]
ggplot(data=results, aes(y=-`log10 p-value`, x=reorder(description, -`log10 p-value`) )) + geom_bar(stat='identity') + coord_flip() +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=15), plot.title = element_text(size=15))+
  labs(y='-log10 p-value', x='', title='Enriched GO Biological Process terms in H9 CEPT vs. Y27632/CE')

# lonza-----
## left off with saving re-made H9 seurat object#
is010.lonza <- MergeSeurat(Lonza_CE, Lonza_CEPT, project = "IS010")
is010.lonza <- MergeSeurat(is010.lonza, Lonza_Y27632, project="IS010")
is010.lonza <- NormalizeData(object = is010.lonza, normalization.method = "LogNormalize", scale.factor = 10000)
is010.lonza <- ScaleData(is010.lonza)
is010.lonza <- FindVariableGenes(object = is010.lonza, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is010.lonza <- SetAllIdent(is010.lonza, id="sample")
is010.lonza <- RunPCA(object = is010.lonza, pcs.print = 1:5, genes.print = 5)
is010.lonza <- RunTSNE(object = is010.lonza, dims.use = 1:8, do.fast = TRUE)
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS010.Lonza.PCA.pdf"), width=10, height=10)
print(PCAPlot(is010.lonza))
dev.off()
pdf(file=paste0(pdfPath, "/IS010.Lonza.TSNE.pdf"), width=10, height=10)
print(TSNEPlot(is010.lonza))
dev.off()
save(is010.lonza, file="IS010.Lonza.Seurat.RData")
lonza.CE.markers <- FindMarkers(object = is010.lonza, ident.1 = "Lonza_CE", min.pct = 0.25)
lonza.CEPT.markers <- FindMarkers(object = is010.lonza, ident.1 = "Lonza_CEPT", min.pct = 0.25)
lonza.Y27632.markers <- FindMarkers(object = is010.lonza, ident.1 = "Lonza_Y27632", min.pct = 0.25)

lonza.CE.markers$GeneId <- row.names(lonza.CE.markers)
lonza.CEPT.markers$GeneId <- row.names(lonza.CEPT.markers)
lonza.Y27632.markers$GeneId <- row.names(lonza.Y27632.markers)

markers <- merge(lonza.CE.markers, lonza.CEPT.markers, all=T, by=c(names(lonza.CE.markers)))
markers <- merge(markers, lonza.Y27632.markers, all=T, by=c(names(lonza.CE.markers)))

markers <- as.data.table(markers)
markers <- markers[order(markers$p_val_adj),]
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS010.Lonza.heatmap.pdf"), width=10, height=10)
print(DoHeatmap(is010.lonza, genes.use=c(markers$GeneId), col.low = "green", col.high="red", remove.key = T, slim.col.label = TRUE))
dev.off()
fwrite(markers, 'Lonza.merged.markers.csv')
names(lonza.CE.markers)[3:4] <- c('pct.Lonza_CE', 'pct.Lonza_rest')
fwrite(lonza.CE.markers, 'Lonza.CE.markers.csv')

names(lonza.CEPT.markers)[3:4] <- c('pct.Lonza_CEPT', 'pct.Lonza_rest')
fwrite(lonza.CEPT.markers, 'Lonza.CEPT.markers.csv')

names(lonza.Y27632.markers)[3:4] <- c('pct.Lonza_Y27632', 'pct.Lonza_rest')
fwrite(lonza.Y27632.markers, 'Lonza.Y27632.markers.csv')




#setwd("/Volumes/ncatssctl/NGS_related/Chromium/IS010/expression_matrices/gene_symbol")
#load("IS010.Seurat.RData")


# Tao's nociceptor D28 merge with IS009 (replace IS009 D28)-----
# note- did not redo this section.
setwd('/data/NCATS_ifx/iPSC/IS009/gene_symbol')
H9_D0 <- fread('IS009_Noc_H9_D0_dense_expression_matrix_genesymbol.csv')
H9_D4 <- fread('IS009_Noc_H9_D4_dense_expression_matrix_genesymbol.csv')
H9_D12 <- fread('IS009_Noc_H9_D12_dense_expression_matrix_genesymbol.csv')
H9_D28 <- fread('/data/NCATS_ifx/iPSC/IS010/IS010_H9_Nociceptor_D28_gene_symbol.csv')

Lonza_D0 <- fread('IS009_Noc_Lonza_D0_dense_expression_matrix_genesymbol.csv')
Lonza_D4 <- fread('IS009_Noc_Lonza_D4_dense_expression_matrix_genesymbol.csv')
Lonza_D12 <- fread('IS009_Noc_Lonza_D12_dense_expression_matrix_genesymbol.csv')
Lonza_D28 <- fread('/data/NCATS_ifx/iPSC/IS010/IS010_Lonza_Nociceptor_D28_gene_symbol.csv')

H9_D0
H9_D4
H9_D12
H9_D28

Lonza_D0
Lonza_D4
Lonza_D12
Lonza_D28

H9_D0 <- reformat_for_seurat(H9_D0, "H9_D0")
H9_D4 <- reformat_for_seurat(H9_D4, "H9_D4")
H9_D12 <- reformat_for_seurat(H9_D12, "H9_D12")
H9_D28 <- reformat_for_seurat(H9_D28, "H9_D28")

Lonza_D0 <- reformat_for_seurat(Lonza_D0, "Lonza_D0")
Lonza_D4 <- reformat_for_seurat(Lonza_D4, "Lonza_D4")
Lonza_D12 <- reformat_for_seurat(Lonza_D12, "Lonza_D12")
Lonza_D28 <- reformat_for_seurat(Lonza_D28, "Lonza_D28")

is009.10.H9 <- MergeSeurat(H9_D0, H9_D4, project = "IS009.10.H9")
is009.10.H9 <- MergeSeurat(is009.10.H9, H9_D12, project="IS009.10.H9")
is009.10.H9 <- MergeSeurat(is009.10.H9, H9_D28, project="IS009.10.H9")
is009.10.H9 <- NormalizeData(object = is009.10.H9, normalization.method = "LogNormalize", scale.factor = 10000)
is009.10.H9 <- ScaleData(is009.10.H9)
is009.10.H9 <- FindVariableGenes(object = is009.10.H9, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is009.10.H9 <- SetAllIdent(is009.10.H9, id="sample")
is009.10.H9 <- RunPCA(object = is009.10.H9, pcs.print = 1:5, genes.print = 5)
is009.10.H9 <- RunTSNE(object = is009.10.H9, dims.use = 1:8, do.fast = TRUE)
save(is009.10.H9, file='IS009-10-H9-nociceptor-seurat.RData')
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS009-10.H9.PCA.pdf"), width=10, height=10)
print(PCAPlot(is009.10.H9))
dev.off()
pdf(file=paste0(pdfPath, "/IS009-10.H9.TSNE.pdf"), width=10, height=10)
print(TSNEPlot(is009.10.H9))
dev.off()


is009.10.Lonza <- MergeSeurat(Lonza_D0, Lonza_D4, project = "IS009.10.Lonza")
is009.10.Lonza <- MergeSeurat(is009.10.Lonza, Lonza_D12, project="IS009.10.Lonza")
is009.10.Lonza <- MergeSeurat(is009.10.Lonza, Lonza_D28, project="IS009.10.Lonza")
is009.10.Lonza <- NormalizeData(object = is009.10.Lonza, normalization.method = "LogNormalize", scale.factor = 10000)
is009.10.Lonza <- ScaleData(is009.10.Lonza)
is009.10.Lonza <- FindVariableGenes(object = is009.10.Lonza, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is009.10.Lonza <- SetAllIdent(is009.10.Lonza, id="sample")
is009.10.Lonza <- RunPCA(object = is009.10.Lonza, pcs.print = 1:5, genes.print = 5)
is009.10.Lonza <- RunTSNE(object = is009.10.Lonza, dims.use = 1:8, do.fast = TRUE)
save(is009.10.Lonza, file='IS009-10-Lonza-nociceptor-seurat.RData')
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS009-10.Lonza.PCA.pdf"), width=10, height=10)
print(PCAPlot(is009.10.Lonza))
dev.off()
pdf(file=paste0(pdfPath, "/IS009-10.Lonza.TSNE.pdf"), width=10, height=10)
print(TSNEPlot(is009.10.Lonza))
dev.off()

# top genes------
load('IS009-10-Lonza-nociceptor-seurat.RData')
is009.10.Lonza <- SetAllIdent(is009.10.Lonza, id="sample")
day28markers <- FindMarkers(object = is009.10.Lonza, ident.1 = 'Lonza_D28', min.pct = 0.25)
day0markers <- FindMarkers(object = is009.10.Lonza, ident.1 = 'Lonza_D0', min.pct = 0.25)
day4markers <- FindMarkers(object = is009.10.Lonza, ident.1 = 'Lonza_D4', min.pct = 0.25)
day12markers <- FindMarkers(object = is009.10.Lonza, ident.1 = 'Lonza_D12', min.pct = 0.25)
fwrite(day0markers, 'IS009-10.Lonza.D0.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(day4markers, 'IS009-10.Lonza.D4.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(day12markers, 'IS009-10.Lonza.D12.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(day28markers, 'IS009-10.Lonza.D28.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
topgenes <- c(row.names(day0markers)[1:10], row.names(day4markers)[1:10], row.names(day12markers)[1:10], row.names(day28markers)[1:10])
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS009-10.Lonza.heatmap.pdf"), width=10, height=10)
print(DoHeatmap(is009.10.Lonza, genes.use=c(topgenes), col.low = "green", col.high="red", remove.key = T, slim.col.label = TRUE))
dev.off()

load('IS009-10-H9-nociceptor-seurat.RData')
is009.10.H9 <- SetAllIdent(is009.10.H9, id="sample")
day28markers <- FindMarkers(object = is009.10.H9, ident.1 = 'H9_D28', min.pct = 0.25)
day0markers <- FindMarkers(object = is009.10.H9, ident.1 = 'H9_D0', min.pct = 0.25)
day4markers <- FindMarkers(object = is009.10.H9, ident.1 = 'H9_D4', min.pct = 0.25)
day12markers <- FindMarkers(object = is009.10.H9, ident.1 = 'H9_D12', min.pct = 0.25)
fwrite(day0markers, 'IS009-10.H9.D0.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(day4markers, 'IS009-10.H9.D4.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(day12markers, 'IS009-10.H9.D12.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(day28markers, 'IS009-10.H9.D28.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
topgenes <- c(row.names(day0markers)[1:10], row.names(day4markers)[1:10], row.names(day12markers)[1:10], row.names(day28markers)[1:10])
pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS009-10.H9.heatmap.pdf"), width=10, height=10)
print(DoHeatmap(is009.10.H9, genes.use=c(topgenes), col.low = "green", col.high="red", remove.key = T, slim.col.label = TRUE))
dev.off()


# plots for specific gene list----
setwd('/data/NCATS_ifx/iPSC/IS010/')
load('IS009-10-H9-nociceptor-seurat.RData')
load('IS009-10-Lonza-nociceptor-seurat.RData')
genelist <- c("POU5F1", "NANOG", "PAX6", "OTX2", "SOX21", "NEUROG1", "SOX10", "POU4F1", "ISL1", "PRPH", "TUBB3", "SCN9A", "SCN10A", "TRPV1", "P2RX3", "CACNA2D1")

pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS009-10.H9.tSNE.markers.pdf"), width=10, height=10)
print(FeaturePlot(is009.10.H9, c(genelist), cols.use = c("darkgreen","red")))
dev.off()

pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS009-10.Lonza.tSNE.markers.pdf"), width=10, height=10)
print(FeaturePlot(is009.10.Lonza, c(genelist), cols.use = c("darkgreen","red")))
dev.off()

pdf(file=paste0(pdfPath, "/IS009-10.H9.VlnPlot.markers.pdf"), width=20, height=20)
VlnPlot(is009.10.H9, c(genelist), remove.legend = T)
dev.off()

pdf(file=paste0(pdfPath, "/IS009-10.Lonza.VlnPlot.markers.pdf"), width=20, height=20)
VlnPlot(is009.10.Lonza, c(genelist), remove.legend = T)
dev.off()


# october 22 analysis of "PT"-----
setwd('/data/NCATS_ifx/iPSC/IS010/')
load('IS010.H9.Seurat.RData')
load('IS010.Lonza.Seurat.RData')

h9.CE.v.CEPT.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25)
h9.CEPT.v.Y27.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CEPT', ident.2='H9_Y27632', min.pct = 0.25)
h9.CE.v.Y27.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_Y27632', min.pct = 0.25)


lonza.CE.v.CEPT.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_CEPT', min.pct = 0.25)
lonza.CEPT.v.Y27.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CEPT', ident.2='Lonza_Y27632', min.pct = 0.25)
lonza.CE.v.Y27.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_Y27632', min.pct = 0.25)

fwrite(h9.CE.v.CEPT.markers, 'IS010.H9.CE.v.CEPT.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(h9.CEPT.v.Y27.markers, 'IS010.H9.CEPT.v.Y27.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(h9.CE.v.Y27.markers , 'IS010.H9.CE.v.Y27.markers .csv', row.names=T, col.names=T, quote=F, sep=',')

fwrite(lonza.CE.v.CEPT.markers, 'IS010.Lonza.CE.v.CEPT.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(lonza.CEPT.v.Y27.markers, 'IS010.Lonza.CEPT.v.Y27.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(lonza.CE.v.Y27.markers , 'IS010.Lonza.CE.v.Y27.markers.csv', row.names=T, col.names=T, quote=F, sep=',')

geneset <- c("FTH1", "GCLC", "GCLM", "GSR", "GSTP1", "HMOX1", "NQO1", "PRDX1", "SQSTM1", "TXN", "TXNRD1", "CASP1", "FAS", "MCL1", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "CDKN1A", "CHEK1", "CHEK2", "DDIT3", "HUS1", "MRE11", "NBN", "RAD17", "RAD9A", "ATM", "ATR", "DDB2", "GADD45A", "GADD45G", "RAD51", "TP53", "XPC", "ATF4", "ATF6", "ATF6B", "BBC3", "BID", "CALR", "DDIT3", "DNAJC3", "HSP90AA1", "HSP90B1", "HSPA4", "HSPA5")

h9.CE.v.CEPT.stress.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(geneset))
h9.CEPT.v.Y27.stress.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CEPT', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(geneset))
h9.CE.v.Y27.stress.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(geneset))

lonza.CE.v.CEPT.stress.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_CEPT', min.pct = 0.25, genes.use = c(geneset)) # no genes pass logfoldchange threshold
lonza.CEPT.v.Y27.stress.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CEPT', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(geneset))
lonza.CE.v.Y27.stress.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(geneset))

cellcycle <- unique(c(unlist(fread('/Volumes/ncatssctl/NGS_related/marker_sets/G2M.genes.txt', header=F), use.names=F), unlist(fread('/Volumes/ncatssctl/NGS_related/marker_sets/S.genes.txt', header=F), use.names=F)))
cat(cellcycle, sep='\n')


h9.CE.v.CEPT.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(cellcycle[c(1:13,15:245)])) # took out: C2ORF69

h9.CE.v.CEPT.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(cellcycle[c(50:58)])) # HN1 --> JPT1
cellcycle[56] <- 'JPT1'

h9.CE.v.CEPT.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(cellcycle[c(58)])) # HRSP12 --> RIDA
h9.CE.v.CEPT.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(cellcycle[c(64)])) # KIAA1524--> CIP2A
cellcycle[145] # C5ORF42 --> CPLANE1
cellcycle[146] # C11ORF82 --> DDIAS 
cellcycle[179] # FAM178A --> SLF2
cellcycle[194] # KIAA1598 --> SHTN1
cellcycle[202] # MLF1IP --> CENPU
cellcycle[207] # NRD1 --> NRDC

h9.CE.v.CEPT.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(cellcycle[c(1:12,15:245)])) 
h9.CEPT.v.Y27.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CEPT', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(cellcycle[c(1:12,15:245)])) # no genes pass logfoldchange threshold
h9.CE.v.Y27.cc.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(cellcycle[c(1:12,15:245)])) 

autophagy <- c("ATG7", "GABARAP", "GABARAPL2", "ATG4A", "ATG4B", "GABARAPL1", "ULK3", "PIK3R4", "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21", "IFNG", "INS", "BECN1P1", "PIK3C3", "PRKAA1", "PRKAA2", "ATG3", "ULK1", "ATG4C", "ATG4D", "BECN1", "ATG12", "ATG5", "ULK2")

h9.CE.v.CEPT.autph.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(autophagy)) 
autophagy[24]<- 'BECN2' #BECN1P1 -->BECN2

# none of the autophagy markers pass logfoldchange threshold.
h9.CEPT.v.Y27.autph.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CEPT', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(autophagy)) #None
h9.CE.v.Y27.autph.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(autophagy))  #none

dnadamagerepair <- c("ABL1", "PARP1", "PARP4", "ADPRTL2", "ADPRTL3", "ALKBH1", "APEX1", "APEX2", "APTX", "ATM", "ATR", "ATRX", "BLM", "BRCA1", "BRCA2", "BRIP1", "CCNH", "CDK7", "CDKN1A", "CETN2", "CETN3", "CHAF1A", "CHAF1B", "CHEK1", "CHEK2", "CIB1", "ERCC8", "CRY1", "CRY2", "CSNK1D", "CSNK1E", "DCLRE1A", "DCLRE1B", "DCLRE1C", "DDB1", "DDB2", "ALKBH3", "DMC1", "DUT", "EME1", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "EXO1", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCL", "FANCM", "FEN1", "ENDOV", "MTOR", "XRCC6", "GADD45A", "GADD45G", "GTF2H1", "GTF2H2", "GTF2H3", "GTF2H4", "GTF2H5", "H2AFX", "HELQ", "HMGB1", "HMGB2", "HUS1", "XRCC6BP1", "LIG1", "LIG3", "LIG4", "MAD1L1", "MAD2L1", "MAD2L2", "MBD1", "MBD2", "MBD3", "MBD4", "MBD5", "ALKBH2", "MGMT", "MLH1", "MLH3", "MMS19", "MNAT1", "MPG", "MRE11A", "MRE11B", "MSH2", "MSH3", "MSH4", "MSH5", "MSH6", "MUS81", "MUTYH", "N4BP2", "NLRP2", "NEIL1", "NEIL2", "NEIL3", "NTHL1", "NUDT1", "NUDT3", "OGG1", "PCNA", "PMS1", "PMS2", "PMS2P1", "PMS2L2", "PMS2P3", "PMS2P4", "PMS2P5", "PMS2P1", "PNKP", "POLA1", "POLA2", "POLB", "POLD1", "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POLE3", "POLE4", "POLG", "POLG2", "POLH", "POLI", "POLK", "POLL", "POLM", "POLN", "POLQ", "PAPD7", "PRKDC", "PTTG1", "RAD1", "RAD17", "RAD18", "RAD21", "RAD23A", "RAD23B", "RAD50", "RAD51", "RAD51AP1", "RAD51C", "RAD51B", "RAD51D", "RAD52", "RDM1", "RAD54B", "RAD54L", "RAD9A", "RAD9B", "RBBP4", "RBBP8", "RECQL", "RECQL4", "RECQL5", "REV1", "REV3L", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5", "RIF1", "RPA1", "RPA2", "RPA3", "RPA4", "RRM1", "RRM2", "RRM2B", "RUVBL2", "SHFM1", "SMC1A", "SMC1B", "SMC2", "SMC3", "SMC4", "SMUG1", "SPO11", "SSRP1", "SUMO1", "SUPT16H", "TDG", "TDP1", "TEP1", "TERF1", "TERF2", "TERT", "TINF2", "TOP1", "TOP2A", "TOP2B", "TOP3A", "TOP3B", "TOPBP1", "TP53", "TP53BP1", "TP53I3", "TP73", "TP63", "TREX1", "TREX2", "UBA1", "UBE2A", "UBE2B", "UBE2D2", "UBE2D3", "UBE2I", "UBE2L3", "UBE2N", "UBE2V2", "UNG", "WDR33", "WRN", "WRNIP1", "XAB2", "XPA", "XPC", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "XRCC5") 


h9.CE.v.CEPT.dna.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_CEPT', min.pct = 0.25, genes.use = c(dnadamagerepair)) 
dnadamagerepair[4] <- 'PARP2' # ADPRTL2
dnadamagerepair[5] <- 'PARP3'
dnadamagerepair[73] <- 'ATP23' # XRCC6BP1
dnadamagerepair[92] <- 'MRE11' #..A
dnadamagerepair[93] <- 'MRE11' #..B
dnadamagerepair[113] <- 'PMS2' # PMS2P1
dnadamagerepair[114] <- 'PMS2' # PMS2L2
dnadamagerepair[115] <- 'PMS2' # PMS2P2
dnadamagerepair[113:118] <- c(rep('PMS2', 6))
dnadamagerepair[140] <- c('TENT4A') #PAPD7
dnadamagerepair[182] <- 'SEM1' # SHFM1

h9.CEPT.v.Y27.dna.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CEPT', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(dnadamagerepair))
h9.CE.v.Y27.dna.markers <- FindMarkers(object = is010.h9, ident.1 = 'H9_CE', ident.2='H9_Y27632', min.pct = 0.25, genes.use = c(dnadamagerepair)) # no genes pass logfoldchange threshold.

lonza.CE.v.CEPT.cc.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_CEPT', min.pct = 0.25, genes.use = c(cellcycle[c(1:12,15:245)])) 
lonza.CEPT.v.Y27.cc.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CEPT', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(cellcycle[c(1:12,15:245)]))
lonza.CE.v.Y27.cc.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(cellcycle[c(1:12,15:245)])) 

lonza.CE.v.CEPT.autph.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_CEPT', min.pct = 0.25, genes.use = c(autophagy)) 
lonza.CEPT.v.Y27.autph.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CEPT', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(autophagy))
lonza.CE.v.Y27.autph.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(autophagy))

lonza.CE.v.CEPT.dna.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_CEPT', min.pct = 0.25, genes.use = c(dnadamagerepair)) 
lonza.CEPT.v.Y27.dna.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CEPT', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(dnadamagerepair))
lonza.CE.v.Y27.dna.markers <- FindMarkers(object = is010.lonza, ident.1 = 'Lonza_CE', ident.2='Lonza_Y27632', min.pct = 0.25, genes.use = c(dnadamagerepair))

# only lonza stress, cc passed logfoldchange thresholds

fwrite(h9.CE.v.CEPT.stress.markers, 'H9.CE.v.CEPT.stress.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(lonza.CE.v.CEPT.stress.markers, 'Lonza.CE.v.CEPT.stress.markers.csv', row.names=T, col.names=T, quote=F, sep=',') # doesnt exist
fwrite(h9.CE.v.Y27.stress.markers, 'H9.CE.v.Y27.stress.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(lonza.CE.v.Y27.stress.markers, 'Lonza.CE.v.Y27.stress.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(h9.CEPT.v.Y27.stress.markers, 'H9.CEPT.v.Y27.stress.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(lonza.CEPT.v.Y27.stress.markers, 'Lonza.CEPT.v.Y27.stress.markers.csv', row.names=T, col.names=T, quote=F, sep=',')

fwrite(h9.CE.v.CEPT.cc.markers, 'H9.CE.v.CEPT.cc.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(lonza.CE.v.CEPT.cc.markers, 'Lonza.CE.v.CEPT.cc.markers.csv', row.names=T, col.names=T, quote=F, sep=',') #doesnt exist
fwrite(h9.CE.v.Y27.cc.markers, 'H9.CE.v.Y27.cc.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(lonza.CE.v.Y27.cc.markers, 'Lonza.CE.v.Y27.cc.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(h9.CEPT.v.Y27.cc.markers, 'H9.CEPT.v.Y27.cc.markers.csv', row.names=T, col.names=T, quote=F, sep=',') #doesnt exist
fwrite(lonza.CEPT.v.Y27.cc.markers, 'Lonza.CEPT.v.Y27.cc.markers.csv', row.names=T, col.names=T, quote=F, sep=',')

# none of the autophagy genes were significant.

fwrite(h9.CE.v.CEPT.dna.markers, 'H9.CE.v.CEPT.dna.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(lonza.CE.v.CEPT.dna.markers, 'Lonza.CE.v.CEPT.dna.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(h9.CE.v.Y27.dna.markers, 'H9.CE.v.Y27.dna.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(lonza.CE.v.Y27.dna.markers, 'Lonza.CE.v.Y27.dna.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
fwrite(h9.CEPT.v.Y27.dna.markers, 'H9.CEPT.v.Y27.dna.markers.csv', row.names=T, col.names=T, quote=F, sep=',')
#fwrite(lonza.CEPT.v.Y27.dna.markers, 'Lonza.CEPT.v.Y27.dna.markers.csv', row.names=T, col.names=T, quote=F, sep=',')

# locally, summarize pathway DE tests-----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS010/DE/pathways_DE/')
ccfiles <- Sys.glob('*cc.markers.csv')
dnafiles <- Sys.glob('*dna.markers.csv')
stressfiles <- Sys.glob('*stress.markers.csv')

data <- data.table()

for (i in stressfiles){
    
  if (nrow(data)>0){
    data.temp <- fread(i, header=T, select=c(1,3,6))
    sample1 <- str_split_fixed(i, '\\.', Inf)[3]
    sample2 <- str_split_fixed(i, '\\.', Inf)[5]
    cell_line <- str_split_fixed(i, '\\.', Inf)[2]
    
    names(data.temp) <- c('GeneId', paste0(cell_line, '_', sample1, '_', sample2, '_LFC'), paste0(cell_line, '_', sample1, '_', sample2, '_padj'))
    
    data <- merge(data, data.temp, by=c('GeneId'), all=T)
  }
  
  if (nrow(data)==0){
    data <- fread(i, header=T, select=c(1,3,6))
    sample1 <- str_split_fixed(i, '\\.', Inf)[3]
    sample2 <- str_split_fixed(i, '\\.', Inf)[5]
    cell_line <- str_split_fixed(i, '\\.', Inf)[2]
    
    names(data) <- c('GeneId', paste0(cell_line, '_', sample1, '_', sample2, '_LFC'), paste0(cell_line, '_', sample1, '_', sample2, '_padj'))
  }
}

ccdata <- data
dnadata <- data
stressdata <- data

ccdata[,pathway:= 'Cell cycle']
dnadata[,pathway := 'DNA damage and repair']
stressdata[,pathway:= 'Stress panel']

data <- data.table()

data <- merge(ccdata, dnadata, all=T, by=c(intersect(names(ccdata), names(dnadata))))

data <- merge(data, stressdata, all=T, by=c(intersect(names(data), names(stressdata))))

data <- data[order(pathway)]
data <- data[,c("GeneId","pathway", "H9_CE_CEPT_LFC", "H9_CE_CEPT_padj", "H9_CE_Y27_LFC", "H9_CE_Y27_padj", "Lonza_CE_Y27_LFC", "Lonza_CE_Y27_padj", "Lonza_CEPT_Y27_LFC", "Lonza_CEPT_Y27_padj", "H9_CEPT_Y27_LFC", "H9_CEPT_Y27_padj")]

fwrite(data, 'CC.DNA.Stress.DE.csv', sep=',', col.names = T, row.names = F, quote=F)

# make heatmap of all the pathway DE genes -----
load('IS010.H9.Seurat.RData')
load('IS010.Lonza.Seurat.RData')

markers <-c("KPNA2", "MALAT1", "NEAT1", "TUBA1A", "TUBB2A", "DUT", "XRCC5", "FTH1", "GSTP1", "HSP90AA1", "HSP90B1", "HSPA5", "NQO1", "PRDX1", "SQSTM1")

pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS010.H9.pathway.heatmap.pdf"), width=10, height=10)
DoHeatmap(is010.h9, genes.use=c(markers), col.low = "green", col.high="red", remove.key = F, slim.col.label = TRUE)
VlnPlot(is010.h9, c(markers))

dev.off()


# locally----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS010/')
load('IS010.H9.Seurat.RData')

AverageExpression(is010.h9, 'FTH1')
FindMarkers(is010.h9, ident.1='H9_CE', ident.2 = 'H9_Y27632', c('FTH1'))

AverageExpression(is010.h9, c(markers))

FindMarkers(is010.h9, ident.1='H9_CEPT', ident.2 = 'H9_Y27632', c('GSTP1'))

stress <- c("FTH1", "GSTP1", "HSP90AA1", "HSP90B1", "HSPA5", "NQO1", "PRDX1", "SQSTM1")

FindMarkers(is010.h9, ident.1='H9_CEPT', ident.2 = 'H9_Y27632', c(stress))


is010.h9 <- CellCycleScoring(is010.h9, s.genes = s.genes, g2m.genes = g2.m.genes, set.ident = TRUE)
#is010.h9 <- ScaleData(object = is010.h9, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = FALSE)


#is010.h9 <- RunTSNE(is010.h9)
#TSNEPlot(object = is010.h9)
is010.h9 <- SetAllIdent(is010.h9, id="Phase")
PCAPlot(is010.h9)
is010.h9 <- RunTSNE(is010.h9)
TSNEPlot(object = is010.h9)
is010.h9 <- SetAllIdent(is010.h9, id="sample")
TSNEPlot(object = is010.h9)

# no real difference by phase.

pdata <- as.data.frame(is010.h9@meta.data)



ggplot(data=pdata, aes(color=sample)) + geom_violin(aes(x=sample, y=G2M.Score, group=sample)) + geom_jitter(aes(x=sample, y=G2M.Score), width=0.1)+ 
  labs(title='G2M cell cycle scores for IS010 H9', x='Sample')+scale_color_discrete(name = "Sample")

ggplot(data=pdata, aes(color=sample)) + geom_violin(aes(x=sample, y=S.Score, group=sample)) + geom_jitter(aes(x=sample, y=S.Score), width=0.1)+ 
  labs(title='S cell cycle scores for IS010 H9', x='Sample')+scale_color_discrete(name = "Sample")

is010.h9 <- FindClusters(is010.h9)
is010.h9@meta.data

is010.h9 <- SetAllIdent(is010.h9, id="res.0.8")
TSNEPlot(is010.h9)

cluster3.v.rest <- FindMarkers(is010.h9, ident.1='3')
fwrite(cluster3.v.rest, 'DE/Cluster3.v.rest.csv', row.names = T, col.names = T, quote=F, sep=',')

AverageExpression(is010.h9, 'DDIT4')

#
