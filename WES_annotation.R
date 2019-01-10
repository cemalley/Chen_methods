# Analysis of WES SNPs and indels after annotation with ANNOVAR.
# Claire Malley 2018
# NCATS NIH

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(ggcorrplot)
source("http://www.sthda.com/upload/rquery_cormat.r")
library(mvnormtest)
library(ggrepel)
library(factoextra)

setwd('/Volumes/ncatssctl/NGS_related/WES/annotation/')

anno.files <- Sys.glob('CT_ALL*SNP4*multianno*csv')
anno.files

samples <- c('CT_ALL')

snps <- data.table()

snp.header <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "avsnp150", "1000g2015aug_all", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_pred", "FATHMM_score", "FATHMM_pred", "RadialSVM_score", "RadialSVM_pred", "LR_score", "LR_pred", "VEST3_score", "CADD_raw", "CADD_phred", "GERP++_RS", "phyloP46way_placental", "phyloP100way_vertebrate", "SiPhy_29way_logOdds", "cosmic70", "Otherinfo")

for (i in samples){
  if (nrow(snps)>0){
    snps.temp <- fread(paste0(i, '_SNP4_anno.hg38_multianno.csv'), sep=",", header=F, stringsAsFactors=F, skip = 1L)
    names(snps.temp) <- snp.header
    snps <- merge(snps, snps.temp, by=c(snp.header), all=T)
  }
  if (nrow(snps)==0){
    snps <- fread(paste0(i, '_SNP4_anno.hg38_multianno.csv'), sep=",", header=F, stringsAsFactors=F, skip=1L)
    names(snps) <- snp.header
  }
}

snp.hotspots <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/cancer_hotspots_SNPs.txt', header=F)
names(snp.hotspots) <- 'gene'


snps[,c('X1','X2','X3',"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late") := tstrsplit(Otherinfo, '\t')]

snps <- unique(snps)
mosaic.muts <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/hg38-iPSC-mosaic-damaging-mutations.txt', header=F, sep=':', sep2 = '-')
mosaic.muts[,c('start', 'end'):= tstrsplit(V2, '-')]
names(mosaic.muts)[1] <- 'chrom'
mosaic.muts[,V2 := NULL]
mosaic.muts
snps.mosaic <- subset(snps, snps$Chr %in% mosaic.muts$chrom)
snps.mosaic <- subset(snps.mosaic, snps.mosaic$Start %in% mosaic.muts$start) #0

snps.tp53 <- snps[Gene.refGene =='TP53',] #3

# rs1042522 : 'Homozygotes live 3 years longer. Chemotherapy is more effective. Increased breast cancer risk. https://www.snpedia.com/index.php/Rs1042522'


# only 'PASS' SNPs----
snps <- snps[FILTER =='PASS',] #74428

# subset data to either all SNPs or only those in hotspots----

#data <- snps.hotspots[,c('Chr', 'Start', 'End','Ref', 'Alt','avsnp150', "Gene.refGene", 'Func.refGene','1000g2015aug_all','ExAC_ALL',"Polyphen2_HDIV_pred",'Polyphen2_HDIV_score', 'cosmic70', "Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late")]

data.allsnps <- snps[,c('Chr', 'Start', 'End','Ref', 'Alt','avsnp150', "Gene.refGene", 'Func.refGene','ExonicFunc.refGene','1000g2015aug_all','ExAC_ALL',"Polyphen2_HDIV_pred",'Polyphen2_HDIV_score', 'cosmic70', "Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late")]

data <- data.allsnps

# simplify genotypes----
simplifygeno <- function(n) {
  df <- as.data.frame(n)
  
  genotype.cols <- length(df)-6
  
  
  for(i in names(df[,-c(1:genotype.cols)])){
    df[[i]] <- sub(':.*', '', df[[i]])
  }
  
  for(i in names(df[,c((genotype.cols+1):length(df))])){
    df[[i]] <- sub('0/1', 'het', df[[i]])
    df[[i]] <- sub('0/0', 'non', df[[i]])
    df[[i]] <- sub('1/0', 'het', df[[i]])
    df[[i]] <- sub('1/1', 'hom', df[[i]])
    df[[i]] <- sub('./0', 'miss', df[[i]])
    df[[i]] <- sub('0/.', 'miss', df[[i]])
    df[[i]] <- sub('1/.', 'miss', df[[i]])
    df[[i]] <- sub('./1', 'miss', df[[i]])
    df[[i]] <- sub('\\./\\.', 'miss', df[[i]])
    df[[i]] <- sub('\\.', 'miss', df[[i]])
    df[[i]] <- sub('0', 'non', df[[i]])
    df[[i]] <- sub('1', 'het', df[[i]])
  }
  x <- (length(df))
  
  df$miss <- rowSums(df == "miss")
  
  # remove rows with all missing (low quality) genotypes
  df <- subset(df, df$miss < 6)
  df <- df[,-(x+1)]
  
  # make sure there are no duplicate rows
  df <- unique(df)
  
  # done
  return(df)
}

data <- simplifygeno(data) #73571
data <- as.data.table(data)

# get counts of filtered variants----
snp.hotspots.list <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/cancer_hotspots_SNPs.txt', header=F)
snp.hotspots.list <- unique(snp.hotspots.list)
names(snp.hotspots.list) <- 'gene'
snps.in.hotspots <- c()

for (b in 1:nrow(snp.hotspots.list)){
  snps.in.hotspots <- c(snps.in.hotspots, grep(snp.hotspots.list$gene[b], data$Gene.refGene))
}

snps.hotspots <- unique(data[snps.in.hotspots,]) #4336
snps.hotspots <- snps.hotspots[Func.refGene=='exonic' & ExonicFunc.refGene=='nonsynonymous SNV',] #612

snps.hotspots.cosmic <- snps.hotspots[cosmic70 !='.',] #21
snps.hotspots.damaging <- snps.hotspots[Polyphen2_HDIV_pred=='D' & ((`1000g2015aug_all` <= 0.05) |(`1000g2015aug_all` == '.')),] #35

snps.hotspots.cosmic.damaging <- snps.hotspots.cosmic[Polyphen2_HDIV_pred=='D' & ((`1000g2015aug_all` <= 0.05) |(`1000g2015aug_all` == '.')),] #0

#View(snps.hotspots.cosmic.damaging)

snps.tp53 <- snps[Gene.refGene =='TP53',] #3

snps.cosmic <- snps[cosmic70 !='.',] #681
snps.cosmic.exonic.nonsyn.hotspot <- snps.hotspots[cosmic70 !='.' & Func.refGene=='exonic' & ExonicFunc.refGene =='nonsynonymous SNV',] #21

#

data.melt <- melt(data, measure.vars=c("Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late"))
data.melt
data.melt <- as.data.table(data.melt)
data.melt[,Passages := variable]
data.melt$Passages <- gsub('Lonza_P31_No_CEPT_early', 0, data.melt$Passages)
data.melt$Passages <- gsub('Lonza_P20_CEPT_late', 20, data.melt$Passages)
data.melt$Passages <- gsub('H9_P26_No_CEPT_early', 0, data.melt$Passages)
data.melt$Passages <- gsub('H9_P10_CEPT_mid', 10, data.melt$Passages)
data.melt$Passages <- gsub('H9_P20_CEPT_late', 20, data.melt$Passages)
data.melt$Passages <- gsub('Lonza_P10_CEPT_mid', 10, data.melt$Passages)

data.melt <- data.melt[value != '2/3' & value != '2/2',]

data.melt[,Genotype := value]
data.melt$Genotype <- gsub('het', 0.5, data.melt$Genotype)
data.melt$Genotype <- gsub('hom', 1, data.melt$Genotype)
data.melt$Genotype <- gsub('miss', 0, data.melt$Genotype)
data.melt$Genotype <- gsub('non', 0, data.melt$Genotype)
data.melt$Genotype <- as.numeric(data.melt$Genotype)

data.melt$Passages <- as.numeric(data.melt$Passages)

data.melt[,TGP:= as.numeric(`1000g2015aug_all`)]


#ggplot(data.melt, aes(x=Passages, y=value)) + geom_jitter()

#data.melt[,c('Polyphen2_HDIV_pred')]

data.melt[,Polyphen2 := Polyphen2_HDIV_pred]
data.melt$Polyphen2 <- gsub('B', '0', data.melt$Polyphen2)
data.melt$Polyphen2 <- gsub('P', '0.5', data.melt$Polyphen2)
data.melt$Polyphen2 <- gsub('D', '1', data.melt$Polyphen2)
data.melt$Polyphen2 <- gsub('//.', '0', data.melt$Polyphen2)
data.melt$Polyphen2 <- as.numeric(data.melt$Polyphen2)

data.melt[,Polyphen2_HDIV_score := as.numeric(Polyphen2_HDIV_score)]

data.melt[,COSMIC := ifelse(cosmic70 != '.', 1, 0)]

# saved data.melt containing all snps regardless of hotspot/cosmic presence--
#fwrite(data.melt, 'WES.SNPs.melted.csv', col.names = T, row.names = T, quote=T, sep=',')
#data.melt <- fread('WES.SNPs.melted.csv')

# exploratory correlation testing----
p.mat <- cor_pmat(data.melt[,c('Passages', 'Genotype', 'TGP', 'COSMIC', 'Polyphen2')])
#head(p.mat)
library(RColorBrewer)
ggcorrplot(p.mat, hc.order = TRUE, type = "lower", outline.color = "black", p.mat = p.mat, method='square', colors=rev(brewer.pal(n=3, name="RdBu")))

p.mat <- cor_pmat(data.melt[,c('Passages', 'Genotype', 'TGP', 'COSMIC', 'Polyphen2_HDIV_score')])
ggcorrplot(p.mat, hc.order = TRUE, type = "lower", outline.color = "black", p.mat = p.mat, method='square', colors=rev(brewer.pal(n=3, name="RdBu")))

ggplot(data=data.melt, aes(x=factor(Passages), y=TGP)) + geom_boxplot()
ggplot(data=data.melt, aes(x=Passages, y=Genotype, group=Passages )) + geom_boxplot()
ggplot(data=data.melt, aes(x=Passages, y=Polyphen2, group=Passages )) + geom_jitter()

library(ggpubr)
ggballoonplot(p.mat, fill = "value")+
  scale_fill_viridis_c(option = "C")


ggscatter(data.melt, x = "Passages", y = "Genotype",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Passages", ylab = "Genotype")

ggscatter(data.melt, x = "Passages", y = "Polyphen2",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Passages", ylab = "Polyphen2")

ggplot(data=data.melt, aes(x=factor(Passages), y=Genotype)) + geom_boxplot()



res.pca <- prcomp(na.omit(data.melt[,c('Passages', 'Genotype', 'TGP', 'COSMIC', 'Polyphen2')]))
fviz_eig(res.pca)

# correlation testing, different style graph----

data.melt[,`Cell line`:= variable]
data.melt$`Cell line` <- gsub('Lonza_P31_No_CEPT_early', 1, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('H9_P26_No_CEPT_early', 0, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('H9_P10_CEPT_mid', 0, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('H9_P20_CEPT_late', 0, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('Lonza_P10_CEPT_mid', 1, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('Lonza_P20_CEPT_late', 1, data.melt$`Cell line`)
data.melt$`Cell line` <- as.numeric(data.melt$`Cell line`)

rquery.cormat(data.melt[,c('Passages', 'Genotype', 'TGP', 'COSMIC', 'Polyphen2', 'Cell line')], type='upper')

M <- cor(na.omit(data.melt[,c('Passages', 'Genotype', 'TGP', 'COSMIC', 'Polyphen2', 'Cell line')]))

res1 <- cor.mtest(M, conf.level = .95)
res2 <- cor.mtest(M, conf.level = .99)

corrplot(M, p.mat = res1$p, sig.level = .2, type='upper', bg='lightgrey', col = rev(brewer.pal(11,"RdBu")))


#hotspots genotype barplot----
data.barplot <- data.melt %>% dplyr::group_by(Passages, Genotype) %>% dplyr::summarize(count = n())

#all snps version of barplot----
ggplot(data.barplot, aes(x=factor(Passages), y=count, fill=factor(Genotype), label=count)) + geom_bar(stat='identity')+ 
  theme_light() + theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=1, color='#e0e0e0'), panel.grid.major = element_line(size=1, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+
  geom_text(position = position_stack(vjust = 0.5), size=5)+
  labs(x='Passages', y='Frequency', title='Genotype frequencies by passage number for all SNPs') + scale_fill_manual(values=c('0' = "#deebf7", '0.5' = "#9ecae1", '1' = "#3182bd"), name = "Genotype", breaks=c("0", "0.5", "1"), labels=c("Noncarrier", "Heterozygous carrier", "Homozygous carrier"))

# cancer hotspots version of barplot----
data.barplot <- subset(data.melt, data.melt$Gene.refGene %in% snp.hotspots$gene) %>% dplyr::group_by(Passages, Genotype) %>% dplyr::summarize(count = n())

ggplot(data.barplot, aes(x=factor(Passages), y=count, fill=factor(Genotype), label=count)) + geom_bar(stat='identity')+ 
  theme_light() + theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=1, color='#e0e0e0'), panel.grid.major = element_line(size=1, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+
  geom_text(position = position_stack(vjust = 0.5), size=5)+
  labs(x='Passages', y='Frequency', title='Genotype frequencies by passage number for cancer hotspots') + scale_fill_manual(values=c('0' = "#deebf7", '0.5' = "#9ecae1", '1' = "#3182bd"), name = "Genotype", breaks=c("0", "0.5", "1"), labels=c("Noncarrier", "Heterozygous carrier", "Homozygous carrier"))

#
#tileplot----
data.melt.cosmic <- data.melt[COSMIC ==1,]
data.melt.cosmic <- data.melt.cosmic[Polyphen2_HDIV_pred !='.',]
 # ggplot(data.melt.cosmic, aes(x=factor(Passages), y=avsnp150, fill=factor(Polyphen2_HDIV_pred))) + geom_tile(width=.875, height=.875) +theme_light()+
 #   scale_fill_discrete(name='Polyphen2 \npredicted effect', breaks=c('B', 'P', 'D'), labels=c('Benign', 'Probably\ndamaging', 'Damaging'))

data.melt.cosmic[,Rarity := ifelse( (TGP <= 0.05 | is.na(TGP)), 'Rare', 'Common')]

cols <- c('Common-het'='#9ecae1','Common-hom'='#3182bd','Common-non'='#deebf7', 'Rare-het'='#fc9272','Rare-hom'='#de2d26','Rare-non'='#fee0d2')

data.melt.cosmic[,Polyphen2_HDIV_pred_ordered :=Polyphen2_HDIV_pred]
data.melt.cosmic$Polyphen2_HDIV_pred_ordered <- gsub('B', 'Benign', data.melt.cosmic$Polyphen2_HDIV_pred_ordered)
data.melt.cosmic$Polyphen2_HDIV_pred_ordered <- gsub('P', 'Possibly damaging', data.melt.cosmic$Polyphen2_HDIV_pred_ordered)
data.melt.cosmic$Polyphen2_HDIV_pred_ordered <- gsub('D', 'Probably damaging', data.melt.cosmic$Polyphen2_HDIV_pred_ordered)
data.melt.cosmic$Polyphen2_HDIV_pred_ordered = factor(data.melt.cosmic$Polyphen2_HDIV_pred_ordered, levels=c('Benign','Possibly damaging','Probably damaging'))

data.melt.cosmic[,Gene_rsID := paste0(Gene.refGene, ', ', avsnp150)]
data.melt.cosmic$Gene_rsID <- reorder(data.melt.cosmic$Gene_rsID, desc(data.melt.cosmic$Gene_rsID))

# 
ggplot(data.melt.cosmic, aes(x=factor(Passages), y=Gene_rsID, fill=paste0(factor(Rarity), '-', factor(value)))) + geom_tile() +theme_light()+
  scale_fill_manual(name='Rarity and genotype', breaks=c('Common-non', 'Common-het','Common-hom', 'Rare-non', 'Rare-het', 'Rare-hom'), labels=c('Common, noncarrier', 'Common, heterozygous carrier', 'Common, homozygous carrier', 'Rare, noncarrier', 'Rare, heterozygous carrier', 'Rare, homozygous carrier'), values=c(cols))+
  facet_grid(cols = vars(factor(data.melt.cosmic$Polyphen2_HDIV_pred_ordered)), drop=T, as.table = F) +
  labs(title='COSMIC mutations by passage number\nExonic, nonsynonymous, in cancer hotspots', x='Passages', y='SNP (avsnp150)')

#split tileplot by cell line AND passage----

data.melt.cosmic[,variable2 := factor(variable, levels=c("H9_P26_No_CEPT_early","H9_P10_CEPT_mid",  "H9_P20_CEPT_late", "Lonza_P31_No_CEPT_early", "Lonza_P10_CEPT_mid","Lonza_P20_CEPT_late"))]

ggplot(data.melt.cosmic, aes(x=factor(variable2), y=Gene_rsID, fill=paste0(factor(Rarity), '-', factor(value)))) + geom_tile() +theme_light()+
  scale_fill_manual(name='Rarity and genotype', breaks=c('Common-non', 'Common-het','Common-hom', 'Rare-non', 'Rare-het', 'Rare-hom'), labels=c('Common, noncarrier', 'Common, heterozygous carrier', 'Common, homozygous carrier', 'Rare, noncarrier', 'Rare, heterozygous carrier', 'Rare, homozygous carrier'), values=c(cols))+
  facet_grid(cols = vars(factor(data.melt.cosmic$Polyphen2_HDIV_pred_ordered)), drop=T, as.table = F) +
  labs(title='COSMIC mutations by passage number\nExonic, nonsynonymous', x='Cell line and passages', y='SNP (avsnp150)')+
  theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text.x = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=15), axis.title = element_text(size=15),axis.text.y=element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=0), panel.grid.major = element_line(size=0.5, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'),  plot.margin = margin(0,0,0, 2, "in"))

fwrite(data.melt.cosmic, 'COSMIC_all_tileplot_data.csv', quote=T, col.names = T, row.names = F)


# MANOVA ----

res.man <- manova(cbind(Passages, Genotype, Polyphen2) ~ `Cell line`, data = data.melt.cosmic)
summary(res.man)

# Kruskal-Wallis test-----
res.krs <- kruskal.test(Polyphen2 ~ `Cell line`, data = data.melt.cosmic)
summary(res.krs)
res.krs

res.krs <- kruskal.test(TGP ~ Genotype, data = data.melt.cosmic)
res.krs

res.krs <- kruskal.test(Polyphen2 ~ Passages, data = data.melt.cosmic)
res.krs

res.krs <- kruskal.test(Genotype ~ Passages, data = data.melt.cosmic)
res.krs

res.krs <- kruskal.test(TGP ~ Polyphen2, data = data.melt.cosmic)
res.krs

res.krs <- kruskal.test(COSMIC ~ `Cell line`, data = data.melt)
res.krs

# counts of cells per polyphen annotation----
h9.cosmic <- unique(data.melt.cosmic[`Cell line`==0 & (Genotype ==0.5 | Genotype ==1),])
lonza.cosmic <- unique(data.melt.cosmic[`Cell line`==1 & (Genotype ==0.5 | Genotype ==1),])

nrow(h9.cosmic[Polyphen2 ==1 & (Genotype ==0.5 | Genotype ==1),]) #86
nrow(h9.cosmic[Polyphen2 ==1 & (Genotype ==1),]) #6
nrow(h9.cosmic[Polyphen2 ==1 & (Genotype ==0.5),]) #80

nrow(lonza.cosmic[Polyphen2 ==1 & (Genotype ==0.5 | Genotype ==1),]) #92; 89
nrow(lonza.cosmic[Polyphen2 ==1 & (Genotype ==1),])  #16; 15
nrow(lonza.cosmic[Polyphen2 ==1 & (Genotype ==0.5),]) #76; 74

# look for zscan10----
anno <- fread('/Volumes/ncatssctl/NGS_related/WES/annotation/CT_ALL_SNP4_anno.hg38_multianno.csv')
str(anno)

str(data.melt)

# barplot of functional annotations by passage number----

data.barplot.0 <- data.melt[Passages==0 & (value =='het' | value =='hom'),]
data.barplot.10 <- data.melt[Passages==10 & (value =='het' | value =='hom'),]
data.barplot.20 <- data.melt[Passages==20 & (value =='het' | value =='hom'),]

data.barplot.0 <-  data.barplot.0 %>% dplyr::group_by(Passages, Func.refGene) %>% dplyr::summarize(count = n())
data.barplot.10 <-  data.barplot.10 %>% dplyr::group_by(Passages, Func.refGene) %>% dplyr::summarize(count = n())
data.barplot.20 <-  data.barplot.20 %>% dplyr::group_by(Passages, Func.refGene) %>% dplyr::summarize(count = n())

data.barplot.0 <- as.data.table(data.barplot.0)
data.barplot.10 <- as.data.table(data.barplot.10)
data.barplot.20 <- as.data.table(data.barplot.20)

data.barplot <- rbind(data.barplot.0, data.barplot.10, data.barplot.20)

#data.barplot <- data.melt %>% group_by(Passages, Func.refGene) %>% summarize(count = n())
data.barplot <- as.data.table(data.barplot)
ggplot(data.barplot, aes(x=reorder(Func.refGene, -count), y=count))+geom_bar(stat='identity')+theme_light()+facet_grid(rows=vars(factor(data.barplot$Passages)))+
  theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), axis.text.x=element_text(angle=45, hjust=1), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=1, color='#e0e0e0'), panel.grid.major = element_line(size=1, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+xlim('intronic','exonic','UTR3','intergenic','ncRNA_intronic','ncRNA_exonic','UTR5','upstream', 'downstream')+labs(x='Functional annotation',y='Count', title='Functional SNP annotations by passage number')


ggplot(data.barplot, aes(x=Passages, y=count, group=Func.refGene, color=Func.refGene)) + geom_line(size=1)+
  theme_light()+
  geom_text_repel(data = data.barplot[Passages==0,], aes(label = Func.refGene, colour = Func.refGene, x = 10, y = count+1000), segment.size=0.8, hjust = -.1, size=5)+
  scale_colour_discrete(guide = 'none')+
  theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=0.5, color='#e0e0e0'), panel.grid.major = element_line(size=0.5, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+
  labs(x='Passages',y='Count',title='Functional SNP annotations by passage number')+
  scale_x_continuous(breaks=c(0, 10,20))


fwrite(data.barplot, 'WES.barplot.data.csv', quote=T, col.names = T, row.names = F)

# exonic functional annotation barplot -----
data.barplot.0 <- data.melt[Passages==0 & (value =='het' | value =='hom'),]
data.barplot.10 <- data.melt[Passages==10 & (value =='het' | value =='hom'),]
data.barplot.20 <- data.melt[Passages==20 & (value =='het' | value =='hom'),]

data.barplot.0 <-  data.barplot.0 %>% dplyr::group_by(Passages, ExonicFunc.refGene) %>% dplyr::summarize(count = n())
data.barplot.10 <-  data.barplot.10 %>% dplyr::group_by(Passages, ExonicFunc.refGene) %>% dplyr::summarize(count = n())
data.barplot.20 <-  data.barplot.20 %>% dplyr::group_by(Passages, ExonicFunc.refGene) %>% dplyr::summarize(count = n())

data.barplot.0 <- as.data.table(data.barplot.0)
data.barplot.10 <- as.data.table(data.barplot.10)
data.barplot.20 <- as.data.table(data.barplot.20)

data.barplot <- rbind(data.barplot.0, data.barplot.10, data.barplot.20)

#data.barplot <- data.melt %>% group_by(Passages, Func.refGene) %>% summarize(count = n())
data.barplot <- as.data.table(data.barplot)

data.barplot <- as.data.table(data.barplot)

data.barplot <- data.barplot[ExonicFunc.refGene!='.',]

ggplot(data.barplot, aes(x=reorder(ExonicFunc.refGene, -count), y=count))+geom_bar(stat='identity')+theme_light()+facet_grid(rows=vars(factor(data.barplot$Passages)))+
  theme(legend.text=element_text(size=13), legend.title = element_text(size=15), axis.text = element_text(size=15), axis.text.x=element_text(angle=45, hjust=1), plot.title = element_text(size=15), axis.title = element_text(size=15), panel.border = element_rect(fill=NA, colour = "black", size=1), panel.grid.minor = element_line(size=1, color='#e0e0e0'), panel.grid.major = element_line(size=1, color='#e0e0e0'), axis.ticks = element_line(size=1, color='#e0e0e0'))+labs(x='Exonic functional annotation', y='Count', title='Exonic functional SNP annotations by passage number')


# INDELS-----------------------------------

setwd('/Volumes/ncatssctl/NGS_related/WES/annotation/')

anno.files <- Sys.glob('CT_ALL_INDEL*clinvar*csv')
anno.files

samples <- c('CT_ALL')

indels <- data.table()

indel.header <- c('Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'avsnp150', '1000g2015aug_all', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'cosmic70', 'clinvar_20150629','Otherinfo')

for (i in samples){
  if (nrow(indels)>0){
    indels.temp <- fread(paste0(i, '_INDEL4_anno_clinvar.hg38_multianno.csv'), sep=",", header=F, stringsAsFactors=F, skip = 1L)
    names(indels.temp) <- indel.header
    indels <- merge(indels, indels.temp, by=c(indel.header), all=T)
  }
  if (nrow(indels)==0){
    indels <- fread(paste0(i, '_INDEL4_anno_clinvar.hg38_multianno.csv'), sep=",", header=F, stringsAsFactors=F, skip=1L)
    names(indels) <- indel.header
  }
}





indel.hotspots <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/cancer_hotspots_INDELs.txt', header=F)
names(indel.hotspots) <- 'gene'


indels[,c('X1','X2','X3',"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late") := tstrsplit(Otherinfo, '\t')]

indels <- unique(indels) #10833


# only 'PASS' INDELs----
indels <- indels[FILTER =='PASS',] #8554

indels[,]



indel.hotspots <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/cancer_hotspots_INDELs.txt', header=F)
names(indel.hotspots) <- 'gene'

indels.tp53 <- indels[Gene.refGene =='TP53',] #1


indels[,c('CLINSIG', 'CLNDBN', 'CLNREVSTAT', 'CLNACC', 'CLNDSDB', 'CLNDSDBID'):= tstrsplit(clinvar_20150629, ';')]

indels[,c('CLINSIG')]

indels[,CLINSIG := tstrsplit(CLINSIG, "=")[2]]
unique(indels[,c('CLINSIG')])


indels.in.hotspots <- c()

for (b in 1:nrow(indel.hotspots)){
  indels.in.hotspots <- c(indels.in.hotspots, grep(indel.hotspots$gene[b], indels$Gene.refGene))
}

indels.hotspots <- indels[indels.in.hotspots]
indels.hotspots <- indels.hotspots[Func.refGene=='exonic' & (ExonicFunc.refGene=='frameshift insertion' | ExonicFunc.refGene =='frameshift deletion'),] #1

indels.hotspots <- unique(indels.hotspots) # there's only one variant here, chr6 152697965 152697968 

indels.hotspots.cosmic <- indels.hotspots[cosmic70 !='.',] #0
#indels.hotspots.damaging <- indels.hotspots[CLINSIG=='pathogenic' & ((`1000g2015aug_all` <= 0.05) |(`1000g2015aug_all` == '.')),]
indels.hotspots.rare <- indels.hotspots[((`1000g2015aug_all` <= 0.05) |(`1000g2015aug_all` == '.')),] #0
# 
# indels.hotspots.cosmic.rare <- indels.hotspots.cosmic[((`1000g2015aug_all` <= 0.05) |(`1000g2015aug_all` == '.')),]

#View(indels.hotspots.cosmic.damaging)

indels.pathogenic <- indels[CLINSIG=='pathogenic',] #2



# INDEL graphs for paper----

#data <- indels.hotspots[,c('Chr', 'Start', 'End','Ref', 'Alt','avsnp150', "Gene.refGene", 'Func.refGene','1000g2015aug_all','ExAC_ALL','cosmic70', 'CLINSIG', "Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late")]

data.allindels <- indels[,c('Chr', 'Start', 'End','Ref', 'Alt','avsnp150', "Gene.refGene", 'Func.refGene','ExonicFunc.refGene','1000g2015aug_all','ExAC_ALL','cosmic70', 'CLINSIG', "Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late")]

data <- data.allindels

# INDEL simplify genotypes----
data <- simplifygeno(data)
data

data.melt <- melt(data, measure.vars=c("Lonza_P31_No_CEPT_early", "H9_P26_No_CEPT_early", "H9_P10_CEPT_mid", "H9_P20_CEPT_late", "Lonza_P10_CEPT_mid", "Lonza_P20_CEPT_late"))
data.melt
data.melt <- as.data.table(data.melt)

data.melt <- data.melt[value =='het' | value == 'hom' | value=='miss' | value=='non',] #50715

length(unique(data.melt$avsnp150)) #7138 unqiue indels passing QC.

data.melt[,Passages := variable]
data.melt$Passages <- gsub('Lonza_P31_No_CEPT_early', 0, data.melt$Passages)
data.melt$Passages <- gsub('Lonza_P20_CEPT_late', 20, data.melt$Passages)
data.melt$Passages <- gsub('H9_P26_No_CEPT_early', 0, data.melt$Passages)
data.melt$Passages <- gsub('H9_P10_CEPT_mid', 10, data.melt$Passages)
data.melt$Passages <- gsub('H9_P20_CEPT_late', 20, data.melt$Passages)
data.melt$Passages <- gsub('Lonza_P10_CEPT_mid', 10, data.melt$Passages)

data.melt[,Genotype := value]
data.melt$Genotype <- gsub('het', 0.5, data.melt$Genotype)
data.melt$Genotype <- gsub('hom', 1, data.melt$Genotype)
data.melt$Genotype <- gsub('miss', 0, data.melt$Genotype)
data.melt$Genotype <- gsub('non', 0, data.melt$Genotype)
data.melt$Genotype <- as.numeric(data.melt$Genotype)

data.melt$Passages <- as.numeric(data.melt$Passages)

data.melt[,TGP:= as.numeric(`1000g2015aug_all`)]
data.melt$TGP[is.na(data.melt$TGP)] <- 0

data.melt[,CLINSIG_value := CLINSIG]
data.melt$CLINSIG_value <- gsub('non-pathogenic,non-pathogenic\\|untested', 0, data.melt$CLINSIG_value)
data.melt$CLINSIG_value <- gsub('non-pathogenic\\|non-pathogenic', 0, data.melt$CLINSIG_value)
data.melt$CLINSIG_value <- gsub('untested', 0, data.melt$CLINSIG_value)
data.melt$CLINSIG_value <- gsub('unknown', 0, data.melt$CLINSIG_value)
data.melt$CLINSIG_value <- gsub('probable-non-pathogenic', 0, data.melt$CLINSIG_value)
data.melt$CLINSIG_value <- gsub('non-pathogenic', 0, data.melt$CLINSIG_value)
data.melt$CLINSIG_value <- gsub('pathogenic', 1, data.melt$CLINSIG_value)
data.melt$CLINSIG_value[is.na(data.melt$CLINSIG_value)] <- 0
data.melt$CLINSIG_value <- as.numeric(data.melt$CLINSIG_value)

data.melt[,COSMIC := ifelse(cosmic70 != '.', 1, 0)]

data.melt[,`Cell line`:= variable]
data.melt$`Cell line` <- gsub('Lonza_P31_No_CEPT_early', 1, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('H9_P26_No_CEPT_early', 0, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('H9_P10_CEPT_mid', 0, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('H9_P20_CEPT_late', 0, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('Lonza_P10_CEPT_mid', 1, data.melt$`Cell line`)
data.melt$`Cell line` <- gsub('Lonza_P20_CEPT_late', 1, data.melt$`Cell line`)
data.melt$`Cell line` <- as.numeric(data.melt$`Cell line`)

subset_for_M <- na.omit(data.melt[,c('Passages', 'Genotype', 'TGP', 'COSMIC', 'CLINSIG_value', 'Cell line')])
names(subset_for_M)[5] <- 'CLINSIG'

M <- cor(subset_for_M)

library(corrplot)

res1 <- cor.mtest(M, conf.level = .95)
res2 <- cor.mtest(M, conf.level = .99)

corrplot(M, p.mat = res1$p, sig.level = .2, type='upper', bg='lightgrey', col = rev(brewer.pal(11,"RdBu")))


#

indels.writeout <- indels[avsnp150=='rs59758982' | avsnp150=='rs3841162',]
fwrite(indels.writeout, 'Indels-pathogenic.csv', sep=',', col.names = T, row.names = F, quote=T)
