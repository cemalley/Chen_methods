library(enrichR)
library(data.table)

# Query a gene list against databases for enrichment analysis, e.g. GO Biological Pathways
# Thanks to John Braisted for the following code
# Claire Malley
# NCATS NIH


combineResults <-function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                        oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  print(length(libs))
  for (i in 1:length(libs)) {
    res<-data.frame(enr[i])
    print(dim(res)[1])
    for(j in 1:(dim(res)[1])) {
      currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
      names(currdf) <- col_names
      fullRes <- rbind(fullRes, currdf)
    }
  }
  return(fullRes)
}
getSigTerms <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullSigRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                           oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  print(length(libs))
  for (i in 1:length(libs)) {
    res<-data.frame(enr[i])
    print(dim(res)[1])
    for(j in 1:(dim(res)[1])) {
      if(res[j,4]<0.1) {				
        currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
        names(currdf) <- col_names
        fullSigRes <- rbind(fullSigRes, currdf)
      }
    }
  }
  return(fullSigRes)
}
getTop10SigTerms <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  topSigRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                          oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  for (i in 1:length(libs)) {
    #loop over result, append lib name to front, save only adjP <0.1
    res<-data.frame(enr[i])
    rowCnt<-dim(res)[2]
    # only consider the top 10 scores, only take adjP < 0.1
    for(j in 1:min(10,rowCnt)) {
      print(paste("**",res[j,4],"**"))
      if((!is.na(res[j,4])) & (res[j,4]<0.1)) {
        currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
        names(currdf) <- col_names
        topSigRes <- rbind(topSigRes, currdf)
      }
    }
  }
  return(topSigRes)
}
unpivotGeneLists <- function(fsr) {
  geneFrame<- data.frame(libName=character(), term=character(), adjP=double(), comb_score=double(), gene_symbol=character(), stringsAsFactors=FALSE)
  currRow <- 1
  mainRow <- 1
  for(r in 1:(dim(fsr)[1])){
    genes <- unlist(strsplit(as.character(fsr[r,"gene_list"]),';'))		
    for(gene in 1:length(genes)){
      geneFrame[mainRow,1]<-as.vector(fsr[currRow,"libName"])
      geneFrame[mainRow,2]<-as.vector(fsr[currRow,"term"])
      geneFrame[mainRow,3]<-fsr[currRow,"adjPval"]
      geneFrame[mainRow,4]<-fsr[currRow,"score"]
      geneFrame[mainRow,5]<-genes[gene]
      mainRow <- mainRow + 1
    }
    currRow <- currRow + 1
  }
  return(geneFrame)
}

libFile <- file.choose()
geneFile <- file.choose()

outputSigFile <- paste(gsub(".txt","",geneFile),"_SIG_EnrichR_Output.csv", sep=",")
outputSigUnpivotFile<- paste(gsub(".txt","",geneFile),"_SIG_EnrichR_UNPIVOT.csv", sep=",")
outputTopFile <- paste(gsub(".txt","",geneFile),"_TOP10_EnrichR_Output.csv", sep=",")
outputTopUnpivot<- paste(gsub(".txt","",geneFile),"_Top10_EnrichR_UNPIVOT.csv", sep=",")

geneSet <- read.table(geneFile, header=FALSE)
libSet <- read.table(libFile, header=FALSE)
genes<-as.vector(geneSet[,1])
libs<-as.vector(libSet[,1])

enr<-enrichr(genes, libs)

fsr <- getSigTerms(enr)

if(dim(fsr)[1] > 0) {
  fwrite(fsr, outputSigFile)
  #fsrTop10 <- getTop10SigTerms(enr)
  #fwrite(fsrTop10, outputTopFile)
} else {
  fwrite("No Significant Enrichments", outputSigFile)
  fwrite("No Significant Enrichments", outputTopFile)
}

