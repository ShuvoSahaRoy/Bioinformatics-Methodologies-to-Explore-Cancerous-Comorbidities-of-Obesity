
# 1. Set your work folder
setwd("C:/Users/SSROY/Desktop/bio")

# 2. Load libraries for script one
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

#3. Download the GEO dataset 
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60403/matrix/"
dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
dataset <- unlist(strsplit(dataset, "\r\n"))
for (ds in dataset) {
  download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
}
 
# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE60403_series_matrix.txt.gz")

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 8),rep('OB',8)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCTRL-dOB, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv2(subtable_result,"GSE60403_table(1).csv")
nrow(subtable_result[subtable_result$P.Value<0.05,])

write.table(table_result,file="GSE60403 table.csv",sep = ",",row.names = FALSE)
write.table(subtable_result,file="GSE60403subtable.csv",sep = ",",row.names = FALSE)


# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
OB_GOdata <- new("topGOdata",
                 description = "OB study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 5)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(OB_GOdata)
ug <- usedGO(OB_GOdata)

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(OB_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(OB_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
pvalFis <- score(resultFisher)
pvalKS <- score(resultKS,whichGO = names(pvalFis))
cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(OB_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 22)
showSigOfNodes(OB_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(OB_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "OB1_GSE60403", useInfo = "all", pdfSW = TRUE)

# 14. Create text file for the correspondence GO terms - genes (this file is mandatory for the script two)
terms <- allRes$GO.ID
genes <- genesInTerm(OB_GOdata,terms)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "OB1_GSE60403_correspondence.txt" )
}
