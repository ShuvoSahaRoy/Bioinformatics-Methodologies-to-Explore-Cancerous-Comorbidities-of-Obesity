# 1. Set your work folder
setwd("C:/Users/SSROY/Desktop/bio/all datasets/thyroid")

# 2. Load libraries for script one
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65144/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }
# 
# # 4. Convert download dataset in usable class
# gse <- getGEO(filename = "GSE65144_series_matrix.txt.gz")

# 3. Download the GEO dataset only if not downloaded
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65144/matrix/"
filename <- "GSE65144_series_matrix.txt.gz"

# If file is missing, then download
if (!file.exists(filename)) {
  dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  dataset <- unlist(strsplit(dataset, "\r\n"))
  for (ds in dataset) {
    download.file(paste0(url, ds),
                  paste0(getwd(), "/", ds))
  }
}

# 4. Convert download dataset in usable class
gse <- getGEO(filename = filename)

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CD',12),rep('CTRL',13)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv2(subtable_result,"GSE65144_table.csv")
nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
CD_GOdata <- new("topGOdata",
                 description = "cd study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 5)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(CD_GOdata)
ug <- usedGO(CD_GOdata)

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(CD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(CD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
pvalFis <- score(resultFisher)
pvalKS <- score(resultKS,whichGO = names(pvalFis))
cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(CD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 22)
showSigOfNodes(CD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(CD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "TC2_GSE65144", useInfo = "all", pdfSW = TRUE)

# 14. Create text file for the correspondence GO terms - genes (this file is mandatory for the script two)
terms <- allRes$GO.ID
genes <- genesInTerm(CD_GOdata,terms)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "TC2_GSE65144_correspondence.txt" )
}
