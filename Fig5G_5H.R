###Figure 5G-5H

##note: This script analyzes sequencing data from Bioengineered brain organoids (BENOS)treated with 3-miR-mimic-mix or scramble RNA
options(warn = -1)

options(connectionObserver = NULL)
library("Homo.sapiens")
library("org.Hs.eg.db")
library("clusterProfiler")
library(rrvgo)
library(DESeq2)
library(RUVSeq)
library(vsn)

#rm(list = ls())


load("./data/datMeta_BENO.RData")

#datMeta contains data related to 3-miR-mix overexpression for 24 hours in 
#Bioengineered brain organoids (BENOS)
#dataExpr contains raw gene counts in each sample
#condition contains information about the groups
#samples from both groups were equally distributed in different lanes for sequencing
#to avoid technical effects

zfGenes=dataExpr
nSamples=ncol(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=nSamples/2) #remove genes with low reads.
filtered <- zfGenes[filter,]

genes <- rownames(filtered)

set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(condition,  row.names=colnames(filtered)))

##differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = pData(set), design=~condition)
dds <- DESeq(dds)



meanSdPlot(assay(rld)) #plot mean SD PLOT


res<- results(dds, contrast=c("condition","B","A")) 
#B = 3-miR mimic mix treated cells, A = scramble RNA treated cells
res <- res[order(res$padj), ]
resdata <- data.frame(res)


#volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20,  xlim=c( (min(res$log2FoldChange, na.rm=T)-1), (max(res$log2FoldChange, na.rm=T)+1) ), ylim=c( (min(-log10(res$padj), na.rm=T)-1),7 ))) #(max(-log10(res$padj), na.rm=T+1)
with(subset(res, padj<0.05 & log2FoldChange >0), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, padj<0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=20, col="turquoise"))
abline(h=-log10(0.05), lwd = 2, lty = 2)


####convert ensembl id to gene symbol
resdata$symbol <- mapIds(Homo.sapiens,
                         keys = rownames(res),
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

resdata$entrez <- mapIds(Homo.sapiens,
                         keys = rownames(res),
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

gene.list=resdata

##prepare the up- and down-regulated gene list for GO analysis
gene.list.sig=gene.list[gene.list$padj<0.05,]
gene.list.sig=na.omit(gene.list.sig)
gene.list.sig.fc=gene.list.sig[,c(1:3)]
gene.list.sig.fc$entrez=gene.list.sig$entrez
gene.list.sig.fc$group="upregulated"
gene.list.sig.fc$group[gene.list.sig.fc$log2FoldChange<0]="downregulated"
names(gene.list.sig.fc)[1]="genes"
gene.list.sig.up=gene.list.sig.fc[which(gene.list.sig.fc$group=="upregulated"),]
gene.list.sig.down=gene.list.sig.fc[which(gene.list.sig.fc$group=="downregulated"),]



#
##analysis of upregulated genes as example ####
### similar analysis has been performed for downregulated genes #######
##################################
#GO analysis for upregulated genes
up_go_strong <- enrichGO(gene     = as.character(gene.list.sig.up$entrez),
                         OrgDb = 'org.Hs.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)

#calculate semantic simlarity
simMatrix <- calculateSimMatrix(up_go_strong$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(up_go_strong$qvalue), up_go_strong$ID)

#group GO terms based on similarity
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Hs.eg.db")


reduced_terms = go_analysis[match(go_analysis$ID, reducedTerms$go), ] #match the reduced parental GO terms with the initial GO analysis file from line 91-94
reduced_terms = na.omit(reduced_terms) #to remove the terms with NAs, normally not needed step
reducedTermsFull = data.frame(cbind(reducedTerms, reduced_terms)) #combine two files


#calculate cumulative group p value for the parent term considering all the member GO terms within the given parental term  by fisher method

term = data.frame(levels(as.factor(reducedTermsFull$parentTerm)))
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)

result = vector()
fishersumMetap = vector()

for (i in term[,1]){
  analysis = reducedTermsFull[c(grep(i,reducedTermsFull$parentTerm)), ]
  print(analysis$parentTerm)
  analysis$fishersum = fishersMethod(analysis$p.adjust)
  result[i] <- analysis$fishersum[1]
  
}

result_GO_fishersum = data.frame(result)
names(result_GO_fishersum) = "fishersumP"

##adjust the group p value with multiple adjustments
result_GO_fishersum$bonferroniadjP = p.adjust(result_GO_fishersum$fishersumP, method = "bonferroni", n = length(result_GO_fishersum$fishersumP))


##similar analyses was done for the downregulated genes
