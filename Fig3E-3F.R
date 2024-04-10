#####################################################################
############################ load required packages #################
#####################################################################
library(DESeq2)
library(RUVSeq)
library(ggplot2)
library(vsn)


##please take corresponding vignettes as reference for more details of the functions used

rm(list = ls())

#####################################################################
message("differential expression analysis related to miR-146a-5p data")
#####################################################################

###load required data
load("./data/datMeta_mir146.RData")


#datMeta contains data related to miR-146a-5p overexpression in 
#immortalized microglial (iMG) culture. 
#dataExpr contains raw microRNA counts in each sample
#condition contains information about the groups
#samples from both groups were equally distributed in different lanes for sequencing. 


### data preprocessing and differential expression analysis
zfGenes=dataExpr
nSamples=ncol(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=nSamples/2)#removing genes with low counts
filtered <- zfGenes[filter,]
genes <- rownames(filtered)
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(condition,  row.names=colnames(filtered)))

###differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = pData(set), design=~condition)
dds <- DESeq(dds)
rld <- rlog(dds)
meanSdPlot(assay(rld)) #plot mean SD plot



#####################################################################
message("PCA plot related to miR-146a-5p data")
#####################################################################
z=plotPCA(rld, intgroup = "condition", ntop = 500,  returnData = FALSE)
nudge <- position_nudge(y = 1)
print(z + theme_bw() )


### find the number of differentially expressed genes
res<- results(dds, contrast=c("condition","mimic","nc")) 
res <- res[order(res$padj), ]
normalized_counts=as.data.frame(counts(dds, normalized=TRUE))
diffexpr_mir146 <- merge(as.data.frame(res), normalized_counts, by="row.names", sort=FALSE)
names(diffexpr_mir146)[1] <- "Gene"

###number of signficantly deregulated genes
table(diffexpr_mir146$padj<0.05)[[2]]


#print(head(resdata))

#####################################################################
message("Volcano plot related to miR-146a-5p data")
#####################################################################

with(res, plot(log2FoldChange, -log10(padj), pch=20,  xlim=c( (min(res$log2FoldChange, na.rm=T)-1), (max(res$log2FoldChange, na.rm=T)+1) ), ylim=c( (min(-log10(res$padj), na.rm=T)-1), (max(-log10(res$padj), na.rm=T)+1) )))
with(subset(res, padj<0.05 & log2FoldChange >0), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, padj<0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
abline(h=-log10(0.05), lwd = 2, lty = 2)




#####################################################################
message("differential expression analysis related to miR-148a-3p data")
#####################################################################

###load required data
load("./data/datMeta_mir148.RData")
initials = "mir148"

#datMeta contains data related to miR-148a-3p overexpression in 
#primary hippocampal neuronal culture. 
#dataExpr contains raw microRNA counts in each sample
#condition contains information about the groups
#samples from both groups were equally distributed in different lanes for sequencing. 


### data preprocessing and differential expression analysis
zfGenes=dataExpr
nSamples=ncol(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=nSamples/2) #removing genes with low counts
filtered <- zfGenes[filter,]
genes <- rownames(filtered)
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(condition,  row.names=colnames(filtered)))

###differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = pData(set), design=~condition)
dds <- DESeq(dds)
rld <- rlog(dds)
meanSdPlot(assay(rld)) #plot mean SD plot




#####################################################################
message("PCA plot related to miR-148a-3p data")
#####################################################################

z=plotPCA(rld, intgroup = "condition", ntop = 500,  returnData = FALSE)
nudge <- position_nudge(y = 1)
print(z  + theme_bw() )


### find the number of differentially expressed genes
res<- results(dds, contrast=c("condition","mimic","control")) 
res <- res[order(res$padj), ]
normalized_counts=as.data.frame(counts(dds, normalized=TRUE))
diffexpr_mir148 <- merge(as.data.frame(res), normalized_counts, by="row.names", sort=FALSE)
names(diffexpr_mir148)[1] <- "Gene"
###number of signficantly deregulated genes
table(diffexpr_mir148$padj<0.05)[[2]]
#print(head(resdata))

#####################################################################
message("Volcano plot related to miR-148a-3p data")
#####################################################################

with(res, plot(log2FoldChange, -log10(padj), pch=20,  xlim=c( (min(res$log2FoldChange, na.rm=T)-1), (max(res$log2FoldChange, na.rm=T)+1) ), ylim=c( (min(-log10(res$padj), na.rm=T)-1), (max(-log10(res$padj), na.rm=T)+1) )))
with(subset(res, padj<0.05 & log2FoldChange >0), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, padj<0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
abline(h=-log10(0.05), lwd = 2, lty = 2)



#####################################################################
message("differential expression analysis related to miR-181a-5p data")
#####################################################################

###load required data
load("./data/datMeta_mir181.RData")
intials = "mir181"

#datMeta contains data related to miR-181a-5p overexpression in 
#primary hippocampal neuronal culture. 
#dataExpr contains raw microRNA counts in each sample
#condition contains information about the groups
#samples from both groups were equally distributed in different lanes for sequencing. 

### data preprocessing and differential expression analysis
zfGenes=dataExpr
nSamples=ncol(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=nSamples/2) #removing genes with low counts
filtered <- zfGenes[filter,]
genes <- rownames(filtered)
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(condition,  row.names=colnames(filtered)))


###differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = pData(set), design=~condition)
dds <- DESeq(dds)
res<- results(dds, contrast=c("condition","E","D"))
### E = miR-181a-overexpressed cells, D = scramble-RNA treated cells
rld <- rlog(dds)
meanSdPlot(assay(rld))  #plot mean SD plot


#####################################################################
message("PCA plot related to miR-181a-5p data")
#####################################################################
z=plotPCA(rld, intgroup = "condition", ntop = 500,  returnData = FALSE)
nudge <- position_nudge(y = 1)
print(z + theme_bw() )

### find the number of differentially expressed genes
res <- res[order(res$padj), ]
normalized_counts=as.data.frame(counts(dds, normalized=TRUE))
diffexpr_mir181 <- merge(as.data.frame(res), normalized_counts, by="row.names", sort=FALSE)
names(diffexpr_mir181)[1] <- "Gene"
###number of signficantly deregulated genes
table(diffexpr_mir181$padj<0.05)[[2]]
#print(head(resdata))



#####################################################################
message("Volcano plot related to miR-181a-5p data")
#####################################################################

with(res, plot(log2FoldChange, -log10(padj), pch=20,  xlim=c( (min(res$log2FoldChange, na.rm=T)-1), (max(res$log2FoldChange, na.rm=T)+1) ), ylim=c( (min(-log10(res$padj), na.rm=T)-1), (max(-log10(res$padj), na.rm=T)+1) )))
with(subset(res, padj<0.05 & log2FoldChange >0), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, padj<0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
abline(h=-log10(0.05), lwd = 2, lty = 2)


###########################################
### comparative GO analysis ########
###########################################

options(connectionObserver = NULL)
library("Mus.musculus")
library("org.Mm.eg.db")
library("clusterProfiler")


#mir181
####convert ensembl id to gene symbol
diffexpr_mir181$symbol <- mapIds(Mus.musculus,
                         keys = diffexpr_mir181$Gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

diffexpr_mir181$entrez <- mapIds(Mus.musculus,
                         keys = diffexpr_mir181$Gene,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

##separate up and downregulated genes 
gene.list=diffexpr_mir181
gene.list.sig=gene.list[gene.list$padj<0.05,]
gene.list.sig=na.omit(gene.list.sig)
gene.list.sig.fc=gene.list.sig[,c(1:3)]
gene.list.sig.fc$entrez=gene.list.sig$entrez
gene.list.sig.fc$group="upregulated"
gene.list.sig.fc$group[gene.list.sig.fc$log2FoldChange<0]="downregulated"
names(gene.list.sig.fc)[1]="genes"
gene.list.sig.up=gene.list.sig.fc[which(gene.list.sig.fc$group=="upregulated"),]
gene.list.sig.down=gene.list.sig.fc[which(gene.list.sig.fc$group=="downregulated"),]


#GO analysis
up_go_mir181 <- enrichGO(gene     = as.character(gene.list.sig.up$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)

down_go_mir181 <- enrichGO(gene     = as.character(gene.list.sig.down$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)


rm(gene.list, gene.list.sig.up, gene.list.sig.down)



###mir146
####convert ensembl id to gene symbol
diffexpr_mir146$symbol <- mapIds(Mus.musculus,
                         keys = diffexpr_mir146$Gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

diffexpr_mir146$entrez <- mapIds(Mus.musculus,
                         keys = diffexpr_mir146$Gene,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

##separate up and downregulated genes 
gene.list=diffexpr_mir146
gene.list.sig=gene.list[gene.list$padj<0.05,]
gene.list.sig=na.omit(gene.list.sig)
gene.list.sig.fc=gene.list.sig[,c(1:3)]
gene.list.sig.fc$entrez=gene.list.sig$entrez
gene.list.sig.fc$group="upregulated"
gene.list.sig.fc$group[gene.list.sig.fc$log2FoldChange<0]="downregulated"
names(gene.list.sig.fc)[1]="genes"
gene.list.sig.up=gene.list.sig.fc[which(gene.list.sig.fc$group=="upregulated"),]
gene.list.sig.down=gene.list.sig.fc[which(gene.list.sig.fc$group=="downregulated"),]

##GO analysis
up_go_mir146 <- enrichGO(gene     = as.character(gene.list.sig.up$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)

down_go_mir146 <- enrichGO(gene     = as.character(gene.list.sig.down$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)

rm(gene.list, gene.list.sig.up, gene.list.sig.down)





###########################################################

##mir148
####convert ensembl id to gene symbol
diffexpr_mir148$symbol <- mapIds(Mus.musculus,
                         keys = diffexpr_mir148$Gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

diffexpr_mir148$entrez <- mapIds(Mus.musculus,
                         keys = diffexpr_mir148$Gene,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

##separate up and downregulated genes 
gene.list=diffexpr_mir148
gene.list.sig=gene.list[gene.list$padj<0.05,]
gene.list.sig=na.omit(gene.list.sig)
gene.list.sig.fc=gene.list.sig[,c(1:3)]
gene.list.sig.fc$entrez=gene.list.sig$entrez
gene.list.sig.fc$group="upregulated"
gene.list.sig.fc$group[gene.list.sig.fc$log2FoldChange<0]="downregulated"
names(gene.list.sig.fc)[1]="genes"
gene.list.sig.up=gene.list.sig.fc[which(gene.list.sig.fc$group=="upregulated"),]
gene.list.sig.down=gene.list.sig.fc[which(gene.list.sig.fc$group=="downregulated"),]

##GO analysis
up_go_mir148 <- enrichGO(gene     = as.character(gene.list.sig.up$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)

down_go_mir148 <- enrichGO(gene     = as.character(gene.list.sig.down$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)

rm(gene.list, gene.list.sig.up, gene.list.sig.down)




##plot 3F

# comparative GO plotting for upregulated terms

plot_up=merge_result(list(miR_146a_up=up_go_mir146,
                          miR_148a_up=up_go_mir148,
                          miR_181a_up=up_go_mir181))

how.many.items = 5  #plot top 5 from each condition

dotplot(object = plot_up, x = ~Cluster,
                  color = "p.adjust", showCategory = how.many.items, split = NULL,
                  font.size = 12, title = "", by = "count", includeAll = TRUE)



# comparative GO plotting for downregulated terms

plot_down=merge_result(list(miR_146a_down=down_go_mir146,
                          miR_148a_down=down_go_mir148,
                          miR_181a_down=down_go_mir181))

how.many.items = 5 #plot top 5 from each condition

dotplot(object = plot_down, x = ~Cluster,
        color = "p.adjust", showCategory = how.many.items, split = NULL,
        font.size = 12, title = "", by = "count", includeAll = TRUE)



