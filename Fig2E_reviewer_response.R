####longitudinal expression analysis between control 
library(RUVSeq)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(vsn)

load("./data/longiAging.RData")
zfGenes=datExpr
filter <- apply(zfGenes, 1, function(x) length(x[x>1])>=20)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(time, row.names=colnames(filtered)))
#####remove unwanted variation using replicate samples
differences=makeGroups(time)
set2=RUVs(set, genes, k=1, differences)


dds <- DESeqDataSetFromMatrix(countData = counts(set2), colData = pData(set2), design=~W_1+time)
ddsTC <- DESeq(dds, test="LRT", fitType = 'local', reduced = ~W_1+1)


rld <- rlog(ddsTC, blind=FALSE)
meanSdPlot(assay(rld))

resTC <- results(ddsTC)
table(resTC$padj<0.05)
diff_results_homecage=data.frame(resTC)
diff_results_homecage_sig = diff_results_homecage[diff_results_homecage$padj<0.05,]
diff_results_homecage_sig_exp_100 = diff_results_homecage_sig[diff_results_homecage_sig$baseMean>100,]

betas_aging <- data.frame(coef(ddsTC))
log2fc_aging <- betas_aging
log2fc_aging$Intercept = NULL; log2fc_aging$W_1 = NULL


load("./data/longiAging2.RData")

zfGenes=datExpr_learning
filter <- apply(zfGenes, 1, function(x) length(x[x>1])>=16)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(group, row.names=colnames(filtered)))
differences=makeGroups(group)
set2=RUVs(set, genes, k=1, differences)
dds <- DESeqDataSetFromMatrix(countData = counts(set2), colData = pData(set2), design=~W_1 + group)
ddsTC <- DESeq(dds, test="LRT",fitType = 'local', reduced = ~W_1 + 1)
rld <- rlog(ddsTC, blind=FALSE)
meanSdPlot(assay(rld))

resTC <- results(ddsTC)
table(resTC$padj<0.05)
diff_results_learning=data.frame(resTC)
diff_results_learning_sig=diff_results_learning[diff_results_learning$padj<0.05,]
diff_results_learning_sig_exp_100 = diff_results_learning_sig[diff_results_learning_sig$baseMean>100,]


message("common microRNAs:")
print(intersect(rownames(diff_results_homecage_sig), rownames(diff_results_learning_sig)))
print(intersect(rownames(diff_results_homecage_sig_exp_100), rownames(diff_results_learning_sig_exp_100)))


#########
common_mirnas = data.frame(intersect(rownames(diff_results_homecage_sig_exp_100), rownames(diff_results_learning_sig_exp_100)))
names(common_mirnas) = "microRNA"
plot = log2fc_aging[match(common_mirnas$microRNA, rownames(log2fc_aging)),]
plot = as.matrix(plot)
my_palette <- colorRampPalette(c("darkblue","blue","lightblue","white","orange", "red", "darkred"))(n = 256)

heatmap.2(plot, density.info="none", Rowv=TRUE, 
          Colv=FALSE, dendrogram="row", trace="none",
          labRow=rownames(plot),
          labCol=colnames(plot), 
          col=my_palette, 
          cexCol = 0.5, 
          cexRow = 0.20, srtCol = 50,
          margins=c(2,3))

