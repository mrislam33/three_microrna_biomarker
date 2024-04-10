resdata = genes_module

####convert ensembl id to gene symbol
options(connectionObserver = NULL)
# suppressPackageStartupMessages(library(AnnotationDbi))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Mus.musculus")  ### put it to the starting packages

library("Mus.musculus")
library("org.Mm.eg.db")
library("clusterProfiler")
dir.create("figures")
dir.create("results")

resdata$symbol <- mapIds(Mus.musculus,
                         keys = resdata[,1],
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

resdata$entrez <- mapIds(Mus.musculus,
                         keys = resdata[,1],
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

gene.list=resdata



go_strong <- enrichGO(gene     = as.character(gene.list$entrez),
                         OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                         ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE, pool = FALSE)
go=data.frame(go_strong)


#based on
#http://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("rrvgo")

library(rrvgo)
go_analysis <- go
simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")



scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Mm.eg.db")


reduced_terms = go_analysis[match(go_analysis$ID, reducedTerms$go), ] #match the GO with the parental and original file
reduced_terms = na.omit(reduced_terms)
reducedTermsFull = data.frame(cbind(reducedTerms, reduced_terms)) #combine two files




term = data.frame(levels(as.factor(reducedTermsFull$parentTerm)))
#sumlog #combine p value by fisher p value
#sump #Combine p-values using the sum of p
#(Edgington's) method
#sumz #Combine p-values using the sum of z
#(Stouffer's) method
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)

result = vector()
fishersumMetap = vector()

for (i in term[,1]){
  analysis = reducedTermsFull[c(grep(i,reducedTermsFull$parentTerm)), ]
  print(analysis$parentTerm)
  analysis$fishersum = fishersMethod(analysis$p.adjust)
  result[i] <- analysis$fishersum[1]
  
}

print(result)


result_GO_fishersum = data.frame(result)
names(result_GO_fishersum) = "fishersumP"

result_GO_fishersum$bonferroniadjP = p.adjust(result_GO_fishersum$fishersumP, method = "bonferroni", n = length(result_GO_fishersum$fishersumP))
result_GO_fishersum$BHadjP = p.adjust(result_GO_fishersum$fishersumP, method = "BH", n = length(result_GO_fishersum$fishersumP))
result_GO_fishersum$fdradjP = p.adjust(result_GO_fishersum$fishersumP, method = "fdr", n = length(result_GO_fishersum$fishersumP))
result_GO_fishersum$neglogBonferroni = -log10(result_GO_fishersum$bonferroniadjP)

write.csv(result_GO_fishersum, paste0("./results/", initials, "_","parent_GO_BP_processes.csv"))##export parent term

result_GO_fishersum_order = result_GO_fishersum[match(reducedTermsFull$parentTerm, rownames(result_GO_fishersum)), ]

reducedTermsComplete = data.frame(cbind(reducedTermsFull,result_GO_fishersum_order))
write.csv(reducedTermsComplete, paste0("./results/",initials,  "_","detailed_GO_BP_processes.csv"))##export full details

pdf(paste0("./figures/",initials,  "_", "treemap_GO_BP.pdf"),  useDingbats = FALSE)
treemapPlot(reducedTermsComplete)
dev.off()



