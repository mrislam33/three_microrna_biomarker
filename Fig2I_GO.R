rm(list=ls())
library(clusterProfiler)
library(AnnotationDbi)
library(ReactomePA)
library(rrvgo)
library("Mus.musculus")



####get the predicted genes
mir146_target = read.csv("./data/TargetScan7.1__miR-146-5p_predicted_targets.csv",
                         sep = ",")

mir148_target = read.csv("./data/TargetScan7.1__miR-148-3p_predicted_targets.csv",
                         sep = ",")

mir181_target = read.csv("./data/TargetScan7.1__miR-181-5p_predicted_targets.csv",
                         sep = ",")

mir = data.frame(rbind(mir146_target, mir148_target,mir181_target)) #combine the data


###convert gene symbol to entrezids
mir$entrez <- mapIds(Mus.musculus,
                     keys=as.character(mir$Target.gene),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")


predicted_targets = mir #contain the predicted targets for all three microRNAs

## GO analysis
go_strong_enrich <- enrichGO(gene = as.character(predicted_targets$entrez),
                             OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                             ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                             readable = TRUE, pool = FALSE)

#calculate semantic simlarity
simMatrix <- calculateSimMatrix(go_strong_enrich$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_strong_enrich$qvalue), go_strong_enrich$ID)

#group GO terms based on similarity
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Mm.eg.db")

reduced_terms = go_strong_enrich[match(go_strong_enrich$ID, reducedTerms$go), ] #match the reduced parental GO terms with the initial GO analysis file from line 34-37
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


##similar analyses were performed for the predicted genes that overlapped with syngo and inflammatory genes
##dataset related to syngo can be accessed here ("./data/syngo_genes.xlsx")
##dataset related to inflammatory genes can be accessed here ("./data/GO_inflammatory_response_GO-0006954.csv")

