

#this script runs WGCNA analysis and find modules in aging mice model treated with 3-miR inhibitor
#and scramble control. Same analysis is also used to analyze the data from APP/PS1 mice only with a change in softpower as described in Methods and Materials. 

options(stringsAsFactors  =  FALSE)
library("edgeR")
library(WGCNA)
library(flashClust)
library(AnnotationDbi)
library("Mus.musculus")
library("org.Mm.eg.db")
library("clusterProfiler")
library(rrvgo)



load("./data/datMeta_agingWGCNA_rawData.Rdata")

#condition contains the numeric annotation of the experimental groups; 1:young_nc; #2:old_nc #3:old_inhibitor
#details contain the alphabetical annotation of the experimental groups; A:young_nc; #B:old_nc #C:old_inhibitor

###data preprocessing
zfGenes=dataExpr
nSamples=ncol(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=nSamples/2) #remove genes with low counts
filtered <- zfGenes[filter,]
number=as.numeric(nrow(filtered))

genes <- rownames(filtered)


## Run the regression and make the sequencing lane adjusted normalized matrix
regvars <- as.data.frame(cbind(sequencing_lane, condition))

##library size normalization
logCPM=cpm(filtered, normalized.lib.sizes=TRUE, log=FALSE)
for (i in 1:dim(logCPM)[1]) {
  logCPM[i, ]=sapply(logCPM[i, ], function(x) log2(1+x))
}



## regressed expression matrix
datExpr.reg <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(logCPM))
rownames(datExpr.reg) <- rownames(logCPM)
colnames(datExpr.reg) <- colnames(logCPM)
coefmat <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(regvars)+1)


for (i in 1:nrow(logCPM)) {
  lmmod1 <- lm(as.numeric(logCPM[i,])~sequencing_lane + condition,data=regvars)
  coef <- coef(lmmod1)
  datExpr.reg[i,] <- logCPM[i,]  - coef["sequencing_lane"]*regvars[,"sequencing_lane"]
}
datExpr_reg <- datExpr.reg


datNorm = data.frame(datExpr_reg) ##normalized coutns

#save files
save(datNorm, detail,  file ="./data/datMeta_agingWGCNA.Rdata")



load("./data/datMeta_agingWGCNA.Rdata") 


datExpr = as.data.frame(t(datNorm))
datExpr[datExpr == 0] <- 1
datExpr=na.omit(datExpr)
datTraits = detail
names(datTraits) = "condition_short"
datTraits$group = c(rep("aging_young_nc", 6), rep("aging_old_nc", 6), rep("aging_old_inhibitor", 7))
datTraits$group_numeric = as.numeric(as.factor(datTraits$group))
datatrait_numeric=sapply(datTraits, function (x) as.numeric(x))
datatrait_numeric=datatrait_numeric[,-c(1,2)]
rownames(datExpr) = rownames(datTraits)
table(rownames(datTraits)==rownames(datExpr)) 


A = adjacency(t(datExpr),type="signed") 
k = as.numeric(apply(A,2,sum))-1 
Z.k = scale(k)
thresholdZ.k = -2.5 
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
traitColors = data.frame(numbers2colors(datatrait_numeric,signed=FALSE))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram Heatmap")



allowWGCNAThreads() 


# Pick soft threshold
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr,networkType = "signed", corFnc = "bicor",verbose = 5,powerVector = powers)
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=1.0
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

softPower = 17  ## Choose based on fit to scale-free topology at R2 > 0.8
adjacency = adjacency(datExpr, corFnc = "bicor", type = "signed", power = softPower)

###to minimize the time to run the code the output of the following
##silent code can be accessed using load function. 

# TOM = TOMsimilarity(adjacency,TOMType = "signed", verbose = 0)
# dissTOM = 1-TOM

# save(TOM, dissTOM, file ="./data/network_TOMaging.RData")

load("./data/network_TOMaging.RData")

geneTree = flashClust(as.dist(dissTOM), method = "average")




################### Iterate WGCNA parameters #########
# colors <- moduleLabel  <- c()
# 
# for (minModuleSize in c(50,100)) {
#   for (ds in c(1:4)) {
#     
#     cutHeight = 0.99999
#     dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
#                                 deepSplit = ds, pamRespectsDendro = T,pamStage = T,
#                                 minClusterSize = minModuleSize, cutHeight = cutHeight)
#     
#     
#     merged <- mergeCloseModules(exprData = datExpr,colors = labels2colors(dynamicMods))
#     
#     
#     colors = cbind(colors, labels2colors((merged$colors)))
#     moduleLabel <- c(moduleLabel,paste("DS=",ds,
#                                        " mms=\n",minModuleSize))
#     
#     #plot dendrotree
#     plotDendroAndColors(geneTree,colors,
#                         groupLabels=moduleLabel,
#                         addGuide=FALSE,
#                         dendroLabels=FALSE,
#                         main="Dendrogram With Different Module Cutting Parameters",
#                         cex.colorLabels=0.5)
#     
#   }}
# 
# dim = 6 #selected parameter
# 
# c_ref = as.character(colors[,dim]) #module size 100, ds 2
# for (i in 1:dim(colors)[[2]]){
#   ci = as.character(colors[,i])
#   c_new = matchLabels(ci, c_ref)
#   colors[,i] = c_new
# }
# 
# colors = cbind(colors[,dim], colors)
# moduleLabel = c("selected module", moduleLabel)
# 
# plotDendroAndColors(geneTree,colors,groupLabels = moduleLabel,addGuide=T,dendroLabels=F,cex.colorLabels=0.3)




###Run WGCNA with selected parameters

minModuleSize =100; ds = 2;cutHeight = 0.99999
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                            deepSplit = ds, pamRespectsDendro = T,pamStage = T,
                            minClusterSize = minModuleSize, cutHeight = cutHeight)
dynamicColors = labels2colors(dynamicMods)


#calculate eigevalue
MEList = moduleEigengenes(datExpr, colors = dynamicColors,softPower = softPower)
MEs = MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#merge similar clusters
MEDissThres = 0.15
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors = merge$colors
MEs = merge$newMEs


# Calculate module membership for each module

modNames=substring(names(MEs),3)
geneModulecor=corAndPvalue(datExpr, MEs, use = "p")
geneModuleMembership = as.data.frame(geneModulecor$cor)
MMPvalue = as.data.frame(geneModulecor$p)
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneInfo = data.frame(geneID = colnames(datExpr), moduleColor = moduleColors, geneModuleMembership, MMPvalue)



#convert ensembl ID to gene symbols
geneInfo$symbol <- mapIds(Mus.musculus,
                          keys = geneInfo$geneID,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")


##preparing files to plot 7B
group_mice = data.frame(c(rep("3m_nc", 6), rep("16m_nc", 6), rep("16m_in",7)))
data_eigengene = data.frame(cbind(datTraits$group,group_mice[,1],  MEs))
names(data_eigengene)[2] = "plot_pheno.group"

comparison = list(c("3m_nc","16m_nc"),   c("3m_nc", "16m_in"))


#plot 7B

ggboxplot(data_eigengene, x = "plot_pheno.group", y = "MEblue",
ylab = "eigenvalue",
          title = "MEblue") + stat_compare_means(comparisons = comparison) 


###blue module gene ontology
blue_module = geneInfo[which(geneInfo$moduleColor == "blue"),]
genes_module = data.frame(blue_module[,"geneID"])



genes_module$symbol <- mapIds(Mus.musculus,
                         keys = genes_module[,1],
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

genes_module$entrez <- mapIds(Mus.musculus,
                         keys = genes_module[,1],
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

#GO analysis
go_analysis <- enrichGO(gene     = as.character(genes_module$entrez),
                      OrgDb = 'org.Mm.eg.db', keyType = "ENTREZID",
                      ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                      readable = TRUE, pool = FALSE)

#calculate semantic simlarity
simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)


#group GO terms based on similarity
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Mm.eg.db")
reduced_terms = go_analysis[match(go_analysis$ID, reducedTerms$go), ] 
#match the reduced parental GO terms with the initial GO analysis file from line 265-268

reduced_terms = na.omit(reduced_terms)
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
result_GO_fishersum$neglogBonferroni = -log10(result_GO_fishersum$bonferroniadjP)

write.csv(result_GO_fishersum, "GO_analysis_blue_module_results.csv")





