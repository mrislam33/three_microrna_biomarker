#####Figure 1C related

rm(list = ls())

library(edgeR)
library(WGCNA)
library(flashClust)
options(stringsAsFactors  =  FALSE)
library(clusterProfiler)

list.files()

load("./data/Fig1datMeta.Rdata")

dir.create("Fig1C_related")
setwd("Fig1C_related")

###removal of row reads, normalization to the library size
zfGenes = datExpr

nSamples=ncol(zfGenes)
minimum_reads=5
filter <- apply(zfGenes, 1, function(x) length(x[x>minimum_reads])>=nSamples/2) #remove genes low read counts
filtered <- zfGenes[filter,]

##normalization
logCPM=cpm(filtered, normalized.lib.sizes=TRUE, log=FALSE)
for (i in 1:dim(logCPM)[1]) {
  logCPM[i, ]=sapply(logCPM[i, ], function(x) log2(1+x))
}
logCPM = as.matrix(logCPM)


## Run the regression and make the gender adjusted normalized matrix
sex <- as.numeric(as.factor(datTrait[,"sex"]))-1
age=as.numeric(datTrait[,"age"])
regvars <- as.data.frame(cbind(age,sex))

datExpr.reg <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(logCPM))
rownames(datExpr.reg) <- rownames(logCPM)
colnames(datExpr.reg) <- colnames(logCPM)
coefmat <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(regvars)+1)


for (i in 1:nrow(logCPM)) {
  lmmod1 <- lm(as.numeric(logCPM[i,])~age+sex,data=regvars)
  coef <- coef(lmmod1)
  datExpr.reg[i,] <- logCPM[i,]  - coef["sex"]*regvars[,"sex"] 
}

datExpr_reg <- datExpr.reg #regressed expression matrix


# Manipulate file format so it matches the format WGCNA needs 

datExpr = as.data.frame(t(datExpr_reg)) # now samples are rows and genes are columns
dim(datExpr) # 132 samples and 456 genes 
rownames(datExpr) = colnames(datExpr_reg)
datExpr[datExpr == 0] <- 1 #replace counts 0 to 1


# check of gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK 

filtered_datTraits = datTrait
filtered_datTraits$age = NULL
datatrait_numeric=sapply(filtered_datTraits, function (x) as.numeric(x))




#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datExpr) = rownames(filtered_datTraits)

table(rownames(filtered_datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly



##Clustering of samples
A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datatrait_numeric,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(filtered_datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),cex.dendroLabels = 0.5, cex.rowText = 0.3,
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")


#####find and remove outlying samples##
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr=datExpr[!remove.samples, ]
filtered_datTraits=filtered_datTraits[!remove.samples,]
datatrait_numeric=sapply(filtered_datTraits, function (x) as.numeric(x))


###WGCNA analysis

allowWGCNAThreads()

# Pick soft threshold
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr,networkType = "signed", corFnc = "bicor",verbose = 5,powerVector = powers)
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=1.0

#pdf("./Fig1C_related/soft_threshold.pdf", useDingbats = FALSE)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
#dev.off()

#pdf("./Fig1C_related/mean_connectivity.pdf", useDingbats = FALSE)

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
#dev.off()

softPower = 9 # Choose based on fit to scale-free topology at R = 0.9
adjacency = adjacency(datExpr, corFnc = "bicor", type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType = "signed", verbose = 0)
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method = "average")



################### Iterate WGCNA parameters #########
colors <- moduleLabel  <- c()

for (minModuleSize in c(5,10,15,20)) {
  for (ds in c(1:4)) {
  
    cutHeight = 0.99999
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                                deepSplit = ds, pamRespectsDendro = T,pamStage = T,
                                minClusterSize = minModuleSize, cutHeight = cutHeight)
    
    
    merged <- mergeCloseModules(exprData = datExpr,colors = labels2colors(dynamicMods))
    
    
    colors = cbind(colors, labels2colors((merged$colors)))
    moduleLabel <- c(moduleLabel,paste("DS=",ds,
                                       " mms=\n",minModuleSize))
    
  
    plotDendroAndColors(geneTree,colors,
                        groupLabels=moduleLabel,
                        addGuide=FALSE,
                        dendroLabels=FALSE,
                        main="Dendrogram With Different Module Cutting Parameters",
                        cex.colorLabels=0.5)
    
  }}

dim = 15

c_ref = as.character(colors[,dim]) #module size 20, ds 3
for (i in 1:dim(colors)[[2]]){
  ci = as.character(colors[,i])
  c_new = matchLabels(ci, c_ref)
  colors[,i] = c_new
}

colors = cbind(colors[,dim], colors)
moduleLabel = c("selected module", moduleLabel)

plotDendroAndColors(geneTree,colors,groupLabels = moduleLabel,addGuide=T,dendroLabels=F,cex.colorLabels=0.3)




###Run WGCNA with selected parameters
minModuleSize = 20  
ds = 3
cutHeight = 0.99999
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                            deepSplit = ds, pamRespectsDendro = T,pamStage = T,
                            minClusterSize = minModuleSize, cutHeight = cutHeight)
dynamicColors = labels2colors(dynamicMods)

#pdf("./Fig1C_related/moduledetection.pdf", useDingbats = FALSE)
plotDendroAndColors(geneTree,colors = dynamicColors,
                    groupLabels="moduleLabel",
                    addGuide=FALSE,
                    dendroLabels=FALSE,
                    main="",
                    cex.colorLabels=0.5)
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors,softPower = softPower)
MEs = MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")


#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
# Merge modules whose module eigengenes are highly correlated
MEDissThres = 0.15 #merging had no effect but ran as routine analysis
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors = merge$colors
MEs = merge$newMEs

MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")


# Calculate module membership before filtration
modNames=substring(names(MEs),3)
geneModulecor=corAndPvalue(datExpr, MEs, use = "p")
geneModuleMembership = as.data.frame(geneModulecor$cor)
MMPvalue = as.data.frame(geneModulecor$p)
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

mirornainfo = data.frame(miRID = colnames(datExpr), moduleColor = moduleColors, geneModuleMembership, MMPvalue)


###find reliable microRNA members by defining module membership based on the score

modules=data.frame(sub("ME","", names(MEs)))
dim(modules)

dir.create("filtered_mirs")
setwd("filtered_mirs")

dir.create("filtered_mirs")

threshold = 0.60 #membership correlation score


for (i in 1:nrow(modules)){
  module = modules[i,]
  modulemm = paste0("MM",module)
  print(modulemm)
  mirornainfo_select = mirornainfo[which(mirornainfo$moduleColor == module), ]  
  for (j in c(3:6)){
    filtered_genes<-subset(mirornainfo_select, mirornainfo_select[,modulemm]>threshold) 
    write.csv(filtered_genes, paste0("./filtered_mirs/", module, "_", "filtered_mirs.csv"))
    
  }
}

setwd("filtered_mirs/")
list.files()

temp = list.files(pattern="*.csv") #import all csv files
myfiles = lapply(temp, read.csv, header = TRUE)
combined.df <- do.call(rbind , myfiles)
combined.df = na.omit(combined.df)
table(combined.df$moduleColor)
summary_filtered = data.frame(table(combined.df$moduleColor))


setwd("../")
write.csv(summary_filtered, "summary_filtered_mirs.csv")

###Calculate the eignvalue for each module considering the member microRNAs with high membership score
##and update eigenvalue for each module
dynamicColors_filtered = combined.df$moduleColor
moduleColors_filtered = combined.df$moduleColor

datExpr_filtered = datExpr[,match(combined.df$miRID, names(datExpr))]

MEList = moduleEigengenes(datExpr_filtered, colors = dynamicColors_filtered,softPower = softPower)
MEs = MEList$eigengenes


### Correlate module and pheontype traits. 

nGenes = ncol(datExpr_filtered)
nSamples = nrow(datExpr_filtered)
moduleTraitcor = cor(MEs, datatrait_numeric, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitcor, nSamples)
textMatrix= paste(signif(moduleTraitcor, 2), "\n(", 
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitcor)

#plot Fig 1C
#plot module-phenotype correlation in a heatmap
labeledHeatmap(Matrix= moduleTraitcor, 
               xLabels= names(filtered_datTraits), 
               yLabels= names(MEs), 
               ySymbols= names(MEs), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.9, 
               zlim= c(-1,1), 
               plotLegend = TRUE,
               keep.dendro=TRUE,
               cex.lab.x = 0.9,
               cex.lab.y = 0.9,
               
               main= paste("microRNA modules- memory"))


# final module membership calculation
modNames=substring(names(MEs),3)
geneModulecor=corAndPvalue(datExpr_filtered, MEs, use = "p")
geneModuleMembership = as.data.frame(geneModulecor$cor)
MMPvalue = as.data.frame(geneModulecor$p)
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

mirornainfo = data.frame(miRID = colnames(datExpr_filtered), moduleColor = moduleColors_filtered, geneModuleMembership, MMPvalue)


### find blue module member microRNAs
blue_module_mirnas = mirornainfo[which(mirornainfo$moduleColor == "blue"),]
blue_module_mirnas = data.frame(blue_module_mirnas[,"miRID"])


### find brown module member microRNAs
brown_module_mirnas = mirornainfo[which(mirornainfo$moduleColor == "brown"),]
brown_module_mirnas = data.frame(brown_module_mirnas[,"miRID"])

### find turquoise module member microRNAs
turquoise_module_mirnas = mirornainfo[which(mirornainfo$moduleColor == "turquoise"),]
turquoise_module_mirnas = data.frame(turquoise_module_mirnas[,"miRID"])


### find yellow module member microRNAs
yellow_module_mirnas = mirornainfo[which(mirornainfo$moduleColor == "yellow"),]
yellow_module_mirnas = data.frame(yellow_module_mirnas[,"miRID"])



###### pathway analysis #####
### blue module 

mirna.list = data.frame(blue_module_mirnas[,1])
mirtar=read.csv("../../data/miRTarBase_MTI_v.7.csv", header = TRUE) 

#filter the database file based on the  mirna list
mir=mirtar[which(mirtar$miRNA %in% mirna.list[,1]), ]
mir_strong_validated_targets=mir[(mir$Support.Type=="Functional MTI"),] ##filtering target genes based on strong experimental evidence

  
################################### analyzing KEGG terms ##########
blue_kegg_strong <- enrichKEGG(gene = as.character(mir_strong_validated_targets[,5]),
                                   organism     = 'hsa',
                                   pvalueCutoff = 0.05, 
                                   pAdjustMethod = "BH")
blue_kegg_strong= blue_kegg_strong[(grep("pathway", blue_kegg_strong$Description)),] #select pathway related terms
write.csv(as.data.frame(blue_kegg_strong), "blue_module_kegg.csv")
 
  

#####brown module
  
  mirna.list = data.frame(brown_module_mirnas[,1])
  mirtar=read.csv("../../data/miRTarBase_MTI_v.7.csv", header = TRUE) 

  #filter the database file based on the  mirna list 
  mir=mirtar[which(mirtar$miRNA %in% mirna.list[,1]), ]

  mir_strong_validated_targets=mir[(mir$Support.Type=="Functional MTI"),] ##filtering target genes based on strong experimental evidence

  ################################### analyzing KEGG terms ##########
brown_kegg_strong <- enrichKEGG(gene = as.character(mir_strong_validated_targets[,5]),
                                   organism     = 'hsa',
                                   pvalueCutoff = 0.05, 
                                   pAdjustMethod = "BH")
  brown_kegg_strong= brown_kegg_strong[(grep("pathway", brown_kegg_strong$Description)),] #select pathway related terms
  
  write.csv(as.data.frame(brown_kegg_strong), "brown_module_kegg.csv")

  
  
  
  
  ####### turquoise module
  mirna.list = data.frame(turquoise_module_mirnas[,1])
  mirtar=read.csv("../../data/miRTarBase_MTI_v.7.csv", header = TRUE) 
 
  #filter the database file based on the  mirna list 
  mir=mirtar[which(mirtar$miRNA %in% mirna.list[,1]), ]
  mir_strong_validated_targets=mir[(mir$Support.Type=="Functional MTI"),] ##filtering target genes based on strong experimental evidence
  

  ################################### analyzing KEGG terms ##########
  turquoise_kegg_strong <- enrichKEGG(gene = as.character(mir_strong_validated_targets[,5]),
                                   organism     = 'hsa',
                                   pvalueCutoff = 0.05, 
                                   pAdjustMethod = "BH")
  turquoise_kegg_strong= turquoise_kegg_strong[(grep("pathway", turquoise_kegg_strong$Description)),] #select pathway related terms
  
  write.csv(as.data.frame(turquoise_kegg_strong), "turquoise_module_kegg.csv")
  
  
  
  
  ####### yellow module
  mirna.list = data.frame(yellow_module_mirnas[,1])
  mirtar=read.csv("../../data/miRTarBase_MTI_v.7.csv", header = TRUE) 

  #filter the database file based on the  mirna list 
  mir=mirtar[which(mirtar$miRNA %in% mirna.list[,1]), ]
  
  mir_strong_validated_targets=mir[(mir$Support.Type=="Functional MTI"),] ##filtering target genes based on strong experimental evidence
  
  ################################### analyzing KEGG terms ##########
  yellow_kegg_strong <- enrichKEGG(gene = as.character(mir_strong_validated_targets[,5]),
                                   organism     = 'hsa',
                                   pvalueCutoff = 0.05, 
                                   pAdjustMethod = "BH")
  yellow_kegg_strong= yellow_kegg_strong[(grep("pathway", yellow_kegg_strong$Description)),] #select pathway related terms
  
  write.csv(as.data.frame(yellow_kegg_strong), "yellow_module_kegg.csv")
  
  
#### plot 23 common pathways among blue, brown and turquoise modules
  
  pathways_of_interest = data.frame(c(
    "AGE-RAGE signaling pathway in diabetic complications",
    "FoxO signaling pathway",
    "PI3K-Akt signaling pathway",
    "MAPK signaling pathway",
    "TNF signaling pathway",
    "Signaling pathways regulating pluripotency of stem cells",
    "Neurotrophin signaling pathway",
    "PD-L1 expression and PD-1 checkpoint pathway in cancer",
    "p53 signaling pathway",
    "Prolactin signaling pathway",
    "Toll-like receptor signaling pathway",
    "HIF-1 signaling pathway",
    "ErbB signaling pathway",
    "T cell receptor signaling pathway",
    "Rap1 signaling pathway",
    "Longevity regulating pathway",
    "TGF-beta signaling pathway",
    "Hippo signaling pathway",
    "Longevity regulating pathway - multiple species",
    "Chemokine signaling pathway",
    "JAK-STAT signaling pathway",
    "mTOR signaling pathway",
    "Insulin signaling pathway"
  ))
  
  # blue_kegg = read.csv("./blue_GO_analysis/blue_KEGG_pathway.csv", header = T, sep = ",")
  # blue_kegg$Cluster = "blue"
  
  pathways_in_blue = blue_kegg_strong[match(pathways_of_interest[,1], blue_kegg_strong$Description),]
  pathways_in_brown = brown_kegg_strong[match(pathways_of_interest[,1], brown_kegg_strong$Description),]
  pathways_in_turquoise = turquoise_kegg_strong[match(pathways_of_interest[,1], turquoise_kegg_strong$Description),]
  

  all_kegg=merge_result(list(MEblue=pathways_in_blue,
                             MEbrown = pathways_in_brown,
                             MEturquoise = pathways_in_turquoise
                          ))
  
  
## plot Fig 1D
 dotplot(object = all_kegg, x = ~Cluster,
                    color = "p.adjust", showCategory = 18, split = NULL,
                    font.size = 12, title = "", by = "count", includeAll = TRUE)
  
  

  


