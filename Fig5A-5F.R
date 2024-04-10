##prepare Figure 5A-5F

rm (list=ls())
library(WGCNA)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(plyr)
library(sva)
library(flashClust)
library(effsize)
library(meta)
library(dplyr)
library(metafor) 
options(scipen=999)
options(warn=-1)



memory_mirnas_human = data.frame(c("hsa-miR-148a-3p",
                                   "hsa-miR-146a-5p",
                                   "hsa-miR-181a-5p"))


########signature analysis in plasma, Fig 5A #############
  load("./data/datExprgse90828.RData")
  
  #"datPheno contains phenotypic information"
#"datExpr contains qPCR array data curated for seven microRNAs (see Fig 2G) from GSE90828"
#"memory_mirnas_human contains list of 3 microRNA signature"
  
  
  
  ### whole network connectivity and visualization ###
  
  A = adjacency(datExpr,type="signed") 
  k = as.numeric(apply(A,2,sum))-1 
  Z.k = scale(k)
  thresholdZ.k = -2.5 
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
  sampleTree = flashClust(as.dist(1-A), method = "average")
  datColors = data.frame(outlier = outlierColor)
  
  plotDendroAndColors(sampleTree,groupLabels=names(datColors),cex.dendroLabels = 1, cex.rowText = 0.5,
  colors=datColors,main="Sample Dendrogram and Trait Heatmap")
  
  
  ####outlier detection and remvoal ####
  
  remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
  datExpr=datExpr[,!remove.samples]
  datPheno=data.frame(datPheno[!remove.samples,])
  names(datPheno) = "condition"
  datPheno = na.omit(datPheno)
  
  ################################################################
  ####calculation of 3 microRNA signature eigenvalue####
  ################################################################
  my_mirnas_in_mci_control_plasma_data_gse90828 = datExpr
  my_mirnas_in_mci_control_plasma_data_gse90828[is.na(my_mirnas_in_mci_control_plasma_data_gse90828)] = 0
  
  control_samples = my_mirnas_in_mci_control_plasma_data_gse90828[,c(1:27)]
  rowMeans(control_samples)
  
  mir181a_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[1,], function(x) x-rowMeans(control_samples)[[1]]))
  mir146a_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[2,], function(x) x-rowMeans(control_samples)[[2]]))
  mir192_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[3,], function(x) x-rowMeans(control_samples)[[3]]))
  mir130b_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[4,], function(x) x-rowMeans(control_samples)[[4]]))
  mir30a_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[5,], function(x) x-rowMeans(control_samples)[[5]]))
  mir7b_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[6,], function(x) x-rowMeans(control_samples)[[6]]))
  mir148a_normalized_to_control = data.frame(sapply(my_mirnas_in_mci_control_plasma_data_gse90828[7,], function(x) x-rowMeans(control_samples)[[7]]))
  
  
  mirdata = data.frame(rbind(t(mir181a_normalized_to_control), t(mir146a_normalized_to_control), t(mir192_normalized_to_control), 
                             t(mir130b_normalized_to_control), t(mir30a_normalized_to_control), t(mir7b_normalized_to_control),
                             t(mir148a_normalized_to_control)))
  
  rownames(mirdata) = rownames(my_mirnas_in_mci_control_plasma_data_gse90828)
  data = mirdata [match(memory_mirnas_human[,1], rownames(mirdata)), ]
  
  try(source("./scripts/eigenvalue_calculation.R"))
  
    eig.val$diagnosis = datPheno$condition
    
    plot_data = eig.val
    control_samples = plot_data[which(plot_data$diagnosis == "control"),]
    mean_control = mean(control_samples$eigenval)
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenval, function(x) x-mean_control)
    
#find outlying samples
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    ###Fig 5A
    ggboxplot(plot_data,
                                x = "diagnosis",
                                y = "eigenval",
                                add = "dotplot", 
                                fill = "diagnosis",
                                title = "MCI plasma") + stat_compare_means()
    
    rm(eig.val, plot_data)
  
  
  
  
  ####################### Fig 5B-F ###################

load("./data/delcode_paperIslam.RData")
  
  logCPM=cpm(rawDatExpr, normalized.lib.sizes=TRUE, log=FALSE)
  for (i in 1:dim(logCPM)[1]) {
    logCPM[i, ]=sapply(logCPM[i, ], function(x) log2(1+x))
  }
  
  print("preparing files for human")
  age=as.numeric(datPheno$age)
  sex <- as.numeric(as.factor(datPheno[,"sex"]))-1
  condition=as.numeric(as.factor(datPheno[,"diagnosis"]))-1
 
  regvars <- as.data.frame(cbind(age, sex, condition))
  
  
  
  mod=model.matrix(~as.factor(condition), data=regvars)
  mod0=model.matrix(~1, data=regvars)
  svseq_plot <- svaseq(as.matrix(logCPM), mod = mod, mod0 = mod0)$sv
  regvars <- as.data.frame(cbind(condition,age, sex,  svseq_plot))
  
  ## Run the regression and make the age,sex and other covariates adjusted normalized matrix
  datExpr.reg <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(logCPM))
  rownames(datExpr.reg) <- rownames(logCPM)
  colnames(datExpr.reg) <- colnames(logCPM)
  coefmat <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(regvars)+1)
  
  for (i in 1:nrow(logCPM)) {
    lmmod1 <- lm(as.numeric(logCPM[i,])~condition+age+sex+V4+V5+V6+V7+V8+V9+V10+V11+V12,data=regvars)
    
    coef <- coef(lmmod1)
    datExpr.reg[i,] <- logCPM[i,] - coef["sex"]*regvars[,"sex"] - coef["age"]*regvars[,"age"] - coef["V4"]*regvars[,"V4"] - coef["V5"]*regvars[,"V5"] - coef["V6"]*regvars[,"V6"] -
      coef["V7"]*regvars[,"V7"] - coef["V8"]*regvars[,"V8"] - coef["V9"]*regvars[,"V9"] - coef["V10"]*regvars[,"V10"] - coef["V11"]*regvars[,"V11"] - coef["V12"]*regvars[,"V12"]
  }


  datExpr_reg <- datExpr.reg
  A = adjacency(datExpr_reg,type="signed") # this calculates the whole network connectivity
  k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
  Z.k = scale(k)
  thresholdZ.k = -2.5
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
  sampleTree = flashClust(as.dist(1-A), method = "average")
  datColors = data.frame(outlier = outlierColor)
  
  plotDendroAndColors(sampleTree,groupLabels=names(datColors),cex.dendroLabels = 0.5, cex.rowText = 0.3,
  colors=datColors,main="Sample Dendrogram and Trait Heatmap")
  
  
  remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
  datExpr_reg=datExpr_reg[,!remove.samples]
  datPheno=data.frame(datPheno[!remove.samples,])
  
  delcodeExpr = data.frame(datExpr_reg) #age,gender, and latent factors corrected
  colnames(delcodeExpr) = colnames(datExpr_reg)
  
  
  
  print("in delcode_AD_data")
  my_mirnas_in_delcode_AD_data=delcodeExpr[match(memory_mirnas_human[,1], rownames(delcodeExpr)),]
  rownames(my_mirnas_in_delcode_AD_data) = memory_mirnas_human[,1]
  
  data=data.frame(my_mirnas_in_delcode_AD_data)
  colnames(data) = colnames(my_mirnas_in_delcode_AD_data)
  
  
  
  try(source("./scripts/eigenvalue_calculation.R"))
    
    eig.val$diagnosis = datPheno$diagnosis
    DelcodeData = eig.val
    DelcodeData$factor_diagnosis = as.factor(DelcodeData$diagnosis)
    names(DelcodeData) = c("eigenvalue", "diagnosis", "factor_diagnosis")
    DelcodeData$sex = datPheno$sex
    
    control=DelcodeData[which(DelcodeData$diagnosis == "control"),]
    mean_control = mean(control$eigenval)
    DelcodeData$normalized_eigenvalue = sapply(DelcodeData$eigenvalue, function(x) x-mean_control)
    
    
#find outlying samples
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    
    
    not_outliers = unlist(DelcodeData %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) )
    which(not_outliers == "FALSE")
    
    
    
    
    DelcodeData = DelcodeData[-c(86,137, 138), ] #outliers removed
    
    
    ###Fig 5B
   ggboxplot(DelcodeData,
                                  x = "diagnosis",
                                  y = "normalized_eigenvalue",
                                  add = "dotplot",
                                  #fill = "diagnosis",
                                  short.panel.labs = FALSE,
                                  outlier.shape = 19,
                                  ylab = "eigen-expression",
                                  title = "MCI blood") +  stat_compare_means()
    


    ###take only MCI samples
    MCI = DelcodeData[which(DelcodeData$diagnosis == "MCI"),]
    df_mci = data.frame(cbind(MCI$eigenvalue, MCI$eigenvalue))##df_mci contains the eigenvalue expression for MCI patients
    rownames(df_mci) = rownames(MCI)
    df <- scale(df_mci)

    
    ##Fig 5C
    #determine number of cluster

    # Elbow method
    fig_elbow_mci = fviz_nbclust(df, kmeans, method = "wss") +
      geom_vline(xintercept = 2, linetype = 2)+
      labs(subtitle = "Elbow method")
      print(fig_elbow_mci)
    # Dissimilarity matrix
    d <- dist(df, method = "euclidean")
    hc5 <- hclust(d, method = "complete" )

    # Cut tree into 2 groups
    sub_grp <- cutree(hc5, k = 2)

    # Dissimilarity matrix
    d <- dist(df, method = "euclidean")
    hc5 <- hclust(d, method = "complete" )

    # Cut tree into 2 groups
    sub_grp <- cutree(hc5, k = 2)

    # Number of members in each cluster
    table(sub_grp)
    df_mci %>%
      mutate(cluster = sub_grp) %>%
      head

    cluster_member = df_mci %>% mutate(cluster = sub_grp)
    rownames(cluster_member) = rownames(df_mci)
    cluster_member
    cluster_member_1 = cluster_member[cluster_member$cluster=="1",]
    cluster_member_2 = cluster_member[cluster_member$cluster=="2",]

    ##MCI contains cognitive score and eigenvalue
    MCI_eigenexpr_cluster_1 = datPheno[match(rownames(cluster_member_1), rownames(datPheno)),]
    MCI_eigenexpr_cluster_1$eigenvalue = cluster_member_1$X1
    
    MCI_eigenexpr_cluster_1_male =  MCI_eigenexpr_cluster_1[which(MCI_eigenexpr_cluster_1$sex == "m"), ]
    MCI_eigenexpr_cluster_1_female =  MCI_eigenexpr_cluster_1[which(MCI_eigenexpr_cluster_1$sex == "f"), ]
    
    MCI_eigenexpr_cluster_2 = datPheno[match(rownames(cluster_member_2), rownames(datPheno)),]
    MCI_eigenexpr_cluster_2$eigenvalue = cluster_member_2$X1
    
    ##For Fig 5D
    print(cor.test(MCI_eigenexpr_cluster_1$eigenvalue, MCI_eigenexpr_cluster_1$memory_score))
    print(cor.test(MCI_eigenexpr_cluster_2$eigenvalue, MCI_eigenexpr_cluster_2$memory_score))

    
    print(cor.test(MCI_eigenexpr_cluster_1_male$eigenvalue, MCI_eigenexpr_cluster_1_male$memory_score)) #cor -0.345318,  t = -2.3271, df = 40, p-value = 0.02511
    print(cor.test(MCI_eigenexpr_cluster_1_female$eigenvalue, MCI_eigenexpr_cluster_1_female$memory_score)) #cor -0.345318,  t = -2.3271, df = 40, p-value = 0.02511
    
    
    
    ##preparing data for Fig 5E
    plot_mci_to_ad = datPheno[match(rownames(DelcodeData),rownames(datPheno)), ]
    plot_mci_to_ad$eigenvalue = DelcodeData$normalized_eigenvalue
    MCI_conversion_to_AD = plot_mci_to_ad[-c(which(plot_mci_to_ad$mcitoad == "NA")),]
    
    
    
    ##plot Fig 5E
   ggboxplot(MCI_conversion_to_AD, x = "mcitoad", y = "eigenvalue", add = "dotplot",
                               #fill = "mcitoad",
                               short.panel.labs = FALSE, outlier.shape = NA,
                               ylab = "eigen-expression", title = "MCI to AD") +  stat_compare_means()
    
   
    rm(plot_data, eig.val)
  
  
  
  
  
  
  ######################## CSF samples ANALYSIS ############
  load("./data/CSFmetaData.RData") #data kindly provided by Jain et al. 2019
  
  
  #######################################################################
  ##################       data preprocessing             ###############
  #######################################################################
  
  print("remove low expressed microRNAs")
  zfGenes = datExpr
  nSamples=ncol(zfGenes)
  minimum_reads <- 1
  filter <- apply(zfGenes, 1, function(x) length(x[x>minimum_reads])>=nSamples/4)
  filtered <- zfGenes[filter,]
  
  
  #######################################################################
  ##################       data normalization and ######################
  ############### adjustment of covariates #############################
  #######################################################################
  
  
  logCPM=cpm(filtered, normalized.lib.sizes=TRUE, log=FALSE)
  for (i in 1:dim(logCPM)[1]) {
    logCPM[i, ]=sapply(logCPM[i, ], function(x) log2(1+x))
  }
  
  age=as.numeric(datPheno$age)
  sex <- as.numeric(as.factor(datPheno[,"gender"]))-1
  condition=as.numeric(as.factor(datPheno[,"condition"]))-1
  mod=model.matrix(~as.factor(condition), data=datPheno) 
  mod0=model.matrix(~1, data=datPheno)
  
  n.sv = num.sv(logCPM, mod = mod, method = "leek")
  svseq <- svaseq(as.matrix(logCPM), mod = mod, mod0 = mod0)
  svseq_plot <- svaseq(as.matrix(logCPM), mod = mod, mod0 = mod0)$sv
  regvars <- as.data.frame(cbind(condition,age, sex, svseq_plot))
  datExpr.reg <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(logCPM))
  rownames(datExpr.reg) <- rownames(logCPM)
  colnames(datExpr.reg) <- colnames(logCPM)
  coefmat <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(regvars)+1)
  
  for (i in 1:nrow(logCPM)) {
    lmmod1 <- lm(as.numeric(logCPM[i,])~condition+age+sex+V4,data=regvars)
    coef <- coef(lmmod1)
    datExpr.reg[i,] <- logCPM[i,] - coef["sex"]*regvars[,"sex"] - coef["age"]*regvars[,"age"]  - coef["V4"]*regvars[,"V4"] 
  }
  datExpr_reg <- datExpr.reg 
  
  #######################################################################
  ########## detarmination of network connectivity  ########## 
  #######################################################################
  
  A = adjacency(datExpr_reg,type="signed") # this calculates the whole network connectivity
  k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
  Z.k = scale(k)
  thresholdZ.k = -2.5 # often -2.5
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
  sampleTree = flashClust(as.dist(1-A), method = "average")
  datColors = data.frame(outlier = outlierColor)
  
  #plotDendroAndColors(sampleTree,groupLabels=names(datColors),cex.dendroLabels = 1, cex.rowText = 0.5,
  #colors=datColors,main="Sample Dendrogram and Trait Heatmap")
  
  
  dataExpr = datExpr_reg  ##no outliers removed 
  
  #######################################################################
  ########## Eigenvalue calculation of three microRNA signature ########## 
  #######################################################################
  
  
  my_mirnas_in_data=dataExpr[match(memory_mirnas_human[,1], rownames(dataExpr)),]
  rownames(my_mirnas_in_data) = memory_mirnas_human[,1]
  data=data.frame(my_mirnas_in_data)
  
  
  try(source("./scripts/eigenvalue_calculation.R"))
  
    plot_data = data.frame(cbind(eig.val[,1], data.frame(datPheno$condition)))
    names(plot_data) = c("eigenval", "diagnosis")
    
    control=plot_data[which(plot_data$diagnosis == "control"),]
    mean_control = mean(control$eigenval)
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenval, function(x) x-mean_control)
    
#find outlying samples
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    

###Fig 5F
    ggboxplot(plot_data, x = "diagnosis", y ="normalized_eigenvalue",
                             ylab = "eigengene expression", 
                             fill = "diagnosis",
                             title = "MCI CSF", 
                             add = "dotplot"
    ) + stat_compare_means()
    
    
  rm(eig.val, plot_data)
  
  
  
  
  
  
