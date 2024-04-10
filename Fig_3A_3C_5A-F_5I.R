rm (list=ls())
library(WGCNA)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(plyr)
library(sva)
library(flashClust)
library(meta)
library(dplyr)
library(metafor) 
options(scipen=999)
options(warn=-1)


###three microRNAs of interest
memory_mirnas_mouse = data.frame(c("mmu-miR-148a-3p",
                                   "mmu-miR-146a-5p",
                                   "mmu-miR-181a-5p"))

memory_mirnas_human = data.frame(c("hsa-miR-148a-3p",
                                   "hsa-miR-146a-5p",
                                   "hsa-miR-181a-5p"))
##############################
####analysis of three microRNAs in various datasets
  

##### dataset from Swrup et al. 2019. GEO accession: GSE89983

  load("./data/gse89983_cortex_6m.RData") #load data
  my_mirnas_in_gse89983 =  datExpr[match(memory_mirnas_mouse[,1], rownames(datExpr)),]
  rownames(my_mirnas_in_gse89983) = memory_mirnas_mouse[,1]
  my_mirnas_in_gse89983[is.na(my_mirnas_in_gse89983)] = 0
  
  data = my_mirnas_in_gse89983
  
  
  try(source("./scripts/eigenvalue_calculation.R"))
  
  plot_data = data.frame(cbind(eig.val[,1], datTrait))
    
    
    
    control = plot_data[which(plot_data$Wt.Tg == "Wt"),]
    mean_control = mean(control$eig.val...1.)
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eig.val...1., function(x) x-mean_control)
    
#outlier detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(Wt.Tg) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    #plot
    ggboxplot(plot_data, x = "Wt.Tg", y = "normalized_eigenvalue",
                          palette  = c("grey", "red"), 
                          xlab = "genotype", ylab ="eigengene expression", 
                          outlier.shape = NA, title = "FTLD 6 months",
                          add = "dotplot",
                          fill = "Wt.Tg") + stat_compare_means()
    
  
  rm(eig.val, plot_data)
  

################################################################################
  
##### dataset from longitudinal mouse experiment, from this study
  load("./data/longitudinal_aging_mouse.RData")
  
  
 
  my_mirnas_in_longi_aging=longi_aging[match(memory_mirnas_mouse[,1], rownames(longi_aging)),]
  rownames(my_mirnas_in_longi_aging) = memory_mirnas_mouse[,1]
  my_mirnas_in_longi_aging[is.na(my_mirnas_in_longi_aging)] = 0
  
  data = my_mirnas_in_longi_aging
  
  
  try(source("./scripts/eigenvalue_calculation.R"))

    eig.val$diagnosis = condition.longi_aging$diagnosis
    plot_data_longi_aging = eig.val
    plot_data_longi_aging = plot_data_longi_aging[-c(1:10),] #dont consider 12 months data as no phenotypic information not available
    control=plot_data_longi_aging[which(plot_data_longi_aging$diagnosis == "month13.5"),]
    mean_control = mean(control$eigenval)
    
    #eigevalue normalized to young mice
    plot_data_longi_aging$normalized_eigenvalue = sapply(plot_data_longi_aging$eigenval, function(x) x-mean_control)
    

#outlier detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data_longi_aging %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    #plot
ggboxplot(plot_data_longi_aging, x = "diagnosis", y = "normalized_eigenvalue",
                                fill = "diagnosis",
                                add="dotplot",
                                short.panel.labs = FALSE, outlier.shape = NA,
                                ylab = "eigen-expression", title = "longitudinal aging blood") + stat_compare_means()
    
rm(eig.val, plot_data)
  

##### dataset from aging brains, from this study

################ from CA3 region #############
  
  load ("./data/ca3_aging.RData")
  my_mirnas_in_ca3_aging=ca3_aging[match(memory_mirnas_mouse[,1], rownames(ca3_aging)),]
  rownames(my_mirnas_in_ca3_aging) = memory_mirnas_mouse[,1]
  my_mirnas_in_ca3_aging[is.na(my_mirnas_in_ca3_aging)] = 0
  
  data=my_mirnas_in_ca3_aging
  
  try(source("./scripts/eigenvalue_calculation.R"))

    plot_data = data.frame(cbind(eig.val[,1], condition[,1]))
    names(plot_data) = c("eigenvalue", "Treatment")
    plot_data$eigenvalue = as.numeric(plot_data$eigenvalue)
    
    
    control=plot_data[which(plot_data$Treatment == "juvenile"),]
    mean_control = mean(control$eigenvalue )
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenvalue, function(x) x-mean_control)
    plot_data_ca3 = plot_data
    
#outlier detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data_ca3 %>%
                            group_by(Treatment) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    ggboxplot(plot_data, x = "Treatment", y = "normalized_eigenvalue",
                                     fill = "Treatment", 
                                     palette = c("grey",  "red"),
                                     #facet.by = "month", 
                                     add="dotplot",
                                     short.panel.labs = FALSE, outlier.shape = NA,
                                     ylab = "eigen-expression", title = "mouse ca3 aging") + stat_compare_means()
    
    
rm(eig.val, plot_data)



##### dataset from aging brains, from this study

################ from CA1 region #############
  
  load("./data/ca1_aging.RData")
  my_mirnas_in_ca1_aging=ca1_aging[match(memory_mirnas_mouse[,1], rownames(ca1_aging)),]
  rownames(my_mirnas_in_ca1_aging) = memory_mirnas_mouse[,1]
  my_mirnas_in_ca1_aging[is.na(my_mirnas_in_ca1_aging)] = 0
  
  data=my_mirnas_in_ca1_aging
  
  
  try(source("./scripts/eigenvalue_calculation.R"))
 
    
    plot_data = data.frame(cbind(eig.val[,1], condition[,1]))
    names(plot_data) = c("eigenvalue", "Treatment")
    plot_data$eigenvalue = as.numeric(plot_data$eigenvalue)
    
    
    control=plot_data[which(plot_data$Treatment == "juvenile"),]
    mean_control = mean(control$eigenvalue )
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenvalue, function(x) x-mean_control)
    plot_data_ca1 = plot_data
    
    
#outlier detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data_ca1 %>%
                            group_by(Treatment) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    
    ######## plot ############
    
ggboxplot(plot_data_ca1, x = "Treatment", y = "normalized_eigenvalue",
                                     fill = "Treatment", 
                                     palette = c("grey",  "red"),
                                     #facet.by = "month", 
                                     add="dotplot",
                                     short.panel.labs = FALSE, outlier.shape = NA,
                                     ylab = "eigen-expression", title = "mouse CA1 aging") + stat_compare_means()
   
    
    
    
  rm(data,  condition, eig.val, plot_data)
  
  
  
  
##### dataset from aging brains, from this study

################ from DG region #############
  
  
  load("./data/dg_aging.RData")
  my_mirnas_in_dg_aging=dg_aging[match(memory_mirnas_mouse[,1], rownames(dg_aging)),]
  rownames(my_mirnas_in_dg_aging) = memory_mirnas_mouse[,1]
  my_mirnas_in_dg_aging[is.na(my_mirnas_in_dg_aging)] = 0
  
  data=my_mirnas_in_dg_aging
  
  
  try(source("./scripts/eigenvalue_calculation.R"))

    plot_data = data.frame(cbind(eig.val[,1], condition[,1]))
    names(plot_data) = c("eigenvalue", "Treatment")
    plot_data$eigenvalue = as.numeric(plot_data$eigenvalue)
    
    
    control=plot_data[which(plot_data$Treatment == "juvenile"),]
    mean_control = mean(control$eigenvalue )
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenvalue, function(x) x-mean_control)
    plot_data_dg = plot_data
    
#outlier detection
    
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data_dg %>%
                            group_by(Treatment) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    ########plot ############
    
   ggboxplot(plot_data_dg, x = "Treatment", y = "normalized_eigenvalue",
                                    fill = "Treatment", 
                                    palette = c("grey",  "red"),
                                    #facet.by = "month", 
                                    add="dotplot",
                                    short.panel.labs = FALSE, outlier.shape = NA,
                                    ylab = "eigen-expression", title = "mouse dg aging") + stat_compare_means()
    
 
  rm(data, condition, eig.val, plot_data)
  
  
  
##### dataset from aging brains, from this study

################ from ACC region #############

load("./data/acc_aging.RData")
  my_mirnas_in_acc_aging=acc_aging[match(memory_mirnas_mouse[,1], rownames(acc_aging)),]
  rownames(my_mirnas_in_acc_aging) = memory_mirnas_mouse[,1]
  my_mirnas_in_acc_aging[is.na(my_mirnas_in_acc_aging)] = 0
  
  data=my_mirnas_in_acc_aging
  
  
  try(source("./scripts/eigenvalue_calculation.R")) #calculate eigenvalue

    plot_data = data.frame(cbind(eig.val[,1], condition[,1]))
    names(plot_data) = c("eigenvalue", "Treatment")
    plot_data$eigenvalue = as.numeric(plot_data$eigenvalue)
    
    
    control=plot_data[which(plot_data$Treatment == "juvenile"),]
    mean_control = mean(control$eigenvalue )
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenvalue, function(x) x-mean_control)
    plot_data_acc = plot_data
    
#potential outliers detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data_acc %>%
                            group_by(Treatment) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
####################
    
########plot############
    
  ggboxplot(plot_data_acc, x = "Treatment", y = "normalized_eigenvalue",
                                      fill = "Treatment", 
                                      palette = c("grey",  "red"),
                                      add="dotplot",
                                      #facet.by = "month", 
                                      ylim = c(-0.2,1),
                                      short.panel.labs = FALSE, outlier.shape = NA,
                                      ylab = "eigen-expression", title = "mouse ACC aging") + stat_compare_means()
    
  
    rm(data,condition, eig.val, plot_data)
  
  
  
##### dataset from APP/PS1 brains, from this study

load("./data/apps1_brain.RData")
  
  my_mirnas_in_apps1_brain=apps1_brain[match(memory_mirnas_mouse[,1], rownames(apps1_brain)),]
  rownames(my_mirnas_in_apps1_brain) = memory_mirnas_mouse[,1]
  my_mirnas_in_apps1_brain[is.na(my_mirnas_in_apps1_brain)] = 0
  
  data=my_mirnas_in_apps1_brain
  
  try(source("./scripts/eigenvalue_calculation.R"))
 
    genotype= data.frame((c( rep("transgenic",6), rep("control",6),rep("transgenic",6), rep("control",5))))
    
    plot_data  = data.frame(cbind(eig.val[,1], condition[,1], age[,1], genotype[,1]))
    names(plot_data) = c("eigenvalue", "condition", "age", "genotype")
    plot_data$eigenvalue=as.numeric(plot_data$eigenvalue)
    
#potential outliers detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(genotype) %>%
                            group_map(~ out_z(.x$eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
   #plot 8 months data
    plot_data_8m = plot_data[which(plot_data$age=="8m"), ]
    
    control=plot_data_8m[which(plot_data_8m$genotype == "control"),]
    mean_control = mean(control$eigenvalue)
    
    plot_data_8m$normalized_eigenvalue = sapply(plot_data_8m$eigenvalue, function(x) x-mean_control)
    
    
ggboxplot(plot_data_8m, x = "genotype", y = "normalized_eigenvalue",
                                     fill = "genotype", 
                                     palette = c("grey",  "red"),
                                     add="dotplot",
                                     #facet.by = "month", 
                                     ylim = c(-0.2,1),
                                     short.panel.labs = FALSE, outlier.shape = NA,
                                     ylab = "eigen-expression", title = "APP/PS1 8 months") + stat_compare_means()
    
    
    #plot 4 months data
    plot_data_4m = plot_data[which(plot_data$age =="4m"), ]
    control=plot_data_4m[which(plot_data_4m$genotype == "control"),]
    mean_control = mean(control$eigenvalue)
    
    plot_data_4m$normalized_eigenvalue = sapply(plot_data_4m$eigenvalue, function(x) x-mean_control)
    
    
   #plot
ggboxplot(plot_data_4m, x = "genotype", y = "normalized_eigenvalue",
                                     fill = "genotype", 
                                     palette = c("grey",  "red"),
                                     add="dotplot",
                                     #facet.by = "month", 
                                     ylim = c(-0.2,1),
                                     short.panel.labs = FALSE, outlier.shape = NA,
                                     ylab = "eigen-expression", title = "APP/PS1 4 months") + stat_compare_means()

    
    
  rm(data, condition, eig.val, plot_data)
  
  
  
##### dataset from Psycourse, from this study
  
  
  load("./data/datMeta_psycourse_aging.RData")
  ###normalizing to library size
  zfGenes = datExpr
  library(edgeR)
  nSamples=ncol(zfGenes)
  minimum_reads <- 1
  
  filter <- apply(zfGenes, 1, function(x) length(x[x>minimum_reads])>=nSamples/4)
  filtered <- zfGenes[filter,]
  
  
  logCPM=cpm(filtered, normalized.lib.sizes=TRUE, log=FALSE)
  for (i in 1:dim(logCPM)[1]) {
    logCPM[i, ]=sapply(logCPM[i, ], function(x) log2(1+x))
  }
  print("analysis is succesfully done and results are stored as logCPM")
  
  print("preparing files for human")
  age=as.numeric(datPheno$age)
  sex <- as.numeric(as.factor(datPheno[,"gender"]))-1
  
  condition=as.numeric(as.factor(datPheno$age_group))-1
  regvars <- as.data.frame(cbind(condition,age, sex))
  

  
  # ##finding latent surrogate variables#####
  mod=model.matrix(~as.factor(condition), data=regvars) #interested only in condition
  mod0=model.matrix(~1, data=regvars)
  
  #-- svaseq
  svseq_plot <- svaseq(as.matrix(logCPM), mod = mod, mod0 = mod0)$sv
  regvars <- as.data.frame(cbind(condition,age, sex, svseq_plot))
  
  ## Run the regression and make the age,sex and other covariates adjusted normalized matrix
  datExpr.reg <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(logCPM))
  rownames(datExpr.reg) <- rownames(logCPM)
  colnames(datExpr.reg) <- colnames(logCPM)
  coefmat <- matrix(NA,nrow=nrow(logCPM),ncol=ncol(regvars)+1)
  
  for (i in 1:nrow(logCPM)) {
    lmmod1 <- lm(as.numeric(logCPM[i,])~condition+sex+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14,data=regvars)
    
    coef <- coef(lmmod1)
    datExpr.reg[i,] <- logCPM[i,] - coef["sex"]*regvars[,"sex"]  - coef["V4"]*regvars[,"V4"]- coef["V5"]*regvars[,"V5"] - coef["V6"]*regvars[,"V6"] -
      coef["V7"]*regvars[,"V7"] - coef["V8"]*regvars[,"V8"] - coef["V9"]*regvars[,"V9"] - coef["V10"]*regvars[,"V10"] - coef["V11"]*regvars[,"V11"] - coef["V12"]*regvars[,"V12"] - coef["V13"]*regvars[,"V13"] - coef["V14"]*regvars[,"V14"]  
  }
  
  datExpr_reg <- datExpr.reg 
  
  data = data.frame(datExpr_reg) #sex and other covariates corrected
  
  
  ####calculation of whole network connectivity##
  A = adjacency(data,type="signed") # this calculates the whole network connectivity
  k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
  Z.k = scale(k)
  thresholdZ.k = -2.5
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
  sampleTree = flashClust(as.dist(1-A), method = "average")
  datColors = data.frame(outlier = outlierColor)
  
  try(plotDendroAndColors(sampleTree,groupLabels=names(datColors),cex.dendroLabels = 1, cex.rowText = 0.5,
  colors=datColors,main="Sample Dendrogram and Trait Heatmap"))
  
  
  #####outlier detection based on whole network connectivity##
  
  remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
  data=data[,!remove.samples]
  datPheno=datPheno[!remove.samples,]
  t_data = data.frame(t(data))
  
  my_mirnas_in_data=data[match(memory_mirnas_human[,1], rownames(data)),]
  rownames(my_mirnas_in_data) = memory_mirnas_human[,1]
  my_mirnas_in_data[is.na(my_mirnas_in_data)] = 0
  
  data = my_mirnas_in_data
  
  ##eigenvalue calculation###
  try(source("./scripts/eigenvalue_calculation.R"))


    
    dataEigen = data.frame(cbind(eig.val[,1], data.frame(datPheno$age_group), data.frame(datPheno$age)))
    plot_data = dataEigen
    names(plot_data) = c("eigenvalue", "condition", "age")
    rownames(plot_data) = rownames(datPheno)
    
    
    control=plot_data[which(plot_data$condition == "control"),]
    mean_control = mean(control$eigenvalue)
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenvalue, function(x) x-mean_control)
    plot_data$Treatment = plot_data$condition
    
    
    plot_data = data.frame(rbind(plot_data[which(plot_data$Treatment == "control"), ],
                                 plot_data[which(plot_data$Treatment == "middle"), ],
                                 plot_data[which(plot_data$Treatment == "old"), ]))
    
    ###find outliers based on eigenvalue
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(Treatment) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    plot_data = plot_data[-89,] #outliers removed
    
    my_comparison = list(c("control", "middle"), c("control", "old"))
    
##plot
ggboxplot(plot_data, x = "Treatment", y = "normalized_eigenvalue",
                                     fill = "Treatment", add = "dotplot",
                                     palette = c("grey", "#E7B800", "red"),
                                     #facet.by = "month", 
                                     short.panel.labs = FALSE,
                                     outlier.shape = NA,
                                     ylab = "eigen-expression", title = "human aging, blood") + stat_compare_means(comparisons  = my_comparison)
   
    
    rm(eig.val, plot_data)
  
  
##### dataset from plasma, GEO accession: GSE90828
  ###################################signature analysis in plasma#############
  load("./data/datExprgse90828.RData")
  
  print("datPheno contains phenotypic information")
  print("datExpr contains qPCR array data curated for seven microRNAs (see Fig 2G) from GSE90828")
  print("memory_mirnas_human contains list of 3 microRNA signature")
  
  
  
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
  
  try(source("./scripts/eigenvalue_calculation.R")) ##eigenvalue calculation

    eig.val$diagnosis = datPheno$condition
    
    plot_data = eig.val
    control_samples = plot_data[which(plot_data$diagnosis == "control"),]
    mean_control = mean(control_samples$eigenval)
    
    plot_data$normalized_eigenvalue = sapply(plot_data$eigenval, function(x) x-mean_control)
    
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    
ggboxplot(plot_data,
                                x = "diagnosis",
                                y = "eigenval",
                                add = "dotplot", 
                                fill = "diagnosis",
                                title = "MCI plasma") + stat_compare_means()
    
rm(eig.val, plot_data)
  
  


#######dataset from CSF, kindly provided by Jain et al. 2019   ############

  load("./data/CSFmetaData.RData") #load data
  
  
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
    
#outlier detection
out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    not_outliers = unlist(plot_data %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
  

##plot
 ggboxplot(plot_data, x = "diagnosis", y ="normalized_eigenvalue",
                             ylab = "eigengene expression", 
                             fill = "diagnosis",
                             title = "MCI CSF", 
                             add = "dotplot"
    ) + stat_compare_means()
    
    
rm(eig.val, plot_data)
  
  
##### dataset from DELCODE, from this study
  load("./data/delcode_paperIslam.RData")
  
  logCPM=cpm(rawDatExpr, normalized.lib.sizes=TRUE, log=FALSE)
  for (i in 1:dim(logCPM)[1]) {
    logCPM[i, ]=sapply(logCPM[i, ], function(x) log2(1+x))
  }
  
  print("preparing files for human")
  age=as.numeric(datPheno$age)
  sex <- as.numeric(as.factor(datPheno[,"sex"]))-1
  condition=as.numeric(as.factor(datPheno[,"diagnosis"]))-1
 
  regvars <- as.data.frame(cbind(age, sex,  condition))
  
  
  
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
    
    
#outlier detection
    out_z <- function(x, cutoff = 2.5, na.rm = TRUE) {
      abs(x - mean(x, na.rm = na.rm)) <= cutoff * sd(x, na.rm = na.rm)
    }
    
    
    
    not_outliers = unlist(DelcodeData %>%
                            group_by(diagnosis) %>%
                            group_map(~ out_z(.x$normalized_eigenvalue)) ) 
    which(not_outliers == "FALSE")
    
    
    
    
    DelcodeData = DelcodeData[-c(86,137, 138), ] #outliers removed
    
    
##plot
ggboxplot(DelcodeData,
                                  x = "diagnosis", 
                                  y = "normalized_eigenvalue", 
                                  add = "dotplot",
                                  #fill = "diagnosis", 
                                  short.panel.labs = FALSE, 
                                  outlier.shape = 19,
                                  ylab = "eigen-expression",
                                  title = "MCI blood") +  stat_compare_means()
    

##############
    
    plot_mci_to_ad = datPheno[match(rownames(DelcodeData),rownames(datPheno)), ]
    plot_mci_to_ad$eigenvalue = DelcodeData$normalized_eigenvalue
    MCI_conversion_to_AD = plot_mci_to_ad[-c(which(plot_mci_to_ad$mcitoad == "NA")),]
    
    
    
    
   ggboxplot(MCI_conversion_to_AD, x = "mcitoad", y = "eigenvalue", add = "dotplot",
                               #fill = "mcitoad", 
                               short.panel.labs = FALSE, outlier.shape = NA,
                               ylab = "eigen-expression", title = "MCI to AD") +  stat_compare_means()
    
  
  rm(plot_data, eig.val)
  
  
  
  
##prepare the file in the format necessary to perform the meta-analysis
  
###import file with pre-calculated effect size (d),

#d indicates effect size
#lower for lower interval of effect size
##upper for upper interval of effect size
## n for total numbers nc = number of samples in control group
#ne = number of samples in treatment group
#sd = standard deviation,
#var = variance,
#p = unadjusted p value,
#adjp = p value after multiple adjustments,
#adjp_val = asterisks,
#species = species,
#studyname = name of the experimental study

  
meta_manual = read.csv("./data/meta_analysis_three_mirnas.csv") #prepare the file in the given format
  
##meta analysis
meta.adjusted <- metagen(d,
                           lower = lower,
                           upper = upper,
                           data = meta_manual,
                           studlab=paste0(meta_manual$studyname, 
                                          " ", "(", meta_manual$adjp_val, ")"),
                           comb.fixed = FALSE,
                           comb.random = TRUE, #random.model because the effect is not known
                           method.tau = "SJ", #Sidik-Jonkman estimator
                           TE.tau = TRUE,
                          prediction = FALSE,
                           sm = "SMD"
  )

#species wise separation
species_subgroup = update.meta(meta.adjusted,
                                 byvar = meta_manual$Species,
                                 comb.random = TRUE,
                                 comb.fixed = FALSE)
  
#plot
forest(species_subgroup, layout = "JAMA",
         JAMA.pval = TRUE,
         overall.hetstat = TRUE,
         test.overall.random =  TRUE)

  
  

  
  
  


