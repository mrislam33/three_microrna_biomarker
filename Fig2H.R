####
load("./data/metaDataFigure2H.RData")

####################################################
### Auxiliary function ############################
###################################################


hypergeometric_test=function(list_a,  list_b, gs, initials){
  
  intersect_a_b = intersect(list_a, list_b)
  union_a_b = union(list_a, list_b) 
  gs = as.numeric(gs) 
  
  x = as.numeric(length(setdiff(list_a, list_b))) 
  y = as.numeric(length(setdiff(list_b, list_a))) 
  z = gs - length(union_a_b)   
  i = as.numeric(length(intersect_a_b))
  set_a = as.numeric(length(list_a)) 
  set_b = as.numeric(length(list_b))  
  cat("analysis is for:", initials)
  cat("intersection:", i)
  
  ###do Fisher test
  print(fisher.test(matrix(c(z, x, y, i), nrow = 2))) 
  
  Fold_enrichment = (i/set_a)/(set_b/z)  
  
  cat("fold enrichment: ", Fold_enrichment)
  
}

####################################################################
message("Hypergeometric overlap analysis")
####################################################################

hypergeometric_test(list_a = let_7b_5p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs, initials = "let-7b-5p")
hypergeometric_test(list_a = miR_130b_3p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs, initials = "miR_130b_3p")
hypergeometric_test(list_a = miR_146a_5p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs, initials = "miR_146a-5p")
hypergeometric_test(list_a = miR_148a_3p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs, initials = "miR_148a_3p")
hypergeometric_test(list_a = miR_181a_5p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs,initials = "miR_181a_5p")
hypergeometric_test(list_a = miR_192_5p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs, initials = "miR_192_5p")
hypergeometric_test(list_a = miR_30a_3p[,1], list_b = gwas_davies_memory_genes[,1], gs = gs, initials = "miR_30a_3p")

message ("miR-146a-5p (1.8), 
          miR-148a-5p (1.72) 
          and miR-181a-5p (1.65)
          showed higher overlap and passed the threshold as stated in Figure legend of the manuscript")
