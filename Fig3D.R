library(gplots)

#################################################################
############################ load mouse data and run analysis #############
#################################################################


print("analyzing in this new dataset")

rm(list = ls())

###import file
file=read.csv("./data/mirna_celltype_brain_stem_Hoye_et_al_2017.csv", row.names = 1)

###data retrieved from Hoye et al. 2017,          ##https://pubmed.ncbi.nlm.nih.gov/28416596/,
        ##please note that this dataset was generated in mouse brain stem

heatmap=as.matrix(file)

##analysis
for (i in 1:dim(heatmap)[1]) {
  heatmap[i, ]=sapply(heatmap[i, ], function(x) (x-mean(heatmap[i, ]))/sd(heatmap[i, ]))
}

##inverse sign as higher Ct value represents low expression
for (i in 1:dim(heatmap)[1]) {
  heatmap[i, ]=sapply(heatmap[i, ], function(x) -1*x)
}


my_palette <- colorRampPalette(c("blue","white","orange", "red"))(n = 30)

##plot Fig 3D
heatmap.2(heatmap, density.info="none", Rowv=FALSE, Colv=FALSE,
          dendrogram="none",
          trace="none", 
          labRow=rownames(heatmap),
          labCol=colnames(heatmap), 
          col=my_palette, 
          cexCol = 0.8, 
          cexRow = 1,
          srtCol = 50)






