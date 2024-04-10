##eigenvalue calculation
##based on moduleEigengenes function of WGCNA R package

##data preparation of singular value decomposition analysis
expr = t(data)
maxVarExplained = 3
nPC = 1
nVarExplained = min (nPC, maxVarExplained)
colors = "blue"
modLevels = levels(factor(colors))
PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modLevels)))
averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modLevels)))
varExpl = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modLevels)))

modulename = modLevels[1]
restrict1 = as.character(colors) == as.character(modulename)
datModule1 = as.matrix(t(expr[, restrict1]))
n = dim(datModule1)[1]
p = dim(datModule1)[2]
datModule = t(scale(t(datModule1)))
datModule[is.na(datModule)] <- 0

#singular value decomposition
svd1 = svd(datModule, nu = min(n, p, nPC), nv = min (n, p, nPC))
varExpl[,1] = (svd1$d[1:min(n,p, nVarExplained)])^2/sum(svd1$d^2)
veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], t(datModule), use = "p")
eigenval = svd1$v[,1]
scaledExpr = scale(t(datModule1))
averExpr[,1] = rowMeans(scaledExpr, na.rm = TRUE)

corAve = cor(averExpr[,1], eigenval, use = "p")
corTest = cor.test(averExpr[,1], eigenval, use = "p")

if(!is.finite(corAve)){
  corAve = 0;
}
if (corAve<0){
  eigenval = -1* eigenval
}

eig.val = data.frame(eigenval)
rownames(eig.val) =  colnames(data)
