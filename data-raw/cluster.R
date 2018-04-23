# library(iCluster)
data(simu.datasets)

fit=iCluster(simu.datasets, k=3, lambda=c(0.2,0.2))
plotiCluster(fit=fit, label=rownames(simu.datasets[[2]]))
compute.pod(fit)

#Plot heat map
# library(iCluster)
data(gbm)
data(coord)
chr=chr[,1]
# fit=iCluster2(datasets=gbm, k=3, lambda=list(0.44,0.33,0.28))
for(k in 2:5){
  cat(paste("K=",k,sep=""),'\n')
  cv.fit[[k]]=tune.iCluster2(gbm, k, mc.cores=4)
}
##Reproducibility index (RI) plot
plotRI(cv.fit)

##Based on the RI plot, k=3 is the best solution
best.fit=cv.fit[[3]]$best.fit

##Try different color schemes
plotHeatmap(fit=best.fit,datasets=simu.datasets,
sparse=c(TRUE,TRUE),col.scheme=list(bluered(256), greenred(256)))



plotHeatmap(fit=fit, datasets=datasets, feature.order=c(FALSE,TRUE,TRUE),
sparse=c(FALSE,TRUE,TRUE),plot.chr=c(TRUE,FALSE,FALSE), chr=chr)

cna <- readr::read_tsv("/home/mtr/rfactory/paad_utsw_2015/data_CNA.txt")
gene_names <- cna$Hugo_Symbol
patient_id <- colnames(cna[,3:111])
gene_data <- t(as.matrix(cna[,3:111], ncol=length(gene_names)  ))
colnames(gene_data) <- gene_names
rownames(gene_data) <- patient_id
panc_data <- list(cn = gene_data)
cv.fit=alist()
for(k in 2:5){

  cat(paste("K=",k,sep=""),'\n')
  cv.fit[[k]]=tune.iCluster2(panc_data, k, mc.cores=4)
}

fit=iCluster(panc_data, k=3, lambda=c(0.2,0.2))
plotiCluster(fit=fit, label=rownames(simu.datasets[[2]]))
compute.pod(fit)
fit=iCluster(datasets=gbm, k=3, lambda=list(0.44,0.33,0.28))


cv.fit=alist()
for(k in 2:5){

  cat(paste("K=",k,sep=""),'\n')
  cv.fit[[k]]=tune.iCluster2(panc_data, k, mc.cores=4)
}

cv.fit=alist()
for(k in 2:5){
cat(paste("K=",k,sep=""),'\n')
cv.fit[[k]]=tune.iCluster2(simu.datasets, k, mc.cores=6)
}
#Reproducibility index (RI) plot
plotRI(cv.fit)
#Based on the RI plot, k=3 is the best solution
best.fit=cv.fit[[3]]$best.fit
#Try different color schemes
plotHeatmap(fit=best.fit,datasets=simu.datasets,
sparse=c(TRUE,TRUE),col.scheme=list(bluered(256), greenred(256)))


data(gbm)
# setting the penalty parameter lambda=0 returns non-sparse fit
fit=iCluster2(datasets=gbm, k=3, lambda=list(0.44,0.33,0.28))
plotiCluster(fit=fit, label=rownames(gbm[[1]]))
