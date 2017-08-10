source("r2p.R")
library(reshape2)
gene.cor.long<-function(x, method="spearman") {
	corr<-cor(x, method=method)
	corr[lower.tri(corr, diag=TRUE)]<-NA	
	corr<-melt(corr)
	colnames(corr)<-c("Gene.A", "Gene.B", "correlation.coef")
	corr<-corr[!is.na(corr$correlation.coef),]
	n<-nrow(x)
	corr$p.value<-r2p(corr$correlation.coef, n)
	return(corr)
}


