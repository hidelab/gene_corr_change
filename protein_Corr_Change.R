rm(list=ls())
source("gene.cor.long.R")
source("pdiffcorr.R")
options(stringsAsFactors = F)

# ==== Expression Background ====
setwd("/home/md1wwxx/sharedmd1wwxx/OA2/data/")
sampleInfo<-read.table(file="Proteomics_samples.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
protein.data<-read.table(file="Proteomics_NormAbund_15samples_22052017.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

dim(protein.data)
head(rownames(protein.data))
head(sampleInfo)
sampleInfo<-split(sampleInfo, sampleInfo$Tissue)

x<-protein.data
keep <- rowSums(x >0) >= 3
x <- x[keep,]
nrow(x)

protein.data<-as.matrix(protein.data)
protein.data<-t(protein.data)

control<-protein.data[rownames(protein.data) %in% sampleInfo$Intact$ProtID,]
treat<-protein.data[rownames(protein.data) %in% sampleInfo$Degraded$ProtID,]
rownames(control)
rownames(treat)

corr1<-gene.cor.long(control)
corr2<-gene.cor.long(treat)
all(corr1$Gene.A==corr2$Gene.A)
all(corr1$Gene.B==corr2$Gene.B)
colnames(corr1)[c(3,4)]<-c("r.c", "p.c")
colnames(corr2)[c(3,4)]<-c("r.t", "p.t")
res<-cbind(corr1, corr2[, c("r.t", "p.t")])
res$p.diffcor<-pdiffcor(
  r1=res$r.c,
  r2=res$r.t,
  n1=nrow(control),
  n2=nrow(treat)
)
res$FDR<-p.adjust(res$p.diffcor, method="fdr")
res$FDR.c<-p.adjust(res$p.c, method="fdr")
res$FDR.t<-p.adjust(res$p.t, method="fdr")

res$r.c<-signif(res$r.c, digits = 3)
res$r.t<-signif(res$r.t, digits = 3)
res$p.c<-signif(res$p.c, digits = 3)
res$p.t<-signif(res$p.t, digits = 3)
res$p.diffcor<-signif(res$p.diffcor, digits = 3)
res$FDR<-signif(res$FDR, digits = 3)
res$FDR.c<-signif(res$FDR.c, digits = 3)
res$FDR.t<-signif(res$FDR.t, digits = 3)
res$Gene.A<-as.character(res$Gene.A)
res$Gene.B<-as.character(res$Gene.B)

allResults<-res
sigEdge<-subset(allResults, subset=(p.diffcor < 10^-4))
nrow(sigEdge)
sigEdge$abs_corr_change<-abs(sigEdge$r.t-sigEdge$r.c)
sigEdge2<-sigEdge
x<-sort(table(c(sigEdge2$"Gene.A", sigEdge2$"Gene.B")), decreasing = TRUE)
class(x)
x1<-as.data.frame(x)
colnames(x1)<-c("Gene", "Frequency")
setwd("/shared/hidelab2/user/md1wwxx/OA2/output/res/")
write.table(x1, file="OA2_D_vs_C_protein_corr_change_sigEdge_p1e-04_Gene_frequency.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(sigEdge2, file="OA2_D_vs_C_protein_corr_change_sigEdge_p1e-04.txt", sep="\t", row.names=FALSE, quote = FALSE)
sigGenes<-x1$Gene
length(sigGenes)
