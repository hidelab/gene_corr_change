rm(list=ls())
source("/home/md1wwxx/sharedmd1wwxx/src/gene.cor.long.R")
source("/home/md1wwxx/sharedmd1wwxx/src/pdiffcorr.R")
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

OSGWAS_genes<-c("GLT8D1", "GNL3", "ASTN2", "FILIP1", "SENP6", "KLHDC5", "PTHLH", "CHST11", "TP63", "FTO", "SUPT3H", "CDC5L")
OSGWAS_genes<-sub("KLHDC5", "KLHL42", OSGWAS_genes)
setdiff(OSGWAS_genes, rownames(x))

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

intersect(OSGWAS_genes, sigGenes)
sigEdge.GWAS<-subset(sigEdge2, subset=(Gene.A %in% OSGWAS_genes | Gene.B %in% OSGWAS_genes))
library(WriteXLS)
WriteXLS(sigEdge.GWAS, "OA2_D_vs_C_protein_corr_change_sigEdge_GWAS_p1e-04.xlsx")
intact<-as.data.frame(control)
degraded<-as.data.frame(treat)
par(mfrow=c(1,2))
plot(intact$GNL3, intact$CCAR1)
plot(degraded$GNL3, degraded$CCAR1)

c_net<-subset(allResults, subset=(FDR.c < 0.1))
c_net_GWAS_sigGenes<-subset(c_net, subset=(((Gene.A %in% OSGWAS_genes) & (Gene.B %in% sigGenes)) | ((Gene.B %in% OSGWAS_genes) & (Gene.A %in% sigGenes)) ))
nrow(c_net_GWAS_sigGenes)
x<-sort(table(c(c_net_GWAS_sigGenes$"Gene.A", c_net_GWAS_sigGenes$"Gene.B")), decreasing = TRUE)
class(x)
x1<-as.data.frame(x)
colnames(x1)<-c("Gene", "Number of significant correlation in control chondrocytes with genes involved in gene-gene correlation")
OSGWAS_genes_df<-data.frame("Gene"=OSGWAS_genes)
OSGWAS_genes_df<-merge(OSGWAS_genes_df, x1, all.x=TRUE, by.x="Gene", by.y="Gene")

rm(c_net)
t_net<-subset(allResults, subset=(FDR.t < 0.1))
nrow(t_net)
t_net_GWAS_sigGenes<-subset(t_net, subset=(((Gene.A %in% OSGWAS_genes) & (Gene.B %in% sigGenes)) | ((Gene.B %in% OSGWAS_genes) & (Gene.A %in% sigGenes)) ))
nrow(t_net_GWAS_sigGenes)
x<-sort(table(c(t_net_GWAS_sigGenes$"Gene.A", t_net_GWAS_sigGenes$"Gene.B")), decreasing = TRUE)
class(x)
x1<-as.data.frame(x)
colnames(x1)<-c("Gene", "Number of significant correlation in diseased chondrocytes with genes involved in gene-gene correlation")
OSGWAS_genes_df<-merge(OSGWAS_genes_df, x1, all.x=TRUE, by.x="Gene", by.y="Gene")
o<-order(OSGWAS_genes_df$Gene, decreasing = FALSE)
OSGWAS_genes_df<-OSGWAS_genes_df[o,]
write.table(OSGWAS_genes_df, file="OA2_GWAS_gene_sigdiffcor_protein_correlation_frequency.txt", sep="\t", row.names=FALSE, quote = FALSE)
