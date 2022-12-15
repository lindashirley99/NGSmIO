library("clusterProfiler")
library("DOSE")
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(RColorBrewer)
require("biomaRt")
library(ggplot2)
library(ggpubr)
library(pheatmap)
library("ggrepel")

###################################################################################
# Load and process NGSmIO data
###################################################################################
expr = read.csv("mIO_RSEM_TPM.csv", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
genes=row.names(expr)
symbols=sub("\\(.+$","",genes)
ids=sub("^.+\\(","",genes)
ids=sub("\\)$","",ids)
row.names(expr)=ids

anno = read.table("Sample anno.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
anno = anno[colnames(expr),]
anno = anno[order(row.names(anno)),]
samples=row.names(anno)
expr=expr[,samples]

io=read.table("mIO panel signatures.txt",header = T, sep="\t",check.names = F,stringsAsFactors = F)

genes=c()
for(i in 1:nrow(expr)){
  if(mean(as.numeric(expr[i,])) > 0){
    genes=c(genes,row.names(expr)[i])
  }
}
dat=expr[row.names(expr) %in% io[,1],]
genes=c()
for(i in 1:nrow(dat)){
  if(mean(as.numeric(dat[i,])) > 0){
    genes=c(genes,row.names(dat)[i])
  }
}

expr_g=cbind(expr,symbols)

x=expr_g[duplicated(expr_g[,"symbols"]),]
expr_g=expr_g[!duplicated(expr_g[,"symbols"]),]
for(i in 1:nrow(x)){
  symbol=x[i,"symbols"]
  for(j in 1:(ncol(x)-1)){
    expr_g[expr_g$symbols==symbol,j]=expr_g[expr_g$symbols==symbol,j]+x[i,j]
  }
}
row.names(expr_g)=expr_g[,"symbols"]
expr_g[,"symbols"]=NULL
#expr_g=log(expr_g+1,2)


expr_g_filtered=expr_g[row.names(expr_g) %in% io[,1],]
  
dat=expr_g_filtered[,samples]
dat=dat[rowSums(dat) != 0,]
pca <- prcomp(as.matrix(t(dat)), center = TRUE, scale. = TRUE)
eigs <- pca$sdev^2
proVar = eigs/sum(eigs)*100
scores = as.data.frame(pca$x)
scores$Model=anno[samples,"Model"]
scores$Tissue=anno[samples,"Tissue"]
scores$Treatment=anno[samples,"Condition"]
pdf("PCA_mIOpanel.pdf",10,8)
p <- ggplot(data = scores, aes(x = PC1, y = PC2))+ labs(x = paste("PC1(",format(proVar[1], digits=2, nsmall=2),"%)",sep = ""), y = paste("PC2(",format(proVar[2], digits=2, nsmall=2),"%)",sep = ""))
p <- p + geom_point(aes(x = PC1, y = PC2,shape= Model,size=Tissue,color=Treatment), stat = "identity") 
p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

###################################################################################
# Load and process Nanostring data
###################################################################################
nano360=read.table("NanoString360 signature genes.txt",header = T,sep = "\t", stringsAsFactors = F,check.names = F)
nano_norm=read.table("NanoString_NormalizedData.txt",header = T,sep ="\t",
                     row.names = 1,stringsAsFactors = F,check.names = F)
nano_norm=log(nano_norm,2)

anno_nano=read.table("NanoString_sample anno.txt",sep = "\t", header = T,row.names = 1,stringsAsFactors = F,check.names = F)
anno_nano = anno_nano[order(row.names(anno_nano)),]
samples=row.names(anno_nano)
nano_norm=nano_norm[,samples]

pca <- prcomp(as.matrix(t(nano_norm)), center = TRUE, scale. = TRUE)
eigs <- pca$sdev^2
proVar = eigs/sum(eigs)*100
scores = as.data.frame(pca$x)
scores$Model=anno_nano$Model
scores$Tissue=anno_nano$Tissue
scores$Treatment=anno_nano$Condition
pdf("PCA_nanostring.pdf",10,8)
p <- ggplot(data = scores, aes(x = PC1, y = PC2))+ labs(x = paste("PC1(",format(proVar[1], digits=2, nsmall=2),"%)",sep = ""), y = paste("PC2(",format(proVar[2], digits=2, nsmall=2),"%)",sep = ""))
p <- p + geom_point(aes(x = PC1, y = PC2,shape= Model,size=Tissue,color=Treatment), stat = "identity") 
p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

###################################################################################
# Comparison between NGSmIO and NanoString data
###################################################################################


common=nano360[nano360[,1] %in% io[,2],]
gsc=list()
names=unique(common[,2])
for(x in names){
  gsc[[x]]=common[common[,2]==x,1]
}

gsc_io=list()
names=unique(io[,2])
for(x in names){
  gsc_io[[x]]=io[io[,2]==x,1]
}

expr_g_filtered=expr_g[row.names(expr_g) %in% common[,1],]

samples=intersect(colnames(expr),colnames(nano_norm))
samples=samples[order(samples)]

ssgsea= gsva(as.matrix(expr_g_filtered),gsc,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
write.table(ssgsea,"ssGSEA common genes_crown.txt",sep = "\t",row.names = T)
pdf("Heatmap-ssgsea_crown.pdf",13,10)
pheatmap(ssgsea,cluster_rows=F, cluster_cols=T,treeheight_col=20, legend=T,annotation_col = anno[,"Group",drop=F],
         show_rownames = T,show_colnames = T, fontsize=10, fontsize_row=10,fontsize_col=10,cellheight =10,
         cellwidth=10,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()
ssgsea= gsva(as.matrix(expr_g_filtered[,samples]),gsc,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
pdf("Heatmap-ssgsea-woCluster_crown.pdf",13,10)
pheatmap(ssgsea,cluster_rows=F, cluster_cols=F,legend=T,annotation_col = anno[,"Group",drop=F],
         show_rownames = T,show_colnames = T, fontsize=10, fontsize_row=10,fontsize_col=10,cellheight =10,
         cellwidth=10,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

ssgsea= gsva(as.matrix(nano_norm),gsc,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
write.table(ssgsea,"ssGSEA common genes_nanostring.txt",sep = "\t",row.names = T)
pdf("Heatmap-ssgsea_nanostring.pdf",13,10)
pheatmap(ssgsea,cluster_rows=F, cluster_cols=T,treeheight_col=20, legend=T,annotation_col = anno_nano[,"Group",drop=F],
         show_rownames = T,show_colnames = T, fontsize=10, fontsize_row=10,fontsize_col=10,cellheight =10,
         cellwidth=10,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

ssgsea= gsva(as.matrix(nano_norm[,samples]),gsc,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
pdf("Heatmap-ssgsea-woCluster_nanostring.pdf",13,10)
pheatmap(ssgsea[,samples],cluster_rows=F, cluster_cols=F,legend=T,annotation_col = anno[,"Group",drop=F],
         show_rownames = T,show_colnames = T, fontsize=10, fontsize_row=10,fontsize_col=10,cellheight =10,
         cellwidth=10,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


io_specific=io[io[,3]!="IO_associated_functions",]
gsc=list()
names=unique(io_specific[,3])
for(x in names){
  gsc[[x]]=io_specific[io_specific[,3]==x,2]
}

#expr_g_filtered=expr_g[row.names(expr_g) %in% io_specific[,2],]

samples=intersect(colnames(expr),colnames(nano_norm))
samples=samples[order(samples)]

ssgsea= gsva(as.matrix(expr_g[,samples]),gsc,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
pdf("Heatmap-ssgsea-io-woCluster_crown.pdf",13,10)
pheatmap(ssgsea,cluster_rows=F, cluster_cols=F,legend=T,annotation_col = anno[,"Group",drop=F],
         show_rownames = T,show_colnames = T, fontsize=10, fontsize_row=10,fontsize_col=10,cellheight =10,
         cellwidth=10,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


ssgsea= gsva(as.matrix(nano_norm[,samples]),gsc,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
pdf("Heatmap-ssgsea-io-woCluster_nanostring.pdf",13,10)
pheatmap(ssgsea[,samples],cluster_rows=F, cluster_cols=F,legend=T,annotation_col = anno[,"Group",drop=F],
         show_rownames = T,show_colnames = T, fontsize=10, fontsize_row=10,fontsize_col=10,cellheight =10,
         cellwidth=10,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

#correlation
pvs=c()
cors=c()
pdf("pairwise correlation.pdf")
for(i in samples){
  x = nano_norm[common,i]
  y = expr_filtered[common,i]
  #y = log(y + 1, 2)
  data=as.data.frame(cbind(x,y))
  colnames(data)=c("x","y")
  pvalue=(cor.test(as.numeric(x),as.numeric(y),method="pearson"))$p.value
  pvs=c(pvs,pvalue)
  cor=as.vector((cor.test(as.numeric(x),as.numeric(y),method="pearson"))$estimate)
  cors=c(cors,cor)
  main=paste(i,"\nR = ",format(cor, format='e', digits=3), ", p-value = ", format(pvalue, format='e', digits=3))
  p <- ggplot(data, aes(x=x, y=y)) +geom_point(colour="black", shape=21, size = 3) +
    ggtitle(main)+xlab("Expression in Microarray (log2 intensity)") +ylab("Expression in NGSmIO (log2 TPM)")
  p <- p + geom_smooth(aes(color = NULL), method='lm', formula= y~x,color="red")+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  print(p)
}
dev.off()
data=cbind(samples,cors,pvs)
write.csv(data,"pairwise correlation table.csv")

#mouse immune reference
ref = read.csv("C:/Linda/Crownbio/Projects/Stroma/mouse immune/_log2FPKM_16 samples.csv",header = T,
               check.names = F)
ref_anno = read.table("C:/Linda/Crownbio/Projects/Stroma/mouse immune/Sample info 16 samples.txt", header = T, 
                      sep = "\t",check.names = F, row.names = 1, stringsAsFactors = F)  
samples = row.names(ref_anno[ref_anno[,"Strain"]=="Balbc",])
ref_filtered = ref[,samples]
genes = sub("\\(.+$","",row.names(ref_filtered))
ref_filtered[,"Symbol"] = genes
ref_filtered = ref_filtered[!duplicated(ref_filtered[,"Symbol"]),]
row.names(ref_filtered) = ref_filtered[,"Symbol"]
ref_filtered[,"Symbol"] = NULL 
colnames(ref_filtered)=ref_anno[samples,"Cell_type"]
ref_filtered = ref_filtered[,order(colnames(ref_filtered))]

#convert FPKM to TPM
library(RNAontheBENCH)
ref_filtered_TPM = sapply(ref_filtered, function(x) (as.numeric(fpkm2tpm(2^x))))
row.names(ref_filtered_TPM) = row.names(ref_filtered)

#EPIC for mouse
library("EPIC")
mouseRef = list()
mouseRef$refProfiles = 2^ref_filtered
temp = read.table("C:/Linda/Crownbio/Projects/Stroma/mouse immune/sigGenes_C57_filtered.txt", header = F, 
                  sep = "\t",check.names = F, stringsAsFactors = F)  

#temp=temp[!temp[,1] %in% c("Hba-a2","Hba-a1","Hbb-bs","Hbb-bt"),]
temp=temp[!temp[,1] %in% c("Ly6c2","Hba-a2","Ms4a6c"),]

add=c("Cd8a","Cd4","Cd3d","Cd3e","Cd8b1","Fcgr1","Cd79a","Klra1","Cd19")
temp=c(temp,add)
mouseRef$sigGenes = unique(temp)

load("mMCPcounter_signatures.RData")

temp = read.table("C:/Linda/Crownbio/Projects/Stroma/mouse immune/sigGenes_Balbc_filtered.txt", header = T, 
                  sep = "\t",check.names = F, stringsAsFactors = F)  
mouseRef$sigGenes = unique(temp[,1])

samples=intersect(colnames(expr),colnames(nano))
samples=samples[order(samples)]
expr_g_filtered=expr_g[row.names(expr_g) %in% io[,2],]
dat=expr_g[,samples]
#dat=dat[rowSums(dat) != (-2)*ncol(dat),]
#dat=dat[row.names(dat) %in% io[1:1090,2],]

dat=expr_g[row.names(expr_g) %in% io[1:1090,2],]

score = EPIC(dat, reference = mouseRef, mRNA_cell = NULL, mRNA_cell_sub = NULL,
             sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
             constrainedSum = TRUE, rangeBasedOptim = FALSE)
score$cellFractions
a=t(score$cellFractions)
pdf("Relative fraction of stromal cell types in 36 samples_crown.pdf",16,12)
par(mar=c(15,6,6,15), xpd=TRUE)
#colors=c("#FF9900", "#FFFFCC", "#CCFF99","#009900","#0066CC","#9933CC","#FF0099","#CC0000","#663300","#CCCCCC")
colors=c("#FF9900", "#FFFFCC", "#CCFF99","#009900","#0066CC","#9933CC","#FF0099","#CC0000","#CCCCCC")
barplot(a, col=colors, width=0.3,space = 0.1,las=3,xlab = NULL,ylab="Relative proportion",
        cex.axis = 1.3,cex.names = 1.3,cex =1.3)
legend("topright", inset=c(-0.2,0),fill=rev(colors), legend=rev(rownames(a)),cex = 1.3)
dev.off()
write.table(a,"Cell fractions in Crown 36 samples.txt",sep = "\t")


score = EPIC(2^nano_norm[,samples], reference = mouseRef, mRNA_cell = NULL, mRNA_cell_sub = NULL,
             sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
             constrainedSum = TRUE, rangeBasedOptim = FALSE)
score$cellFractions
a=t(score$cellFractions)
pdf("Relative fraction of stromal cell types in 12 samples_nanostring.pdf",16,12)
par(mar=c(15,6,6,15), xpd=TRUE)
#colors=c("#FF9900", "#FFFFCC", "#CCFF99","#009900","#0066CC","#9933CC","#FF0099","#CC0000","#663300","#CCCCCC")
colors=c("#FF9900", "#FFFFCC", "#CCFF99","#009900","#0066CC","#9933CC","#FF0099","#CC0000","#CCCCCC")
barplot(a, col=colors, width=0.3,space = 0.1,las=3,xlab = NULL,ylab="Relative proportion",
        cex.axis = 1.3,cex.names = 1.3,cex =1.3)
legend("topright", inset=c(-0.2,0),fill=rev(colors), legend=rev(rownames(a)),cex = 1.3)
dev.off()
write.table(a,"Cell fractions in 12 samples_nanostring.txt",sep = "\t")

genes=c("Cd4","Cd274","Gzmb","Prf1","Ido1","Nos2","Lag3","Vsir","Ptpn11")
pdf("Boxplots for genes of interest_mIO.pdf",12,8)
for(g in genes){
    a=expr_g[g,]
    a=log(a + 1, 2)
    data=data.frame(t(a))
    data=merge(data,anno,by.x=0,by.y=0)
    colnames(data)[2]="Expression"
    p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_boxplot()+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=.4,fill="black")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(g)
    print(p)
}
dev.off()

pdf("Barplots for genes of interest_nanostring.pdf",12,8)
for(g in genes){
  a=nano_norm[g,]
  data=data.frame(t(a))
  data=merge(data,anno,by.x=0,by.y=0)
  colnames(data)[2]="Expression"
  p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_bar(stat="identity",color="black")+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(g)
  print(p)
}
dev.off()

genes=c("Fcgr1","Tap1","Tap2","Parp9","Tapbp","H2-K1","Fcgr4")
pdf("Boxplots for IFN signaling_mIO.pdf",12,8)
for(g in genes){
  a=expr_g[g,]
  a=log(a + 1, 2)
  data=data.frame(t(a))
  data=merge(data,anno,by.x=0,by.y=0)
  colnames(data)[2]="Expression"
  p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_boxplot()+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=.4,fill="black")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(g)
  print(p)
}
dev.off()

pdf("Barplots for IFN signaling_nanostring.pdf",12,8)
for(g in genes){
  a=nano_norm[g,]
  data=data.frame(t(a))
  data=merge(data,anno,by.x=0,by.y=0)
  colnames(data)[2]="Expression"
  p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_bar(stat="identity",color="black")+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(g)
  print(p)
}
dev.off()

####################################################################################
# Compare with RNAseq data
####################################################################################

RNAseq = read.csv("RNAseq_log2FPKM.csv", header = T,row.names = 1,check.names = F,stringsAsFactors = F)

library(RNAontheBENCH)
RNAseq_TPM = sapply(RNAseq, function(x) (as.numeric(fpkm2tpm(2^x))))
row.names(RNAseq_TPM) = row.names(RNAseq)

RNAseq_TPM=log(RNAseq_TPM,2)

genes=row.names(RNAseq)
ids=sub("^.+\\(","",genes)
ids=sub("\\)$","",ids)
symbols=sub("\\(.+$","",genes)

row.names(RNAseq)=ids
RNAseq=cbind(symbols,RNAseq)
RNAseq=RNAseq[order(RNAseq[,2],decreasing = T),]
RNAseq = RNAseq[!duplicated(RNAseq[,1]),]
row.names(RNAseq)=RNAseq[,1]
RNAseq[,1]=NULL


RNAseq_anno = read.table("./RNA-seq/Sample anno.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
#RNAseq_anno=RNAseq_anno[RNAseq_anno$Group %in% c("Hepa16-tumor-G1","Hepa16-tumor-G4","MC38-tumor-G1","MC38-tumor-G4"),]
samples=row.names(RNAseq_anno)
RNAseq=RNAseq[,samples]

RNAseq_g=cbind(RNAseq_TPM,as.character(symbols))
x=RNAseq_g[duplicated(RNAseq_g[,"symbols"]),]
RNAseq_g=RNAseq_g[!duplicated(RNAseq_g[,"symbols"]),]
for(i in 1:nrow(x)){
  symbol=x[i,"symbols"]
  for(j in 1:(ncol(x)-1)){
    RNAseq_g[RNAseq_g$symbols==symbol,j]=RNAseq_g[RNAseq_g$symbols==symbol,j]+x[i,j]
  }
}
row.names(RNAseq_g)=RNAseq_g[,"symbols"]
RNAseq_g[,"symbols"]=NULL
RNAseq_g=log(RNAseq_g+1,2)

dat=log(RNAseq+1,2)
dat=dat[rowSums(dat) != 0,]
pca <- prcomp(as.matrix(t(dat)), center = TRUE, scale. = TRUE)
eigs <- pca$sdev^2
proVar = eigs/sum(eigs)*100
scores = as.data.frame(pca$x)
scores$Model=RNAseq_anno$Model
scores$Treatment=RNAseq_anno$Treatment
pdf("PCA_RNAseq.pdf",10,8)
p <- ggplot(data = scores, aes(x = PC1, y = PC2))+ labs(x = paste("PC1(",format(proVar[1], digits=2, nsmall=2),"%)",sep = ""), y = paste("PC2(",format(proVar[2], digits=2, nsmall=2),"%)",sep = ""))
p <- p + geom_point(aes(x = PC1, y = PC2,shape= Model,color=Treatment), size=4,stat = "identity") 
p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

nano360all=read.table("C:/Linda/Crownbio/Projects/IO/I0470-S2002/NanoString 770 genes.txt",header = T,sep = "\t",
                      stringsAsFactors = F,check.names = F)
common=intersect(io[,2], unique(nano360all[,2]))
common=intersect(common,row.names(expr_g))

RNAseq_filtered=RNAseq_g[common,c(1:5,16:25,36:40)]
#RNAseq_filtered[RNAseq_filtered < -2] = -2
colnames(RNAseq_filtered)=paste0("RNAseq_",colnames(RNAseq_filtered))
colnames(RNAseq_filtered)=sub("G4","G2",colnames(RNAseq_filtered))
expr_g_filtered= expr_g[common,]
expr_g_filtered= log(expr_g_filtered + 1, 2)
colnames(expr_g_filtered)=paste0("CrownIO_",colnames(expr_g_filtered))
nano_norm_filtered= nano_norm[common,]
colnames(nano_norm_filtered)=paste0("NanoString_",colnames(nano_norm_filtered))
mydata=cbind(nano_norm_filtered,expr_g_filtered,RNAseq_filtered)
corMatrix= cor(mydata, method = c("spearman"))
sample_anno=data.frame(colnames(mydata))
row.names(sample_anno)=colnames(mydata)
sample_anno$Platform=sub("_.+$","",colnames(mydata))
temp=sub("^.+_","",colnames(mydata))
sample_anno$Group=sub("-[0-9]+$","",temp)
sample_anno[,1]=NULL
pdf("Heatmap-correlation matrix.pdf",13,12)
pheatmap(corMatrix,cluster_rows=T, cluster_cols=T,treeheight_col=20, legend=T,annotation_col = sample_anno,
         annotation_row = sample_anno,show_rownames = T,show_colnames = T, fontsize=8, fontsize_row=8,fontsize_col=8,cellheight =8,
         cellwidth=8,scale = "none",color = colorRampPalette(c("white", "red"))(100))
dev.off()

write.table(corMatrix,file="correlation matrix-NanoString_CrownIO_RNAseq.txt",sep = "\t")
temp=grep("Tumor",colnames(corMatrix),value=T)
pdf("Heatmap-correlation matrix-tumor.pdf",13,12)
pheatmap(corMatrix[temp,temp],cluster_rows=T, cluster_cols=T,treeheight_col=20, legend=T,annotation_col = sample_anno[temp,],
         annotation_row = sample_anno[temp,],show_rownames = T,show_colnames = T, fontsize=8, fontsize_row=8,fontsize_col=8,cellheight =8,
         cellwidth=8,scale = "none",color = colorRampPalette(c("white", "red"))(100))
dev.off()

samples=intersect(colnames(expr),colnames(nano_norm))
samples=samples[order(samples)]
genes=c()
sample=c()
nanostring=c()
crown=c()
mean=c()
nano_norm_filtered=nano_norm[common,samples]
expr_g_filtered=expr_g[common,samples]
expr_g_filtered= log(expr_g_filtered + 1, 2)
RNAseq_filtered=RNAseq[common,c(1:5,16:25,36:40)]
#RNAseq_filtered[RNAseq_filtered < -2] = -2
colnames(RNAseq_filtered)=sub("G4","G2",colnames(RNAseq_filtered))
for(i in 1:length(samples)){
  x = nano_norm_filtered[,i]
  y = expr_g_filtered[,i]
  for(j in 1:length(common)){
    if(y[j] <= 0 && x[j] > quantile(nano_norm[,samples[1]],0.01)){
      genes=c(genes,common[j])
      sample=c(sample,samples[i])
      nanostring=c(nanostring,x[j])
      crown=c(crown,y[j])
      temp=sub("-[0-9].+$","",samples[i])
      m=mean(as.numeric(RNAseq_filtered[j,grep(temp,colnames(RNAseq_filtered),value = T)]))
      mean=c(mean,m)
    }
  }
}
data=data.frame(genes,sample,nanostring,crown,mean)
colnames(data)=c("Gene","Sample","NanoString","Crown","RNAseq")
write.table(data,"non-expressed genes in crown.txt",sep = "\t",quote = F,row.names = F)
library("plyr")
freq=count(data[,1])
quantile(nano_norm_filtered[,1], 0.05)
quantile(expr_g_filtered[,1], 0.05)

library(ggpubr)
my_comparisons <- list( c("Hepa16-Tumor-G1", "Hepa16-Tumor-G2"), c("MC38-Tumor-G1", "MC38-Tumor-G2"))
# genes=c("Cd4","Cd274","Gzmb","Prf1","Ido1","Nos2","Lag3","Vsir","Ptpn11")
# pdf("Comparison of immune markers in 3 platforms.pdf",16,12)
genes=c("Ifng","Fcgr1","Tap1","Tap2","Parp9","Tapbp","H2-K1","Fcgr4")
pdf("Comparison of IFN signaling genes in 3 platforms.pdf",16,12)
for(g in genes){
  #nanostring
  a=nano_norm[g,]
  data=data.frame(t(a))
  data=merge(data,anno,by.x=0,by.y=0)
  colnames(data)[2]="Expression"
  max=max(data[,2])
  p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_bar(stat="identity",color="black")+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p.nano <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.margin = unit(c(1,1,1,1), "cm"))
  
  #crown IO
  a=expr_g[g,]
  a=log(a + 1, 2)
  data=data.frame(t(a))
  data=merge(data,anno,by.x=0,by.y=0)
  colnames(data)[2]="Expression"
  if(max(data[,2]) > max){
    max=max(data[,2])
  }
  p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_boxplot()+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=.4,fill="black")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.margin = unit(c(1,1,1,1), "cm"))
  p.crown <- p + stat_compare_means(comparisons = my_comparisons,method = "t.test")+ stat_compare_means(label.y = max+1,method = "anova") 
  
  #RNAseq
  a=RNAseq_g[g,c(1:5,16:25,36:40)]
  a=log(a, 2)
  a[a < -2]=-2
  data=data.frame(t(a))
  data=merge(data,RNAseq_anno,by.x=0,by.y=0)
  data[,1]=sub("G4","G2",data[,1])
  data[,6]=sub("G4","G2",data[,6])
  colnames(data)[2]="Expression"
  if(max(data[,2]) > max){
    max=max(data[,2])
  }
  p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_boxplot()+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=.4,fill="black")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.margin = unit(c(1,1,1,1), "cm"))
  p.rna <- p + stat_compare_means(comparisons = my_comparisons,method = "t.test")+ stat_compare_means(label.y = max+1,method = "anova") 
  
  figure <- ggarrange(p.nano + ylim(c(-2,max+2)), p.crown + ylim(c(-2,max+2)), p.rna + ylim(c(-2,max+2)), 
                      labels = c("NanoString", "Crown IO","RNAseq"), ncol = 2, nrow = 2)
  print(annotate_figure(figure,top = text_grob(g, color = "black", face = "bold", size = 20)))
  
}
dev.off()

a=expr_g["Ifng",samples]
a=log(a + 1, 2)
data=data.frame(t(a))
data=merge(data,anno,by.x=0,by.y=0)
colnames(data)[2]="Expression"
pdf("IFNg expression_crown_12samples.pdf",10,8)
p <- ggplot(data,aes(x=Group, y=Expression,fill=Group)) + geom_bar(stat="identity",color="black")+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.margin = unit(c(1,1,1,1), "cm"))+ylim(c(-2,max+2))
dev.off()

#check G1G2 difference using RNAseq as gold standard
samples=intersect(colnames(expr),colnames(nano_norm))
samples=samples[order(samples)]
genes=common
nano_norm_filtered=nano_norm[common,c(5,6,11,12)]
expr_g_filtered=expr_g[common,c(13:18,31:36)]
expr_g_filtered= log(expr_g_filtered + 1, 2)
RNAseq_filtered=RNAseq_g[common,c(1:5,16:25,36:40)]
#RNAseq_filtered=log(RNAseq_filtered,2)
#RNAseq_filtered[RNAseq_filtered < -2] = -2
colnames(RNAseq_filtered)=sub("G4","G2",colnames(RNAseq_filtered))
data=list()
data[["Hepa 1-6 Microarray vs. RNAseq"]]= data.frame(matrix(NA,nrow = length(genes),ncol = 2))
row.names(data[["Hepa 1-6 Microarray vs. RNAseq"]])=genes
colnames(data[["Hepa 1-6 Microarray vs. RNAseq"]])=c("Microarray","RNAseq")
data[["Hepa 1-6 NGSmIO vs. RNAseq"]]= data.frame(matrix(NA,nrow = length(genes),ncol = 2))
row.names(data[["Hepa 1-6 NGSmIO vs. RNAseq"]])=genes
colnames(data[["Hepa 1-6 NGSmIO vs. RNAseq"]])=c("NGSmIO","RNAseq")
data[["MC38 Microarray vs. RNAseq"]]= data.frame(matrix(NA,nrow = length(genes),ncol = 2))
row.names(data[["MC38 Microarray vs. RNAseq"]])=genes
colnames(data[["MC38 Microarray vs. RNAseq"]])=c("Microarray","RNAseq")
data[["MC38 NGSmIO vs. RNAseq"]]= data.frame(matrix(NA,nrow = length(genes),ncol = 2))
row.names(data[["MC38 NGSmIO vs. RNAseq"]])=genes
colnames(data[["MC38 NGSmIO vs. RNAseq"]])=c("NGSmIO","RNAseq")
for(i in genes){
  d1 = nano_norm_filtered[i,2]-nano_norm_filtered[i,1]
  d2 = mean(as.numeric(expr_g_filtered[i,4:6]))-mean(as.numeric(expr_g_filtered[i,1:3]))
  d0 = mean(as.numeric(RNAseq_filtered[i,6:10]))-mean(as.numeric(RNAseq_filtered[i,1:5]))
  data[["Hepa 1-6 Microarray vs. RNAseq"]][i,]=c(d1,d0)
  data[["Hepa 1-6 NGSmIO vs. RNAseq"]][i,]=c(d2,d0)
  d1 = nano_norm_filtered[i,4]-nano_norm_filtered[i,3]
  d2 = mean(as.numeric(expr_g_filtered[i,10:12]))-mean(as.numeric(expr_g_filtered[i,7:9]))
  d0 = mean(as.numeric(RNAseq_filtered[i,16:20]))-mean(as.numeric(RNAseq_filtered[i,11:15]))
  data[["MC38 Microarray vs. RNAseq"]][i,]=c(d1,d0)
  data[["MC38 NGSmIO vs. RNAseq"]][i,]=c(d2,d0)
}

singleScatter <- function(data){
  x=data[,1]
  y=data[,2]
  pvalue=(cor.test(as.numeric(x),as.numeric(y),method="pearson"))$p.value
  cor=as.vector((cor.test(as.numeric(x),as.numeric(y),method="pearson"))$estimate)
  main=paste("\nR = ",format(cor, format='e', digits=3), ", p-value = ", format(pvalue, format='e', digits=3))
  p <- ggplot(data, aes(x=x, y=y)) +geom_point(colour="black", shape=21, size = 3) +
    ggtitle(main)+xlab(paste0("G2 over G1 Fold Change in ",colnames(data)[1])) +ylab(paste0("G2 over G1 Fold Change in ",colnames(data)[2]))
  p <- p + geom_smooth(aes(color = NULL), method='lm', formula= y~x,color="red")+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}
plots = lapply(data, singleScatter)

pdf("Performance-TumorG2vsG2difference1.pdf",10,11)
ggarrange(plotlist = plots,labels=c("Hepa 1-6 Microarray vs. RNAseq","Hepa 1-6 NGSmIO vs. RNAseq",
                                    "MC38 Microarray vs. RNAseq","MC38 NGSmIO vs. RNAseq"),ncol = 2,nrow = 2)
dev.off()

#look at low expression genes

common = intersect(row.names(expr_g_filtered),row.names(nano_norm_filtered))
common = intersect(common,row.names(RNAseq_filtered))

main="Hepa 1-6 Tumor G1"
list=c()
data=data.frame(matrix(NA,ncol = 3,nrow=0))
colnames(data)=c("NanoString","Crown","RNAseq")
for(i in common){
  x=mean(as.numeric(expr_filtered[i,13:15]))
  if(x < 1){
    line=c(nano_norm_filtered[i,5],x,mean(as.numeric(RNAseq_filtered[i,1:5])))
    data=rbind(data,line)
    row.names(data)[nrow(data)]=i
  }
}
colnames(data)=c("NanoString","Crown","RNAseq")
data=data[row.names(data)!="Pgpep1",]
list=c(list,row.names(data))

p <- ggplot(data, aes(x=NanoString, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
      ggtitle(paste0("\n",main))+xlab("Expression in Microarray (log2 intensity)") +ylab("Expression in RNAseq (log2 TPM)")
p1= p +theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#+geom_text_repel(data = data, aes(x=NanoString, y=RNAseq, label =row.names(data)))
p <- ggplot(data, aes(x=Crown, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
      ggtitle(paste0("\n",main))+xlab("Expression in NGSmIO (log2 TPM)") +ylab("Expression in RNAseq (log2 TPM)")
p2=p + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#+   geom_text_repel(data = data, aes(x=Crown, y=RNAseq, label =row.names(data)))
pdf(paste0(main,"_TPM1.pdf"),10,5)
ggarrange(p1,p2,ncol = 2,nrow = 1)
dev.off()

main="Hepa 1-6 Tumor G2"
data=data.frame(matrix(NA,ncol = 3,nrow=0))
colnames(data)=c("NanoString","Crown","RNAseq")
for(i in common){
  x=mean(as.numeric(expr_filtered[i,16:18]))
  if(x < 1){
    line=c(nano_norm_filtered[i,6],x,mean(as.numeric(RNAseq_filtered[i,6:10])))
    data=rbind(data,line)
    row.names(data)[nrow(data)]=i
  }
}
colnames(data)=c("NanoString","Crown","RNAseq")
data=data[row.names(data)!="Pgpep1",]
list=c(list,row.names(data))
p <- ggplot(data, aes(x=NanoString, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
  ggtitle(paste0("\n",main))+xlab("Expression in Microarray (log2 intensity)") +ylab("Expression in RNAseq (log2 TPM)")
p1= p +theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#+geom_text_repel(data = data, aes(x=NanoString, y=RNAseq, label =row.names(data)))
p <- ggplot(data, aes(x=Crown, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
  ggtitle(paste0("\n",main))+xlab("Expression in NGSmIO (log2 TPM)") +ylab("Expression in RNAseq (log2 TPM)")
p2=p + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#+   geom_text_repel(data = data, aes(x=Crown, y=RNAseq, label =row.names(data)))
pdf(paste0(main,"_TPM1.pdf"),10,5)
ggarrange(p1,p2,ncol = 2,nrow = 1)
dev.off()

main="MC38 Tumor G1"
data=data.frame(matrix(NA,ncol = 3,nrow=0))
colnames(data)=c("NanoString","Crown","RNAseq")
for(i in common){
  x=mean(as.numeric(expr_filtered[i,31:33]))
  if(x < 1){
    line=c(nano_norm_filtered[i,11],x,mean(as.numeric(RNAseq_filtered[i,11:15])))
    data=rbind(data,line)
    row.names(data)[nrow(data)]=i
  }
}
colnames(data)=c("NanoString","Crown","RNAseq")
data=data[row.names(data)!="Pgpep1",]
list=c(list,row.names(data))
p <- ggplot(data, aes(x=NanoString, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
  ggtitle(paste0("\n",main))+xlab("Expression in Microarray (log2 intensity)") +ylab("Expression in RNAseq (log2 TPM)")
p1= p +theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  #     geom_text_repel(data = data, aes(x=NanoString, y=RNAseq, label =row.names(data)))
p <- ggplot(data, aes(x=Crown, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
     ggtitle(paste0("\n",main))+xlab("Expression in NGSmIO (log2 TPM)") +ylab("Expression in RNAseq (log2 TPM)")
p2=p + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
   #    geom_text_repel(data = data, aes(x=Crown, y=RNAseq, label =row.names(data)))
pdf(paste0(main,"_TPM1.pdf"),10,5)
ggarrange(p1,p2,ncol = 2,nrow = 1)
dev.off()

main="MC38 Tumor G2"
data=data.frame(matrix(NA,ncol = 3,nrow=0))
colnames(data)=c("NanoString","Crown","RNAseq")
for(i in common){
  x=mean(as.numeric(expr_filtered[i,34:36]))
  if(x < 1){
    line=c(nano_norm_filtered[i,12],x,mean(as.numeric(RNAseq_filtered[i,16:20])))
    data=rbind(data,line)
    row.names(data)[nrow(data)]=i
  }
}
colnames(data)=c("NanoString","Crown","RNAseq")
data=data[row.names(data)!="Pgpep1",]
list=c(list,row.names(data))
p <- ggplot(data, aes(x=NanoString, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
  ggtitle(paste0("\n",main))+xlab("Expression in Microarray (log2 intensity)") +ylab("Expression in RNAseq (log2 TPM)")
p1= p +theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#+geom_text_repel(data = data, aes(x=NanoString, y=RNAseq, label =row.names(data)))
p <- ggplot(data, aes(x=Crown, y=RNAseq)) +geom_point(colour="black", shape=21, size = 3) +xlim(-2,7)+ylim(-2,7)+
  ggtitle(paste0("\n",main))+xlab("Expression in NGSmIO (log2 TPM)") +ylab("Expression in RNAseq (log2 TPM)")
p2=p + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#+   geom_text_repel(data = data, aes(x=Crown, y=RNAseq, label =row.names(data)))
pdf(paste0(main,"_TPM1.pdf"),10,5)
ggarrange(p1,p2,ncol = 2,nrow = 1)
dev.off()
