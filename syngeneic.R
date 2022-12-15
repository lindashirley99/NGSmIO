library(scales)
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
options(java.parameters = "-Xmx8g")
library("xlsx")
library("clusterProfiler")
library("DOSE")
library(org.Mm.eg.db)


expr = read.csv("baseline_IOpanel_expr_Kallisto_TPM.csv", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
genes=row.names(expr)
symbols=sub("\\(.+$","",genes)
ids=sub("^.+\\(","",genes)
ids=sub("\\)$","",ids)
row.names(expr)=symbols
expr = log(expr+1,2)

anno = read.xlsx("Sample anno.xlsx",header = T,sheetName = "anno1",row.names = 1)
anno = anno[colnames(expr),]
#anno[,2] = as.character(anno[,2])
anno = anno[order(anno[,1]),]
samples=row.names(anno)
expr=expr[,samples]

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = symbols , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
expr_h = merge(genesV2,expr,by.x=1,by.y=0)
expr_h = expr_h[!duplicated(expr_h[,2]),]
row.names(expr_h) = expr_h[,2]
expr_h[,1:2] = NULL


expr_filtered = expr[rowSums(expr) != 0,]
#expr_filtered = expr_filtered[rowSums(expr_filtered) != -2*ncol(expr_filtered),]
#expr_filtered = expr_filtered[,apply(expr_filtered, 1, var, na.rm=TRUE) != 0]
pca <- prcomp(as.matrix(t(expr_filtered)), center = TRUE, scale. = TRUE)
eigs <- pca$sdev^2
proVar = eigs/sum(eigs)*100
scores = as.data.frame(pca$x)
scores$Model=anno$Model
scores$Size=anno$`Tumor volume`
scores$Type=anno$`Cancer type`

pdf("PCA.pdf",10,8)
p <- ggplot(data = scores, aes(x = PC1, y = PC2,group=Model))+ geom_line()+
  geom_point(aes(x = PC1, y = PC2,shape= Size), color="grey",size=5,stat = "identity") 
p<- p + labs(x = paste("PC1(",format(proVar[1], digits=2, nsmall=2),"%)",sep = ""), y = paste("PC2(",format(proVar[2], digits=2, nsmall=2),"%)",sep = ""))+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p + ggforce::geom_mark_ellipse(aes(color = Type, fill = Type), alpha = 0.05) 
dev.off()


data=as.data.frame(as.table(as.matrix(expr)))
colnames(data)=c("Gene","Sample","Expression")
data=merge(data,anno,by.x=2,by.y=0)
# require(tidyverse)
# data <- data %>% unite(Group, c("Model", "Tumor volume"),remove = FALSE)

pdf("DensityPlot.pdf",16,12)
p <- ggplot(data, aes(x=Expression,colour=Group)) + geom_density()
p + geom_vline(aes(xintercept=mean(data[,"Expression"])),color="black", linetype="dashed", size=1)
dev.off()


load("C:/Linda/Crownbio/Projects/IO/automated platform/ref_MuIO.RData")
expr_gs = gsva(as.matrix(expr),gsc_io,min.sz=1,max.sz=500,mx.diff=TRUE,method="ssgsea")
pdf("Heatmap-ssgsea.pdf",28,16)
pheatmap(expr_gs,cluster_rows=F, cluster_cols=T,treeheight_col=20, legend=T,annotation_col = anno[,c(1,2,5)],
         show_rownames = T,show_colnames = T, fontsize=8, fontsize_row=8,fontsize_col=8,cellheight =8,
         cellwidth=8,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


#pairwise correlation
expr.cor = cor(expr, method = c("spearman"))

library(corrplot)
pdf("correlation matrix.pdf",16,16)
corrplot(expr.cor, method = 'square',type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,tl.cex=0.7,
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()
pdf("Heatmap-correlation matrix.pdf",22,20)
pheatmap(expr.cor,cluster_rows=T, cluster_cols=T,treeheight_col=20, legend=T,annotation_col = anno[,c(1,2,5)],
         show_rownames = T,show_colnames = T, fontsize=6, fontsize_row=6,fontsize_col=6,cellheight =6,
         cellwidth=6,scale = "none",color = colorRampPalette(c("white", "firebrick3"))(100))
dev.off()


######################################
#RNAseq
######################################

counts_IO =read.csv("./with RNAseq readcounts/IOPanel_readcount.csv", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
# genes = row.names(counts_IO)
# counts_IO = sapply(counts_IO, as.numeric)
# row.names(counts_IO)=genes
counts_RNAseq = read.csv("./with RNAseq readcounts/RNAseq_readcount.csv", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
# genes = row.names(counts_RNAseq)
# counts_RNAseq = sapply(counts_RNAseq, as.numeric)
# row.names(counts_RNAseq)=genes
anno_RNAseq = read.table("./with RNAseq readcounts/Sample anno.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F)

genes = intersect(row.names(counts_IO),row.names(counts_RNAseq))
counts_IO = counts_IO[genes,anno_RNAseq[,1]]
counts_RNAseq = counts_RNAseq[genes,anno_RNAseq[,2]]
data=cbind(colMeans(counts_IO),colMeans(counts_RNAseq))
colnames(data)=c("Mean counts IO","Mean counts RNAseq")
write.xlsx2(data,"mean counts 209 samples 1095 genes.xlsx",row.names = T)


expr_RNAseq = read.csv("./with RNAseq/RNAseq_expr_TPM.csv", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
samples = anno_RNAseq[,2]
expr_RNAseq = expr_RNAseq[,samples]
genes=row.names(expr_RNAseq)
symbols=sub("\\(.+$","",genes)
ids=sub("^.+\\(","",genes)
ids=sub("\\)$","",ids)
expr_RNAseq = cbind(symbols,expr_RNAseq)
expr_RNAseq = expr_RNAseq[!duplicated(expr_RNAseq[,1]),]
row.names(expr_RNAseq)=expr_RNAseq[,1]
expr_RNAseq[,1] = NULL
expr_RNAseq = log(expr_RNAseq+1,2)

samples = anno_RNAseq[,1]
colnames(expr_RNAseq) = samples
expr_filtered = expr[,samples]

genes = row.names(expr_filtered)


##############################################################################
# Compare mIO panel and RNAseq data with FACS data
##############################################################################

facs= read.xlsx("./FACS data/Summary cells per mg.xlsx", header = T, sheetIndex = 1,check.names = F,stringsAsFactors = F)
facs=facs[facs$SampleIO!="NA",]
row.names(facs)=facs$SampleIO
facs=facs[,6:15]
cells = colnames(facs)
cells = cells[c(1,2,3,4,5,7,8,9)]
samples = row.names(facs)
  
anno_filtered = anno[samples,]
a=expr_filtered[,samples]
b=expr_RNAseq[,samples]


io_gs=read.xlsx2("mIO RNAseq panel gene list 1080.xlsx",sheetIndex = 2, header = T,check.names = F,stringsAsFactors = F)
names = unique(io_gs$Annotation)
names = names[names != "Internal_Reference_Genes"]
gsc_io = list()
for(i in names){
  data = io_gs[io_gs$Annotation==i,]
  gsc_io[[i]] = data[,1]
}

ssgsea_IO = gsva(as.matrix(expr_filtered[,samples]),gsc_io,min.sz=1,max.sz=500,mx.diff=TRUE,method="gsva")
ssgsea_RNAseq = gsva(as.matrix(expr_RNAseq[genes,samples]),gsc_io,min.sz=1,max.sz=500,mx.diff=TRUE,method="gsva")

gsva_IO = gsva(as.matrix(expr),gsc_io,min.sz=1,max.sz=500,mx.diff=TRUE,method="gsva")

library(reshape2)
selected = grep("macrophage",names(gsc_io),value = T,ignore.case = T)
gsva_IO_selected = t(gsva_IO[selected,,drop=F])
data=melt(as.matrix(gsva_IO_selected))
colnames(data)=c("Sample","Cell","Score")
data=merge(data,anno,by.x=1,by.y=0)

library("ggplot2")
pdf("Macrophage scores over syngeniec models_IOpanel.pdf",20,10)
p <- ggplot(data = data, aes(x = reorder(Model, Score, FUN = median), y=Score)) 
p <- p + geom_boxplot(aes(fill=Model),outlier.shape = NA)
p <- p + geom_jitter()
p <- p + facet_wrap( ~ Cell, scales="free")
p <- p + xlab("Model") + ylab("Signature Score")
p <- p + guides(fill=guide_legend(title="Model"))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")
dev.off()

library(reshape2)
selected = grep("monocyte",names(gsc_io),value = T,ignore.case = T)
gsva_IO_selected = t(gsva_IO[selected,,drop=F])
data=melt(as.matrix(gsva_IO_selected))
colnames(data)=c("Sample","Cell","Score")
data=merge(data,anno,by.x=1,by.y=0)

library("ggplot2")
pdf("Monocyte scores over syngeniec models_IOpanel.pdf",20,10)
p <- ggplot(data = data, aes(x = reorder(Model, Score, FUN = median), y=Score)) 
p <- p + geom_boxplot(aes(fill=Model),outlier.shape = NA)
p <- p + geom_jitter()
p <- p + facet_wrap( ~ Cell, scales="free")
p <- p + xlab("Model") + ylab("Signature Score")
p <- p + guides(fill=guide_legend(title="Model"))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")
dev.off()

library("ggplot2")
pdf("Neutrophil scores over all 23 models_RNAseq.pdf",14,8)
p <- ggplot(data = data, aes(x=Model, y=Score_RNAseq)) 
p <- p + geom_boxplot(aes(fill=Model),outlier.shape = NA)
p <- p + geom_jitter()
#p <- p + facet_wrap( ~ Cell, scales="free")
p <- p + xlab("Model") + ylab("Neutrophil Signature Score")
p <- p + guides(fill=guide_legend(title="Model"))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# selected = c("T-cells","T_cell_CD4+_surface_marker","T_cell_CD8+_activation",
#              "Treg_upregulated_receptor","B-cells","monocyte_Ly6G+Ly6C-_signature_C57",
#              "monocyte_Ly6G-Ly6C+_signature_C57","NK_signature",
#              "macrophage_surface_marker","M1_macrophages","M2_macrophages")

selected = c("CD45","T_cells","T_cell_CD4+_surface_marker","T_cell_CD8+_surface_marker",
             "Treg_upregulated_receptor","Macrophages","monocyte_marker","NK_signature")
ssgsea_IO_selected=t(ssgsea_IO[selected,])
colnames(ssgsea_IO_selected)=cells
ssgsea_RNAseq_selected=t(ssgsea_RNAseq[selected,])
colnames(ssgsea_RNAseq_selected)=cells

cors1 = c()
pvs1 = c()
cors2 = c()
pvs2 = c()
for(c in cells){
  x=facs[,c]
  y=ssgsea_IO_selected[,c]
  cor1=cor.test(x,y,method="spearman")$estimate
  cors1 = c(cors1,cor1)
  pv1=cor.test(x,y,method="spearman")$p.value
  pvs1 = c(pvs1,pv1)
  
  y=ssgsea_RNAseq_selected[,c]
  cor2=cor.test(x,y,method="spearman")$estimate
  cors2 = c(cors2,cor2)
  pv2=cor.test(x,y,method="spearman")$p.value
  pvs2 = c(pvs2,pv2)
}
table=data.frame(cells,cors1,pvs1,cors2,pvs2)
colnames(table)=c("Cell","Rho IO-FACS", "pv IO-FACS","Rho RNAseq-FACS", "pv RNAseq-FACS")

write.xlsx2(table,"Correlation on immune profiling.xlsx",row.names = F)

library(reshape2)
data1=melt(as.matrix(facs[,cells]))
data2=melt(as.matrix(ssgsea_IO_selected[,cells]))
data3=melt(as.matrix(ssgsea_RNAseq_selected[,cells]))
data=cbind(data1,data2[,3],data3[,3])
colnames(data)=c("Sample","Cell","CellNumber_FACS","Score_IOpanel","Score_RNAseq")
data=merge(data,anno_filtered[,c(1,5)],by.x=1,by.y=0)
library(cowplot)
pdf("Scatter plots-FACS IO RNAseq.pdf",10,14)
for(i in 1:nrow(table)){
  COR_text = paste(c("R =","p ="),signif(as.numeric(table[i,2:3],3),3),collapse=" ")
  cell=table[i,1]
  temp=data[data$Cell==cell,]
  temp$Model <- factor(temp$Model)
  p1 <- ggplot(temp, aes(x=CellNumber_FACS, y=Score_IOpanel)) + geom_point(aes(color = Model,shape=Model),size=4) +
    scale_shape_manual(values=1:nlevels(temp$Model)) +
    theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    labs(title = COR_text, x = "Cell number per mg by FACS", y = "Signature score by NGS IO panel" ) +
    geom_smooth(method='lm')+scale_x_log10()
  COR_text = paste(c("R =","p ="),signif(as.numeric(table[i,4:5],3),3),collapse=" ")
  p2 <- ggplot(temp, aes(x=CellNumber_FACS, y=Score_RNAseq)) + geom_point(aes(color = Model,shape=Model),size=4) +
    scale_shape_manual(values=1:nlevels(temp$Model)) +
    theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    labs(title = COR_text, x = "Cell number per mg by FACS", y = "Signature score by RNAseq") +
    geom_smooth(method='lm')+scale_x_log10()
  p <- ggarrange(p1, p2, ncol=1,nrow=2, common.legend = TRUE, legend="right")
  title <- cowplot::ggdraw() + cowplot::draw_label(cell, fontface='bold')
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.05, 1)))
}
dev.off()
 


#FACS percentage data
facs1= read.xlsx("./Summary percentage in live cells.xlsx", header = T, sheetIndex = 1,check.names = F,stringsAsFactors = F)
facs1=facs1[facs1$SampleIO!="NA",]
row.names(facs1)=facs1$SampleIO
facs1=facs1[,6:17]/100
cells1 = colnames(facs1)
samples = row.names(facs1)
cells1 = cells1[c(3,4,8,11)]

#mouse immune reference
ref = read.csv("C:/Linda/Crownbio/Projects/Stroma/mouse immune/RNAseq_expr_log2FPKM_16 samples.csv",header = T,
               check.names = F)
ref_anno = read.table("C:/Linda/Crownbio/Projects/Stroma/mouse immune/Sample info 16 samples.txt", header = T, 
                      sep = "\t",check.names = F, row.names = 1, stringsAsFactors = F)  
#C57
ref_samples = row.names(ref_anno[ref_anno[,"Strain"]=="Balbc",])
ref_filtered = ref[,ref_samples ]
genes = sub("\\(.+$","",row.names(ref_filtered))
ref_filtered[,"Symbol"] = genes
ref_filtered = ref_filtered[!duplicated(ref_filtered[,"Symbol"]),]
row.names(ref_filtered) = ref_filtered[,"Symbol"]
ref_filtered[,"Symbol"] = NULL 
colnames(ref_filtered)=ref_anno[ref_samples ,"Cell_type"]
ref_filtered = ref_filtered[,order(colnames(ref_filtered))]

#convert FPKM to TPM
library(RNAontheBENCH)
ref_filtered_TPM = sapply(ref_filtered, function(x) (as.numeric(fpkm2tpm(2^x))))
row.names(ref_filtered_TPM) = row.names(ref_filtered)

#EPIC for mouse
library("EPIC")
mouseRef = list()
#mouseRef$refProfiles = ref_filtered_TPM
mouseRef$refProfiles = log2(ref_filtered_TPM+1)
temp = read.table("C:/Linda/Crownbio/Projects/Stroma/mouse immune/sigGenes_balbc_filtered.txt", header = T, 
                  sep = "\t",check.names = F, stringsAsFactors = F)  
temp=temp[,1]
#temp=temp[!temp[,1] %in% c("Hba-a2","Hba-a1","Hbb-bs","Hbb-bt"),]
#temp=temp[!temp %in% c("Ly6c2","Hba-a2","Ms4a6c","Hba-a2","Hba-a1","Hbb-bs","Hbb-bt")]

#add=c("Cd8a","Cd4","Cd3d","Cd3e","Cd8b1","Fcgr1","Cd79a","Klra1","Cd19")
#temp=unique(c(temp,add))
mouseRef$sigGenes = unique(temp)

library("EPIC")
score = EPIC(expr_filtered, reference = mouseRef, mRNA_cell = NULL, mRNA_cell_sub = NULL,
             sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
             constrainedSum = TRUE, rangeBasedOptim = FALSE)
epic_IO=score$cellFractions
write.table(epic_IO,"Cell fractions in IO.txt",sep = "\t")

genes = row.names(expr_filtered)
score = EPIC(expr_RNAseq, reference = mouseRef, mRNA_cell = NULL, mRNA_cell_sub = NULL,
             sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
             constrainedSum = TRUE, rangeBasedOptim = FALSE)
epic_RNAseq=score$cellFractions
write.table(epic_RNAseq,"Cell fractions in RNAseq.txt",sep = "\t")



epic_IO_selected=epic_IO[samples,c(7,8,3,6)]
colnames(epic_IO_selected)=cells1
epic_RNAseq_selected=epic_RNAseq[samples,c(7,8,3,6)]
colnames(epic_RNAseq_selected)=cells1

cors1 = c()
pvs1 = c()
cors2 = c()
pvs2 = c()
for(c in cells1){
  x=facs1[,c]
  y=epic_IO_selected[,c]
  cor1=cor.test(x,y,method="spearman")$estimate
  cors1 = c(cors1,cor1)
  pv1=cor.test(x,y,method="spearman")$p.value
  pvs1 = c(pvs1,pv1)
  
  y=epic_RNAseq_selected[,c]
  cor2=cor.test(x,y,method="spearman")$estimate
  cors2 = c(cors2,cor2)
  pv2=cor.test(x,y,method="spearman")$p.value
  pvs2 = c(pvs2,pv2)
}
table1=data.frame(cells1,cors1,pvs1,cors2,pvs2)
colnames(table1)=c("Cell","Rho IO-FACS", "pv IO-FACS","Rho RNAseq-FACS", "pv RNAseq-FACS")
write.xlsx2(table1,"Correlation on immune cell fractions_Balbc.xlsx",row.names = F)

library(reshape2)
data1=melt(as.matrix(facs1[,cells1]))
data2=melt(as.matrix(epic_IO_selected[,cells1]))
data3=melt(as.matrix(epic_RNAseq_selected[,cells1]))
data=cbind(data1,data2[,3],data3[,3])
colnames(data)=c("Sample","Cell","Fraction_FACS","Fraction_IOpanel","Fraction_RNAseq")
data=merge(data,anno_filtered[,c(1,5)],by.x=1,by.y=0)
library(cowplot)
pdf("Scatter plots-FACS IO RNAseq_percentage_Balbc.pdf",10,14)
for(i in 1:nrow(table1)){
  COR_text = paste(c("R =","p ="),signif(as.numeric(table1[i,2:3],3),3),collapse=" ")
  cell=table1[i,1]
  temp1=data[data$Cell==cell,]
  temp1$Model <- factor(temp1$Model)
  p1 <- ggplot(temp1, aes(x=Fraction_FACS, y=Fraction_IOpanel)) + geom_point(aes(color = Model,shape=Model),size=4) +
    scale_shape_manual(values=1:nlevels(temp1$Model)) +
    theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    labs(title = COR_text, x = "Cell Fractions in live cells by FACS", y = "Cell fractions in all cells by NGS IO panel" ) +
    geom_smooth(method='lm')
  COR_text = paste(c("R =","p ="),signif(as.numeric(table1[i,4:5],3),3),collapse=" ")
  p2 <- ggplot(temp1, aes(x=Fraction_FACS, y=Fraction_RNAseq)) + geom_point(aes(color = Model,shape=Model),size=4) +
    scale_shape_manual(values=1:nlevels(temp1$Model)) +
    theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    labs(title = COR_text, x = "Cell Fractions in live cells by FACS", y = "Cell fractions in all cells by RNAseq") +
    geom_smooth(method='lm')
  p <- ggarrange(p1, p2, ncol=1,nrow=2, common.legend = TRUE, legend="right")
  title <- cowplot::ggdraw() + cowplot::draw_label(cell, fontface='bold')
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.05, 1)))
}
dev.off()

