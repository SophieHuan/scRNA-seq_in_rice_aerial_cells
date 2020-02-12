##Identification of Differential expressed genes (DEGs) induced by protoplasting
#Input the bulk data
TPM_bulk <- read.csv("bulk_RNA_seq_TPM.csv")
rownames(TPM_bulk) <- TPM_bulk$X
TPM_bulk <- TPM_bulk[,-1]
#Chooing expressed genes
TPM_exp <- TPM_bulk[rowSums(TPM_bulk)>0,]
TPM_exp$gene_id <- rownames(TPM_exp)
#DEGsinduced by protoplasting are selected by linear regression
TPM_lm <- reshape2::melt(TPM_exp, id.vars="gene_id", value.name="exp") 
TPM_lm <- TPM_lm[order(TPM_lm$gene_id),]
TPM_lm$batch <- rep(as.factor(c(1,2,2,3,3,1,2,2,2,2,3,3)),dim(TPM_exp)[1])
TPM_lm$treatment <- rep(as.factor(c(rep(0,5),rep(1,7))),dim(TPM_exp)[1])
TPM_p_lm <- data.frame(gene_id=rownames(TPM_exp))
for (i in 1:41420){
  TPM_p_lm[i,2] <- summary(lm(exp~batch+treatment,TPM_lm[TPM_lm$gene_id==rownames(TPM_exp)[i],]))$coefficients[4,4]
}
colnames(TPM_p_lm)[2] <- "p_value"
#Calculation of FDR
TPM_p_lm$FDR <- qvalue::qvalue(TPM_p_lm$p_value)$qvalues
#Calculation of fold change
for (i in 1:41420){
  TPM_p_lm[i,4] <- (mean(as.numeric(TPM_exp[i,6:12]))+0.1)/(mean(as.numeric(TPM_exp[i,1:5]))+0.1)}
colnames(TPM_p_lm)[4] <- "FC"
TPM_p_lm$DEG <- "N"
TPM_p_lm$DEG[abs(log2(TPM_p_lm$FC))>=1 & -log10(TPM_p_lm$FDR)>2] <- "Y"
#Choosing DEGs induced by protoplasting
DEG_protoplast <- TPM_p_lm[TPM_p_lm$DEG=="Y",]$gene_id

##Input sc-RNA-seq matrix: a row representing a gene, a column representing a cell
raw_Ctrl <- as.data.frame(data.table::fread("scRNA_seq_rice_Ctrl.csv"))
raw_HS <- as.data.frame(data.table::fread("scRNA_seq_rice_HS.csv"))
raw_LN <- as.data.frame(data.table::fread("scRNA_seq_rice_LN.csv"))

##Identification of highly variable genes (HVGs)
#Take Ctrl sample for example
rownames(raw_Ctrl) <- raw_Ctrl$V1
raw_Ctrl <- raw_Ctrl[,-1]
colnames(raw_Ctrl) <- gsub("-1","-Ctrl",colnames(raw_Ctrl))
#Defination of expressed genes (at least five UMIs in > 20 cells)
Ctrl <- raw_Ctrl[rowSums(raw_Ctrl >=5)>20,]
#Calculation of mean and cv
Dat_Ctrl <- data.frame(Gid = rownames(Ctrl),Mean = apply(Ctrl,1,mean),CV = apply(Ctrl,1,sd)/apply(Ctrl,1,mean))
#Calculation of distance to median (DM)
Cal.DM <- function (df) {
	df.sort <- df[order(df[ ,'Mean']), ]
	DMs <- data.frame()
	for (i in 1:dim(df.sort)[1]) {
		s <- (i-50):(i+50)
		s <- s[s %in% 1:dim(df.sort)[1]]
		a <- df.sort$CV[i] - median(df.sort$CV[s])
		DMs <- rbind(DMs, data.frame(Gid=df.sort$Gid[i], DM=a))
	}
	return(DMs)
}
Dat_Ctrl <- merge(Dat_Ctrl, Cal.DM(Dat_Ctrl))
#Defination of highly variable genes (DM >0.2)
feature_Ctrl <- Dat_Ctrl$Gid[Dat_Ctrl$DM>.2]
#"feature_LN" and "feature_HS" are also obtained using the above method
#Collection of HVGs from three samples 
feature_union <- union(union(feature_Ctrl, feature_LN), feature_HS)
#Removal of DEGs induced by protoplasting
library(tidyverse)
feature_union <- feature_union[!(feature_union %in% DEG_protoplast)]

##Correct for batch effect
library(SingleCellExperiment)
#Take Ctrl sample for example
u <- as.matrix(raw_Ctrl[feature_union, ])
#Calculation of CPM
v <- t(log2(t(u)*1e6/rowSums(t(u))+1))
sc_Ctrl <- SingleCellExperiment(list(counts=u, logcounts=v))
#"sc_LN" and "sc_HS" are also obtained using the above method 
library("scran")
out <- fastMNN(sc_Ctrl, sc_LN, sc_HS)

##Dimension reduction and cell clustering
library(umap)
umap <- umap(out$corrected[,1:50], n_neighbors=5, metric="cosine", min_dist=0.01, random_state=21)
umap_layout <- as.data.frame(umap$layout)
umap_layout$Sample <- c(rep("Ctrl",2027),rep("LN",1905),rep("HS",783))
umap_layout$cluster <- dbscan::hdbscan(umap_layout[,1:2], minPts = 120)$cluster
#Recorder clusters according to cell numbers, and the order of cluster 1 and 4 are not changed
umap_layout[umap_layout$cluster==2,]$cluster <- "tmp"
umap_layout[umap_layout$cluster==3,]$cluster <- 2
umap_layout[umap_layout$cluster==5,]$cluster <- 3
umap_layout[umap_layout$cluster=="tmp",]$cluster <- 5

##Assignment of clusters to cell types
#Identification of cell-type-specific marker genes from a previous dataset
#Input data from Jiao et al., (DOI: 10.1038/ng.282)
dat <- Jiao_data[,grep("Seedling...Leaf|Seedling...Shoot", colnames(Jiao_data))]
rownames(dat) <- Jiao_data$Gid
dat[is.na(dat)] <- 0
#Identification of marker genes by expression levels in all cell types
dat$Dom <- apply(dat, 1, function(x){ temp <- sort(x, decreasing=T); return((temp[1]+0.1)/(temp[2]+0.1))})
dat$NO1 <- apply(dat[,-dim(dat)[2]], 1, function(x){ temp <- sort(x, decreasing=T); return(names(temp)[1])})
#Choosing marker genes
dat <- dat[rowSums(dat[,1:(dim(dat)[2]-2)]>1)>=1 & dat$Dom>3, ]
Jiao_marker <- dat[,13:14]
#Identification of cluster-specific expressed genes(CSEGs) from our data set
raw_all <- cbind(raw_Ctrl,raw_LN,raw_Salt)
t_all <- as.data.frame(t(raw_all))
#Calculation of the average expression in one cell type and the proportion of cells expressing the gene in the cluster
cseg <- function(x) {
	tmp1 <- data.frame(cell=rownames(t_all), V1=t_all[,x])
	tmp1 <- merge(tmp1, umap_layout[umap_layout$cluster>0,c('cell','cluster')], by='cell')
	tmp2 <- aggregate(tmp1$V1, by=list(tmp1$cluster), mean) %>% dplyr::arrange(desc(x))
	colnames(tmp2) <- c('cluster','mean')
	tmp3 <- tmp1[tmp1$V1>=2, ]
	if(dim(tmp3)[1]>2) {
		tmp4 <- aggregate(tmp3$V1, by=list(tmp3$cluster), length)
		colnames(tmp4) <- c('cluster','n1')
		tmp5 <- aggregate(tmp1$V1, by=list(tmp1$cluster), length)
		colnames(tmp5) <- c('cluster','n2')
		tmp6 <- merge(tmp4, tmp5, all=T)
		tmp6[is.na(tmp6)] <- 0 
		tmp6$proportion <- tmp6$n1/tmp6$n2
		tmp6 <- tmp6[order(tmp6$proportion, decreasing=T), ]
		tmp7 <- data.frame(Gid=colnames(t_all)[x],
		c1=tmp2$cluster[1], 
		fc_mean12=(tmp2$mean[1]+0.1)/(tmp2$mean[2]+0.1), fc_mean13=(tmp2$mean[1]+0.1)/(tmp2$mean[3]+0.1),
		c2=tmp6$cluster[1], 
		fc_proportion12=(tmp6$proportion[1]+0.01)/(tmp6$proportion[2]+0.01), 
		fc_proportion13=(tmp6$proportion[1]+0.01)/(tmp6$proportion[3]+0.01),
		p=wilcox.test(tmp1$V1[tmp1$cluster==tmp2$cluster[1]],tmp1$V1[tmp1$cluster!=tmp2$cluster[1]])$p.value,
		mean=tmp2$mean[1], proportion=tmp6$proportion[1])
		return(tmp7)
	} 
	else {
		return(NULL)
	}
}
res01 <- do.call('rbind', mclapply(1:dim(t_all)[2], function(x){cseg(x)}, mc.cores=100))
res01 <- res01[complete.cases(res01), ]
#Choosing cluster-specific expressed genes
res02 <- res01[res01$c1==res01$c2 & 
	res01$fc_mean12>2 &
	res01$fc_mean13>5 &
	res01$fc_proportion12>2 &
	res01$fc_proportion13>5 &
	res01$p<1e-3, ]
res03 <- merge(res02, unique(anno_TIGR7[,c('locus','annotation')]), by.x='Gid', by.y='locus') %>% dplyr::arrange(desc(c1))
res04 <- merge(res03, Jiao_marker, by.x='Gid',by.y="row.names")[,c(1,2,11,25)]
#cluster1: mesophyll, cluster2: parenchymal, cluster3: epidermal, cluster4: bulliform, cluster5: primordium 