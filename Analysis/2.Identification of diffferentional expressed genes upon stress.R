##Identification of diffferentional expressed genes (DEGs) upon stress
#Gene expression levels in cluster 1 of Ctrl sample
cluster1_Ctrl_raw <- raw_Ctrl[, umap_layout[umap_layout$cluster==1&umap_layout$Sample=="Ctrl",]$cell]#"umap_layout" from script1, including cell clustering information 
cluster1_Ctrl_cpm <- as.data.frame(t(t(cluster1_Ctrl_raw)*1e6/rowSums(t(cluster1_Ctrl_raw))))
gene_select_Ctrl_cluster1 <- rownames(cluster1_Ctrl_cpm[apply(cluster1_Ctrl_cpm,1,function(x){mean(x)})>1,])
#Gene expression levels in cluster 1 of LN sample
cluster1_LN_raw <- raw_LN[, umap_layout[umap_layout$cluster==1&umap_layout$Sample=="LN",]$cell]
cluster1_LN_cpm <- as.data.frame(t(t(cluster1_LN_raw)*1e6/rowSums(t(cluster1_LN_raw))))
gene_select_LN_cluster1 <- rownames(cluster1_LN_cpm[apply(cluster1_LN_cpm,1,function(x){mean(x)})>1,])
#Choose union of genes that expressed in above
gene_select_union_cluster1 <- union(gene_select_Ctrl_cluster1, gene_select_LN_cluster1)
cluster1_Ctrl <- cpm_Ctrl[gene_select_union_cluster1, umap_layout[umap_layout$cluster==1&umap_layout$Sample=="Ctrl",]$cell]
cluster1_LN <- cpm_LN[gene_select_union_cluster1, umap_layout[umap_layout$cluster==1&umap_layout$Sample=="LN",]$cell]
#Defination of DEGS in cluster 1 upon LN stress
deg <- function(x){
	FC <- ((mean(as.numeric(cluster1_LN[x,])))/mean(as.numeric(cluster1_Ctrl[x,])))
	p.value <- wilcox.test(as.numeric(cluster1_Ctrl[x,]), 
												 as.numeric(cluster1_LN[x,]))$p.value
	Diff <- (log2(mean(as.numeric(cluster1_LN[x,]))+1)) - log2(mean(as.numeric(cluster1_Ctrl[x,]))+1)
	return(data.frame(gene=rownames(cluster1_Ctrl)[x], FC = FC, p.value = p.value, Dif =Dif))
}
deg_cluster1 <- do.call('rbind', mclapply(1:dim(cluster1_Ctrl)[1], deg, mc.cores = 100))
deg_cluster1$Group <- "N"
deg_cluster1$Group[abs(log2(heatmap_cluster1$FC))>=1& -log10(deg_cluster1$p.value)> 2] <- "Y"
degs_c1_LN <- deg_cluster1[deg_cluster1$Group=="Y",]$gene
#DEGs in other cell types and upon HS stress are also identified by above method 

##Upset plot
#"1" represents DEGS and "0" represnts Non-DEGs
upset_LN <- data.frame(row.names = rownames(raw_all))
upset_LN$cluster1 <- 0
upset_LN$cluster1[rownames(upset_LN) %in% degs_c1_LN] <- 1
upset_LN$cluster2 <- 0
upset_LN$cluster2[rownames(upset_LN) %in% degs_c2_LN] <- 1
upset_LN$cluster3 <- 0
upset_LN$cluster3[rownames(upset_LN) %in% degs_c3_LN] <- 1
upset_LN$cluster4 <- 0
upset_LN$cluster4[rownames(upset_LN) %in% degs_c4_LN] <- 1
upset_LN$cluster5 <- 0
upset_LN$cluster5[rownames(upset_LN) %in% degs_c5_LN] <- 1
library(UpSetR)
upset(upset_LN, sets = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5"), order.by = "freq")
#The upset plot of DEGs upon HS is also obtained using above method.