#Calculation of total UMI in each cell
umi <- data.frame(cell=colnames(raw_all),umi=colSums(raw_all),cluster=umap_layout$cluster,Sample=umap_layout$Sample)
#Comparison of total UMI differences between samples 
cell_size_ttest <- data.frame() 
for (i in 1:5){
	cell_size_ttest[i,1] <- t.test(umi[umi$cluster==i&umi$Sample=="Ctrl",]$umi,umi[umi$cluster==i&umi$Sample=="LN",]$umi)$p.value
	cell_size_ttest[i+5,1] <-t.test(umi[umi$cluster==i&umi$Sample=="Ctrl",]$umi,umi[umi$cluster==i&umi$Sample=="HS",]$umi)$p.value
}