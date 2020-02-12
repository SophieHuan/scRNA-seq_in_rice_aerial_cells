##Construction of the developmental trajectory from primordium cells to epidermal cells 
#Create a Monocle's main class, CellDataSet
#Choosing primordium cells and epidermal cells
cell_select <- umap_layout[umap_layout$cluster==3 |umap_layout$cluster==5,]$cell#"umap_layout" from script1, including cell clustering information 
raw_select <- raw_all[,cell_select]
sample_sheet <- data.frame(row.names = colnames(raw_select), Library= rep("rice",dim(raw_select)[2]))
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
rice <- newCellDataSet(as.matrix(raw_select),featureData = fd, phenoData = pd)
#Estimate size factors and dispersions
rice <- estimateSizeFactors(rice)
rice <- estimateDispersions(rice)
rice <- detectGenes(rice, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(rice), num_cells_expressed > 10))
#Constructing single cell trajectories
ordering_genes <- feature_union
rice_pseu <- setOrderingFilter(rice, ordering_genes)
#Reduce data dimensionality
rice_pseu <- reduceDimension(rice_pseu, reduction_method = "DDRTree")
#Order cells along the trajectory
rice_pseu <- orderCells(rice_pseu)
plot_cell_trajectory(rice_pseu, color_by="Pseudotime")

##Construction of the developmental trajectory from primordium to parenchymal and mesophyll cells
#Create a Monocle's main class, CellDataSet
#Chooing primordium cells, parenchymal cells and mesophyll cells
cell_select <- umap_layout[umap_layout$cluster==1 |umap_layout$cluster==2 |umap_layout$cluster==5,]$cell
raw_select <- raw_all[,cell_select]
sample_sheet <- data.frame(row.names = colnames(raw_select), Library= rep("rice",dim(raw_select)[2]))
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
#Estimate size factors and dispersions
rice <- newCellDataSet(as.matrix(raw_select),featureData = fd,phenoData = pd)
rice <- estimateSizeFactors(rice)
rice <- estimateDispersions(rice)
rice <- detectGenes(rice, min_expr = 0.1)
#Constructing single cell trajectories
ordering_genes <- feature_union
rice_pseu <- setOrderingFilter(rice, ordering_genes)
#Reduce data dimensionality
rice_pseu <- reduceDimension(rice_pseu, reduction_method = "DDRTree")
#Order cells along the trajectory
rice_pseu <- orderCells(rice_pseu)
plot_cell_trajectory(rice_pseu, color_by="Pseudotime")
#Clustering genes by pseudotemporal expression pattern
diff_test_res <- differentialGeneTest(raw_all[feature_union,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(rice_pseu[sig_gene_names,],num_clusters = 4,cores = 1,show_rownames = F)

##Calculation of the average expression levels of histone-encoding genes in each cell
histone_TIGR7 <- TIGR7[grep("Core histone|histone H",TIGR7$annotation),]  
#Take log of CPM
cpm_all <- as.data.frame(t(t(raw_all)*1e6/rowSums(t(raw_all))))
cpm_all_log <- log10(cpm_all+1)
#Calculation of expression levels of histone-encoding genes 
histone <- cpm_all_log[rownames(cpm_all_log) %in% hist_anno_TIGR7$locus,]
histone_score <- data.frame(score= colSums(histone)/66)