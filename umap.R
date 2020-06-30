
library(Seurat) #this uses Seurat v2.4

#barcodes for cells
wt_cells <- pbmc@cell.names[pbmc@meta.data$rfp == "RFP-"]

#timepoints
tps <- pbmc@meta.data$timepoint[pbmc@meta.data$rfp == "RFP-"]

#greate a new seurat object that is a subset of pbmc input 
wild_type_cells <- CreateSeuratObject(raw.data = data[,wt_cells])

wild_type_cells <- NormalizeData(object = wild_type_cells, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(wild_type_cells@raw.data)

wild_type_cells <- AddMetaData(object = wild_type_cells, metadata = percent.mito, col.name = "percent.mito")

wild_type_cells <- AddMetaData(object = wild_type_cells, metadata = tps, col.name = "timepoint")

wild_type_cells <- ScaleData(object = wild_type_cells, vars.to.regress = c("nUMI", "percent.mito"))

wild_type_cells <- FindVariableGenes(object = wild_type_cells, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

wild_type_cells <- RunPCA(object = wild_type_cells, pc.genes = wild_type_cells@var.genes, do.print = TRUE, pcs.print = 1:5, 
                genes.print = 5)

wild_type_cells <- RunTSNE(object = wild_type_cells, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = wild_type_cells, pt.size = 2)

wild_type_cells <- RunUMAP(object = wild_type_cells, dims.use = 1:10, do.fast = TRUE)
UMAPPlot(object = wild_type_cells, pt.size = 2)
#tsne data


svg(filename = "umap_plot.svg",height=5, width = 6.15, pointsize = 12)
UMAPPlot(object = wild_type_cells, group.by = "condition",  pt.size = 1.6)
dev.off()










