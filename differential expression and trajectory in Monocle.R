library(monocle)
#run monocle size factor and dispersion functios
CMdata <- estimateSizeFactors(CMdata)
CMdata <- estimateDispersions(CMdata)
CMdata <- detectGenes(CMdata, min_expr = 0.1)

#run differential expression gene test
diff_test_res <- differentialGeneTest(CMdata, fullModelFormulaStr = "~timepoint")
#
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
CMdata <- setOrderingFilter(CMdata, ordering_genes)
plot_ordering_genes(CMdata)

#reduce dimension
CMdata <- reduceDimension(CMdata, max_components = 2,
                        method = 'DDRTree')

#order cells
CMdata <- orderCells(CMdata)

#plot trajectory
plot_cell_trajectory(CMdata, color_by = "matt", size=2)

## diff exp for P7

diff_test_res_P7 <- differentialGeneTest(CMdata[, pData(CMdata)$condition %in% c("P07 RFP-","P07 RFP+")], fullModelFormulaStr = "~matt")
sig_genes <- subset(diff_test_res_P7, qval < 0.1)

## Select
sig_genes <- sig_genes[,c("gene_short_name", "pval", "qval")]

#calculate fold change

CMdata_p7 <- CMdata[, pData(CMdata)$condition %in% c("P07 RFP-","P07 RFP+")]
# get average expression for each group
P00wt <-  as.matrix(rowMeans(as.matrix(CMdata_new[, pData(CMdata_new)$condition %in% c("P00 RFP-")]@assayData$exprs_sfnorm)))