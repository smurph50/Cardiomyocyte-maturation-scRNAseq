# calculate size normalization

#select expression data from monocle object
unnormexp <- as.matrix(HSMM@assayData$exprs)

#select size factors
sizef <- as.matrix(HSMM@phenoData@data$Size_Factor)

# normalize expression
normexp <-  sweep(unnormexp, 2,   sizef, "/")

#add to monocle object
HSMM@assayData$exprs_sfnorm <- normexp