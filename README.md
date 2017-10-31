# descend

R package for DESCEND

DESCEND deconvolves the true gene expression distribution across cells for UMI scRNA-seq counts. It provides estimates of several distribution based statistics (five distribution measurements and the coefficients of covariates (such as batches or cell size)). Based on the estimation, DESCEND also can perform highly variable selection and differential testing of dispersion and burstiness measurements between two groups of cells with covariates adjustment.


Examples (for a Windows machine, add the function type = "PSOCK" to each of the function that needs parallization):

1. For a single cell population. DESCEND can be used to find highly variable genes (HVG).

 data(zeisel)
 result <- runDescend(zeisel$count.matrix.small, 
                      scaling.consts = zeisel$library.size, n.cores = 3)
 hvg <- findHVG(result)
 hvg$HVG.genes

2. For two or more cell popluations, DESCEND can perform differential testing of several distribution measurementsbetween any of the two cell groups with covariates adjustment.

  data(zeisel)
  set.seed(1)
  result.multi <- descendMultiPop(zeisel$count.matrix.small,
                                  labels = zeisel$labels,
                                  scaling.consts = zeisel$library.size,
                                  Z0 = log(zeisel$cell.size), verbose = F, show.message = F,
                                  n.cores = 3)
  ## try 100 null genes first
  detest.result <- deTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
                          zeisel$count.matrix.small, zeisel$labels,
                          verbose = F, show.message = F,
                          N.genes.null = 100, n.cores = 3)
  
  ## 100 null genes may not get small enough p-values
  detest.result <- deTest.more(result.multi, detest.result, 
                               c("endothelial-mural", "pyramidal CA1"),
                               zeisel$count.matrix.small, labels = zeisel$labels, 
                               N.more.genes = 200, verbose = F, 
                               n.cores = 3)
  
  layout(matrix(1:4, nrow = 2))
  de.scores1 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
                          detest.result, measurement.name = "Gini", alpha = 0.05)
  de.scores2 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
                          detest.result, measurement.name = "Active Intensity", 
                          alpha = 0.05, log = "xy")
  de.scores3 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
                          detest.result, measurement.name = "Active Fraction", alpha = 0.1)
  de.scores4 <- plotDeTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
                          detest.result, measurement.name = "Adjusted Active Fraction", alpha = 0.1)


