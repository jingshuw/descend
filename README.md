# descend

R package for DESCEND

DESCEND deconvolves the true gene expression distribution across cells for UMI scRNA-seq counts. It provides estimates of several distribution based statistics (five distribution measurements and the coefficients of covariates (such as batches or cell size)). Based on the estimation, DESCEND also can perform highly variable selection and differential testing of dispersion and burstiness measurements between two groups of cells with covariates adjustment.

## Installation:

```{r eval = F}
library(devtools)
install_github("jingshuw/descend")
```

## Examples:
(for a Windows machine, add the argument: type = "PSOCK" to each of the function that needs parallization)


### For a single cell population, DESCEND can be used to find highly variable genes (HVG).

```{r eval = F}
 library(descend)
 data(zeisel)

 result <- runDescend(zeisel$count.matrix.small, 
                      scaling.consts = zeisel$library.size, n.cores = 3)

 hvg <- findHVG(result)

 hvg$HVG.genes
```

```{r eval = F}
 result1 <- runDescend(zeisel$count.matrix.small, 
                       zeisel$ercc.matrix, ercc.trueMol = zeisel$trueMol,
                       Z0 = log(zeisel$cell.size),
                       n.cores = 3)
 
 ests <- getEstimates(result1)
 ests$CV
```

### For two or more cell popluations, DESCEND can perform differential testing of several distribution measurements between any of the two cell groups with covariates adjustment.

```{r eval = F}
  set.seed(1)

  result.multi <- descendMultiPop(zeisel$count.matrix.small,
                                  labels = zeisel$labels,
                                  scaling.consts = zeisel$library.size,
                                  Z0 = log(zeisel$cell.size), verbose = F, show.message = F,
                                  n.cores = 2)

  ##try 100 null genes first:
  detest.result <- deTest(result.multi, c("endothelial-mural", "pyramidal CA1"),
                          zeisel$count.matrix.small, zeisel$labels,
                          verbose = F, show.message = F,
                          N.genes.null = 100, n.cores = 2)
  
  ## 100 null genes may not get small enough p-values
  detest.result <- deTest.more(result.multi, detest.result, 
                               c("endothelial-mural", "pyramidal CA1"),
                               zeisel$count.matrix.small, labels = zeisel$labels, 
                               N.more.genes = 200, verbose = F, 
                               n.cores = 2)
  
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
```
