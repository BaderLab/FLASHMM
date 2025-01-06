
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FLASHMM

<!-- badges: start -->
<!-- badges: end -->

FLASHMM is a package for analysis of single-cell differential expression
(DE) using a linear mixed- effects model (LMM). The mixed-effects model
has become a powerful tool in single-cell studies due to their ability
to model intra-subject correlation and inter-subject variability.

FLASHMM package provides two functions, lmm and lmmfit, for fitting LMM.
The lmm function uses summary-level statistics as arguments. The lmmfit
function is a wrapper function of lmm, which directly uses cell-level
data and computes the summary statistics inside the function. The lmmfit
function is simple to be operated but it has a limitation of memory use.
For large scale data, it is recommended to precompute the summary
statistics and then use lmm function to fit LMM.

In summary, FLASHMM package provides the following functions.

- lmm: fit LMM using summary-level data.
- lmmfit: fit LMM using cell-level data.
- lmmtest: perform statistical tests on fixed effects and the contrasts
  of the fixed effects.
- sslmm: compute the summary-level data using cell-level data.
- simuRNAseq: simulate multi-sample multi-cell-type scRNA-seq dataset
  based on a negative binomial distribution.

## Installation

You can install the development version of FLASHMM from Github:

``` r
library(devtools)
devtools::install_github("https://github.com/Baderlab/FLASHMM")
```

## Example

This is a basic example which shows you how to use FLASHMM to perform
single-cell differential expression analysis.

``` r
library(FLASHMM)
```

### Simulating a scRNA-seq dataset by simuRNAseq

Simulate a multi-sample multi-cell-cluster scRNA-seq dataset that
contains 25 samples and 4 clusters (cell-types) with 2 treatments.

``` r
set.seed(2412)
dat <- simuRNAseq(nGenes = 50, nCells = 1000, 
                  nsam = 25, ncls = 4, ntrt = 2, nDEgenes = 6)
#> Message: the condition B is treated.

str(dat)
#> List of 5
#>  $ ref.mean.dispersion:'data.frame': 50 obs. of  2 variables:
#>   ..$ mu        : num [1:50] 1.077 0.992 2.815 1.945 2.06 ...
#>   ..$ dispersion: num [1:50] 1.602 2.989 0.954 2.937 1.914 ...
#>  $ metadata           :'data.frame': 1000 obs. of  4 variables:
#>   ..$ sam    : chr [1:1000] "B1" "A6" "A2" "B8" ...
#>   ..$ cls    : chr [1:1000] "C4" "C3" "C1" "C1" ...
#>   ..$ trt    : chr [1:1000] "B" "A" "A" "B" ...
#>   ..$ libsize: num [1:1000] 117 75 101 80 123 113 125 95 77 95 ...
#>  $ counts             : num [1:50, 1:1000] 0 2 6 2 5 2 0 0 2 11 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:50] "Gene1" "Gene2" "Gene3" "Gene4" ...
#>   .. ..$ : chr [1:1000] "Cell1" "Cell2" "Cell3" "Cell4" ...
#>  $ DEgenes            :'data.frame': 6 obs. of  3 variables:
#>   ..$ gene   : chr [1:6] "Gene45" "Gene46" "Gene47" "Gene48" ...
#>   ..$ beta   : num [1:6] 0.55 0.365 0.831 -0.948 -0.875 ...
#>   ..$ cluster: chr [1:6] "C1" "C2" "C3" "C3" ...
#>  $ treatment          : chr "B"
##

#counts and meta data
counts <- dat$counts
metadata <- dat$metadata
rm(dat)
```

### DE analysis using LMM

**Model design**

- Y: gene expression profile (log-transformed counts)
- X: design matrix for fixed effects
- Z: design matrix for random effects

``` r
Y <- log(counts + 1) 
X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = metadata)
Z <- model.matrix(~ 0 + sam, data = metadata)
d <- ncol(Z)
```

**LMM fitting**

1)  Fit LMM by lmmfit using cell-level data.

``` r
fit <- lmmfit(Y, X, Z, d = d)
```

2)  Fit LMM by lmm using summary-level data computed as follows.

``` r
#Computing summary statistics
n <- nrow(X)
XX <- t(X)%*%X; XY <- t(Y%*%X)
ZX <- t(Z)%*%X; ZY <- t(Y%*%Z); ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)

#Fitting LMM
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)

identical(fit, fitss)
#> [1] TRUE
```

3)  Fit LMM by lmm using summary-level data computed by sslmm.

``` r
#Computing summary statistics
ss <- sslmm(X, Y, Z)

#Fitting LMM
fitss <- lmm(summary.stats = ss, d = d)

identical(fit, fitss)
#> [1] TRUE
```

**Hypothesis tests**

``` r
test <- lmmtest(fit)
#head(test)

#t-values
all(t(fit$t) == test[, grep("_t", colnames(test))])
#> [1] TRUE
fit$t[, 1:5]
#>                   Gene1      Gene2      Gene3      Gene4      Gene5
#> log(libsize)  3.5564382  6.1353128  3.6098989  4.1239445  2.8179995
#> clsC1        -2.4158938 -4.8697164 -2.3122737 -2.6709945 -1.3220206
#> clsC2        -2.5957105 -4.9918921 -2.2467924 -2.5108410 -1.2387859
#> clsC3        -2.5576263 -5.0249733 -2.1484996 -2.5856150 -1.2535998
#> clsC4        -2.4466344 -5.0065480 -2.2202673 -2.6983076 -1.2058758
#> clsC1:trtB    0.8697778  0.3659413  0.1515357  0.9707563 -0.9577914
#> clsC2:trtB    2.0700456  1.0976568 -0.1112106  0.7423556 -0.5997369
#> clsC3:trtB    1.5066385  1.5001771 -0.3659261  0.8883383 -1.0410003
#> clsC4:trtB   -0.3263183  1.6810047 -0.4559326  0.6462646 -1.7515738
##

#p-values
all(t(fit$p) == test[, grep("_p", colnames(test))])
#> [1] TRUE
fit$p[, 1:5]
#>                     Gene1        Gene2        Gene3        Gene4       Gene5
#> log(libsize) 0.0003936946 1.226867e-09 0.0003216502 4.036515e-05 0.004928509
#> clsC1        0.0158766072 1.300111e-06 0.0209667982 7.686618e-03 0.186466367
#> clsC2        0.0095791682 7.060831e-07 0.0248727531 1.220264e-02 0.215718133
#> clsC3        0.0106867912 5.971329e-07 0.0319158551 9.862323e-03 0.210283172
#> clsC4        0.0145925607 6.556356e-07 0.0266262016 7.087769e-03 0.228153191
#> clsC1:trtB   0.3846324624 7.144869e-01 0.8795840262 3.319065e-01 0.338401544
#> clsC2:trtB   0.0387066712 2.726210e-01 0.9114719020 4.580478e-01 0.548818677
#> clsC3:trtB   0.1322220329 1.338870e-01 0.7144983040 3.745743e-01 0.298129310
#> clsC4:trtB   0.7442524470 9.307711e-02 0.6485383571 5.182577e-01 0.080156444
```
