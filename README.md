
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FLASH-MM

<!-- badges: start -->
<!-- badges: end -->

FLASH-MM is a method (package name: FLASHMM) for analysis of single-cell differential expression
(DE) using a linear mixed- effects model (LMM). The mixed-effects model
has become a powerful tool in single-cell studies due to their ability
to model intra-subject correlation and inter-subject variability.

FLASHMM package provides two functions, lmm and lmmfit, for fitting LMM.
The lmm function uses summary-level statistics as arguments. The lmmfit
function is a wrapper function of lmm, which directly uses cell-level
data and computes the summary statistics inside the function. The lmmfit
function is simple to be operated but it has a limitation of memory use.
For large scale data, it is recommended to precompute and store the
summary statistics and then use lmm function to fit LMM.

In summary, FLASHMM package provides the following functions.

- lmm: fit LMM using summary-level data.
- lmmfit: fit LMM using cell-level data.
- lmmtest: perform statistical tests on fixed effects and the contrasts
  of the fixed effects.
- simuRNAseq: simulate multi-sample multi-cell-type scRNA-seq dataset
  based on a negative binomial distribution.

## Installation

You can install the development version of FLASHMM from Github:

``` r
devtools::install_github("https://github.com/Baderlab/FLASHMM", build_vignettes = TRUE)
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
dat <- simuRNAseq(nGenes = 50, nCells = 1000, nsam = 25, ncls = 4, ntrt = 2, nDEgenes = 6)
#> Message: the condition B is treated.
names(dat)
#> [1] "ref.mean.dispersion" "metadata"            "counts"             
#> [4] "DEgenes"             "treatment"

#counts and meta data
counts <- dat$counts
metadata <- dat$metadata
head(metadata)
#>       sam cls trt libsize
#> Cell1  B1   4   B     117
#> Cell2  A6   3   A      75
#> Cell3  A2   1   A     101
#> Cell4  B8   1   B      80
#> Cell5 B11   4   B     123
#> Cell6  A4   3   A     113
rm(dat)
```

### DE analysis using LMM

**1. Model design**

- Y: gene expression profile (log-transformed counts)
- X: design matrix for fixed effects
- Z: design matrix for random effects

``` r
Y <- log(counts + 1) 
X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = metadata)
Z <- model.matrix(~ 0 + sam, data = metadata)
d <- ncol(Z)
```

**2. LMM fitting**

Option 1: fit LMM by lmmfit function using cell-level data.

``` r
fit <- lmmfit(Y, X, Z, d = d)
```

Option 2: fit LMM by lmm function using summary-level data.

``` r
#(1) Computing summary statistics
n <- nrow(X)
XX <- t(X)%*%X; XY <- t(Y%*%X)
ZX <- t(Z)%*%X; ZY <- t(Y%*%Z); ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)

#(2) Fitting LMM
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)

identical(fit, fitss)
#> [1] TRUE
```

**3. Hypothesis testing**

``` r
##Testing coefficients (fixed effects)
test <- lmmtest(fit)
#head(test)

##The testing t-value and p-values are also provided in the LMM fit.
range(test - cbind(t(fit$coef), t(fit$t), t(fit$p)))
#> [1] 0 0
#fit$coef[, 1:4]
#fit$t[, 1:4]
fit$p[, 1:4]
#>                     Gene1        Gene2        Gene3        Gene4
#> log(libsize) 0.0003936946 1.226867e-09 0.0003216502 4.036515e-05
#> cls1         0.0158766072 1.300111e-06 0.0209667982 7.686618e-03
#> cls2         0.0095791682 7.060831e-07 0.0248727531 1.220264e-02
#> cls3         0.0106867912 5.971329e-07 0.0319158551 9.862323e-03
#> cls4         0.0145925607 6.556356e-07 0.0266262016 7.087769e-03
#> cls1:trtB    0.3846324624 7.144869e-01 0.8795840262 3.319065e-01
#> cls2:trtB    0.0387066712 2.726210e-01 0.9114719020 4.580478e-01
#> cls3:trtB    0.1322220329 1.338870e-01 0.7144983040 3.745743e-01
#> cls4:trtB    0.7442524470 9.307711e-02 0.6485383571 5.182577e-01
##

##Testing contrasts
##We can make comparisons using contrasts. For example, 
##the effects of treatment B vs A in all clusters can be tested 
##using the contrast constructed as follows:
ct <- numeric(ncol(X))
index <- grep("B", colnames(X))
ct[index] <- 1/length(index)
ct
#> [1] 0.00 0.00 0.00 0.00 0.00 0.25 0.25 0.25 0.25
colnames(X)
#> [1] "log(libsize)" "cls1"         "cls2"         "cls3"         "cls4"        
#> [6] "cls1:trtB"    "cls2:trtB"    "cls3:trtB"    "cls4:trtB"

test <- lmmtest(fit, contrast = ct)
head(test)
#>             _coef         _t        _p
#> Gene1  0.09445436  1.4753256 0.1404426
#> Gene2  0.10333114  1.4540794 0.1462409
#> Gene3 -0.02117872 -0.2900354 0.7718498
#> Gene4  0.10281315  0.9531558 0.3407436
#> Gene5 -0.12106061 -1.3918602 0.1642770
#> Gene6  0.06756558  1.1553425 0.2482287
```

## Citation
If you find FLASH-MM useful for your publication, please cite:
Xu & Pouyabahar et al., FLASH-MM: Fast and Scalable Single-Cell Differential Expression Analysis Using Linear Mixed-Effects Models.
