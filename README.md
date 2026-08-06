BatchVaria: A variance-aware framework for evaluating batch correction
================

# Overview

High-throughput transcriptomic datasets frequently contain unwanted
technical variation arising from batch effects. Numerous batch
correction methods exist, but selecting an appropriate method is
challenging because correction can remove biological signal or distort
sample relationships.

BatchVaria quantifies variance components before and after correction and
records every decomposition together with the total variance it is relative
to and the share the experimental design cannot attribute to any single
term. That provenance is what lets it distinguish a correction that removed
unwanted variation from one that merely removed variation.

# Installation

The current release of BatchVaria can be installed with:
``` r
# install.packages("devtools")
devtools::install_github("cmpmolonc/BatchVaria")
```

# Quick start

``` r
library(BatchVaria)

set.seed(1)
bv <- exampleBatchVaria(confounding = 0.6)

bv <- runCorrection(
    bv,
    method = "combat", batch = "batch",
    newAssayName = "corrected"
)

bv <- profileVariance(
    bv, ~ batch + group,
    assays = c("raw", "corrected"),
    methods = c("anova", "variancePartition")
)
```

Summaries describe one decomposition at a time. With more than one engine in
the ledger, `varianceTable()` asks which you mean rather than mixing results
whose terms are not comparable:

``` r
varianceTable(bv, assays = c("raw", "corrected"))
#> Error: The variance ledger holds results for 2 methods (anova,
#> variancePartition). Their terms are not comparable, so specify method =
#> to choose one

varianceTable(bv, assays = c("raw", "corrected"), method = "anova")$percent
#> # A tibble: 4 x 3
#>   component   raw corrected
#> 1 batch      11.0       5.6
#> 2 group      11.3      15.2
#> 3 shared     11.4      -3.5
#> 4 residual   66.3      82.7
```

The same applies to `formula` when one assay has been profiled under more
than one model.

# Reading the table

Those fractions sum to 100 within each assay, so on their own they cannot be
compared across assays. Here the batch fraction halves and the group fraction
*rises*, which reads as a successful correction. It is not one:

``` r
assayVariance(bv, assays = c("raw", "corrected"))[, c("assay", "total_variance")]
#>       assay total_variance
#> 1       raw          1281.
#> 2 corrected           829.
```

Total variance fell by 35%. The group fraction rose because the denominator
collapsed, not because biological signal was preserved: in absolute terms it
fell from 145 to 126. Meanwhile 5.6% of the batch effect remains, and
`shared` — the variance this confounded design cannot attribute to either
term — has gone negative, indicating the correction distorted the design's
attributability.

Multiplying a fraction by its assay's total variance gives a quantity that is
comparable across assays. `basisRetention()` adds a second view, of how much
of the original variance structure a correction preserved.

# Workflow

1.  Load expression data into a `BatchVariaData` container
2.  Apply corrections with `runCorrection()`, naming any variables to
    `preserve`
3.  Decompose variance with `profileVariance()`
4.  Summarise with `varianceTable(method = ...)`, naming the engine and,
    where ambiguous, the formula
5.  Convert fractions to absolute variance using `assayVariance()`
6.  Compare structure with `basisRetention()`, `comparePCA()` and
    `compareCorrelations()`
7.  Visualise with `plotVarianceComposition()`, `plotVarianceDelta()` and
    the other diagnostic plots

Both layers are extensible behind documented contracts:
`registerVarianceEngine()` and `registerCorrectionMethod()` add an engine or
a correction method, and `availableVarianceMethods()` and
`availableCorrectionMethods()` list those registered. `combat` and `limma`
ship as built-ins, registered through the same public entry point a third
party would use.

# Further detail

Further detail and examples can be found in the preprint at
https://doi.org/10.64898/2026.05.07.721996 and in the package vignette:

``` r
vignette("BatchVaria-quickstart")
```

# Citation

Please cite the following preprint:

Moir N, Sherwood K, Simpson TI (2026). "BatchVaria: a variance-aware framework for evaluating batch correction in high-dimensional omics data". bioRxiv.
https://doi.org/10.64898/2026.05.07.721996
