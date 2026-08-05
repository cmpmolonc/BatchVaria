BatchVaria: A variance-aware framework for evaluating batch correction
================

# Overview

High-throughput transcriptomic datasets frequently contain unwanted
technical variation arising from batch effects. Numerous batch
correction methods exist, but selecting an appropriate method is
challenging because correction can remove biological signal or distort
sample relationships. BatchVaria provides a framework for quantifying
variance components, tracking variance changes after correction, and
evaluating correction performance using multiple metrics, enabling
evidence-based selection of batch correction strategies.

# Installation

The current release of BatchVaria can be installed with:
``` r
# install.packages("devtools")
devtools::install_github("cmpmolonc/BatchVaria")

```
# Example 

Further detail and examples can be found in the preprint at https://doi.org/10.64898/2026.05.07.721996
and in the package introductory vignette:

``` r
library(BatchVaria)
vignette("BatchVaria-quickstart")
```
------------------------------------------------------------------------

# BatchVaria Workflow Summary

The typical BatchVaria workflow is:

1.  Load expression data into BatchVariaData
2.  Apply batch correction methods
3.  Profile variance using profileVariance()
4.  Summarise variance using varianceTable()
5.  Compare corrections using delta variance
6.  Visualise variance composition

# Citation

Please cite the following preprint:

Moir N, Sherwood K, Simpson TI (2026). "BatchVaria: a variance-aware framework for evaluating batch correction in high-dimensional omics data". bioRxiv.
https://doi.org/10.64898/2026.05.07.721996

