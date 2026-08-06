.onLoad <- function(libname, pkgname) {
    ## Built-in engines register through the same entry point a third party
    ## would use, so the registry is exercised on every load rather than
    ## being a path only external code takes.
    registerVarianceEngine("anova", .computeAnovaVariance)
    registerVarianceEngine("variancePartition", .computeVariancePartitionVariance)

    registerCorrectionMethod("combat", .computeCombat)
    registerCorrectionMethod("limma", .computeLimma)
}
