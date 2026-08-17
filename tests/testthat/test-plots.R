test_that("plotBasisProjection returns list of ggplots", {
    bv <- exampleBatchVaria()
    p <- plotBasisProjection(bv)

    ## is list
    expect_true(is.list(p))
    expect_true(length(p) > 0)

    ## all elements are ggplot
    expect_true(all(vapply(p, inherits, logical(1), "ggplot")))
})


test_that("plotVarianceComposition works with tidy backend", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")

    p <- plotVarianceComposition(bv)

    expect_s3_class(p, "ggplot")
})


test_that("plotVarianceDelta works with tidy evalResult output", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")

    evalResult <- evaluateCorrections(bv)

    p <- plotVarianceDelta(evalResult)

    expect_s3_class(p, "ggplot")
})


test_that("plotSampleDistance works", {
    bv <- exampleBatchVaria()
    p <- plotSampleDistance(bv)
    expect_s3_class(p, "ggplot")
})

test_that("plotBatchEntropy works", {
    bv <- exampleBatchVaria()
    p <- plotBatchEntropy(bv, batchVar = "batch")
    expect_s3_class(p, "ggplot")
})


test_that("plotBasisProjection validates its inputs", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    expect_error(plotBasisProjection(bv, assays = "nope"), "Assays not found: nope")
    expect_error(
        plotBasisProjection(bv, assays = c("raw", "nope", "nah")),
        "Assays not found: nope, nah"
    )
    ## previously failed with "arguments imply differing number of rows"
    expect_error(plotBasisProjection(bv, colourBy = "nope"), "colourBy not found in colData")

    ## Rejecting a non-SummarizedExperiment is now an explicit guard
    ## rather than a dispatch failure, since this is a plain function.
    expect_error(
        plotBasisProjection(list(a = 1)),
        "must be a SummarizedExperiment"
    )
})


test_that("plotPCATrajectory validates its inputs", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    expect_error(plotPCATrajectory(list(a = 1)), "must be a SummarizedExperiment")
    expect_error(
        plotPCATrajectory(bv, assayBefore = "nope"),
        "Assays not found: nope"
    )
    expect_error(
        plotPCATrajectory(bv, assayAfter = "nope"),
        "Assays not found: nope"
    )
    expect_error(
        plotPCATrajectory(bv, assayAfter = "raw_center", colourBy = "nope"),
        "colourBy not found in colData"
    )

    ## the happy path still works
    expect_s3_class(
        plotPCATrajectory(bv, assayBefore = "raw", assayAfter = "raw_center"),
        "ggplot"
    )
})


test_that("PCA plots require at least two samples", {
    mat <- matrix(
        stats::rnorm(10), ncol = 1,
        dimnames = list(paste0("g", seq_len(10)), "s1")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat, raw_x = mat),
        colData = S4Vectors::DataFrame(batch = "A", row.names = "s1")
    )
    bv <- se

    ## previously: "cannot rescale a constant/zero column to unit variance"
    expect_error(plotBasisProjection(bv), "At least two samples")
    expect_error(
        plotPCATrajectory(bv, assayAfter = "raw_x"),
        "At least two samples"
    )
})


test_that("plotVarianceComposition plots an absolute composition", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
    bv <- profileVariance(bv, ~batch, methods = "anova")

    p <- plotVarianceComposition(bv)
    expect_s3_class(p, "ggplot")

    ## the defining property: each assay's terms sum to 100 per cent
    sums <- tapply(p$data$percent, p$data$assay, sum)
    expect_true(all(abs(sums - 100) < 1e-8))

    ## a composition cannot be negative
    expect_true(all(p$data$percent >= 0))

    ## stacked, not dodged, and keyed by assay rather than term
    expect_s3_class(p$layers[[1]]$position, "PositionStack")
    expect_equal(rlang::as_label(p$mapping$x), "assay")
    expect_equal(rlang::as_label(p$mapping$fill), "term")
})


test_that("plotVarianceComposition needs no baseline assay", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(
        bv, ~batch,
        assays = "raw_combat", methods = "anova"
    )

    ## only a corrected assay is profiled; 'raw' has no variance record.
    ## the old delta-based implementation could not draw this at all.
    p <- plotVarianceComposition(bv, assays = "raw_combat")
    expect_s3_class(p, "ggplot")
    expect_equal(as.character(unique(p$data$assay)), "raw_combat")
    expect_true(abs(sum(p$data$percent) - 100) < 1e-8)

    expect_error(varianceDelta(bv, method = "anova"), "No variance result")
})


test_that("plotVarianceComposition and plotVarianceDelta differ", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")

    pc <- plotVarianceComposition(bv)
    pd <- plotVarianceDelta(evaluateCorrections(bv))

    ## these used to render byte-identical plots
    expect_false(identical(ggplot2::layer_data(pc), ggplot2::layer_data(pd)))
    expect_false(identical(pc$labels$title, pd$labels$title))
    expect_s3_class(pc$layers[[1]]$position, "PositionStack")
    expect_s3_class(pd$layers[[1]]$position, "PositionDodge")
})


test_that("plotVarianceComposition validates method and assays", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
    bv <- profileVariance(bv, ~batch, methods = "anova")

    expect_error(
        plotVarianceComposition(bv, method = "pca"),
        "No variance results recorded for method 'pca'"
    )
    expect_error(
        plotVarianceComposition(bv, assays = "nope"),
        "No variance results recorded for assay\\(s\\) nope"
    )
    expect_error(plotVarianceComposition(list(a = 1)), "must be a SummarizedExperiment")

    bare <- exampleBatchVaria(nGenes = 40, nSamples = 8)
    expect_error(plotVarianceComposition(bare), "No variance history found")
})


test_that("plotCorrelationChange works", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")

    p <- plotCorrelationChange(evaluateCorrections(bv))
    expect_s3_class(p, "ggplot")
    expect_true(all(c("assay", "correlation_change") %in% colnames(p$data)))
    expect_false("raw" %in% as.character(p$data$assay))

    ## the defensive branch when the diagnostic is absent
    expect_error(plotCorrelationChange(list()), "'correlation_change' not found")
    expect_error(
        plotCorrelationChange(list(variance_delta = 1)),
        "Did you run evaluateCorrections"
    )
})


test_that("plotVarianceRadar works", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")

    p <- plotVarianceRadar(bv)
    expect_s3_class(p, "ggplot")
    expect_true(all(c("assay", "method", "term", "delta") %in% colnames(p$data)))
    expect_s3_class(p$coordinates, "CoordPolar")
})
