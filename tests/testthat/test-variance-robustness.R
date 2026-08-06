## Zero-variance features and single-engine failures are routine on real
## data. Neither should be able to abort a profiling run.

makeConstantFeatures <- function(bv, rows) {
    mat <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    for (i in seq_along(rows)) {
        mat[rows[i], ] <- i - 1
    }
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat),
        colData = SummarizedExperiment::colData(bv)
    )
    BatchVariaData(se)
}


test_that("zero-variance features are excluded, not fatal", {
    set.seed(1)
    bv <- makeConstantFeatures(exampleBatchVaria(nGenes = 60), c(3, 7))

    expect_warning(
        res <- suppressMessages(profileVariance(bv, ~batch)),
        "2 of 60 features.*zero or non-finite variance"
    )

    ## all three engines ran despite the constant features
    methods <- vapply(varianceHistory(res), function(x) x$method, character(1))
    expect_setequal(methods, availableVarianceMethods())
})


test_that("the exclusion is reported once per assay, not once per method", {
    set.seed(1)
    bv <- makeConstantFeatures(exampleBatchVaria(nGenes = 60), c(3, 7))

    warnings <- character()
    withCallingHandlers(
        suppressMessages(profileVariance(bv, ~batch)),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    constantWarnings <- grep("zero or non-finite variance", warnings)
    expect_length(constantWarnings, 1)
})


test_that("feature-counting engines report the filtered feature count", {
    set.seed(1)
    bv <- makeConstantFeatures(exampleBatchVaria(nGenes = 60), c(3, 7))

    res <- suppressWarnings(suppressMessages(profileVariance(bv, ~batch)))

    counts <- vapply(
        Filter(function(x) x$method != "pca", varianceHistory(res)),
        function(x) unique(x$result$n),
        numeric(1)
    )

    expect_true(all(counts == 58))
})


test_that("an assay with no usable features is skipped, not fatal", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    live <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    dead <- live
    dead[] <- 3

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = live, dead = dead),
        colData = SummarizedExperiment::colData(bv)
    )

    expect_warning(
        res <- suppressMessages(
            profileVariance(BatchVariaData(se), ~batch, methods = "anova")
        ),
        "Skipping assay 'dead'"
    )

    assays <- vapply(varianceHistory(res), function(x) x$assay, character(1))
    expect_setequal(assays, "raw")
})


test_that("profileVariance errors when every attempt fails", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    mat <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    mat[] <- 3

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat),
        colData = SummarizedExperiment::colData(bv)
    )

    expect_error(
        suppressWarnings(
            suppressMessages(profileVariance(BatchVariaData(se), ~batch))
        ),
        "All 3 variance profiling attempts failed"
    )
})


test_that("one failing engine does not discard the engines that succeeded", {
    skip_if_not_installed("variancePartition")

    ## Four samples in a 2x2 design leaves no residual degrees of freedom,
    ## so variancePartition cannot fit while pca and anova still can.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 4)

    expect_warning(
        res <- suppressMessages(
            profileVariance(bv, ~ batch + group, assays = "raw")
        ),
        "Variance method 'variancePartition' failed for assay 'raw'"
    )

    methods <- vapply(varianceHistory(res), function(x) x$method, character(1))
    expect_setequal(methods, c("pca", "anova"))
})


test_that(".dropConstantFeatures keeps only informative features", {
    mat <- matrix(
        c(rnorm(20), rep(5, 5), rep(0, 5)),
        nrow = 6,
        byrow = TRUE
    )

    expect_warning(
        kept <- BatchVaria:::.dropConstantFeatures(mat),
        "2 of 6 features"
    )
    expect_equal(nrow(kept), 4)

    expect_error(
        BatchVaria:::.dropConstantFeatures(matrix(1, nrow = 3, ncol = 4)),
        "no variance to decompose"
    )
})


test_that("the anova engine guards against Inf ratios when called directly", {
    ## profileVariance() screens constant features out first, so this
    ## exercises the engine's own guard. Before it, a single constant
    ## feature drove residual_fraction to Inf and the model fraction to
    ## -Inf, surfacing as "variance_fraction must be non-negative".
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 50)
    mat <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    mat[3, ] <- 5

    modelMatrix <- stats::model.matrix(
        ~batch,
        data = as.data.frame(SummarizedExperiment::colData(bv))
    )

    res <- BatchVaria:::.computeAnovaVariance(
        assayMatrix = mat,
        modelMatrix = modelMatrix
    )

    expect_true(all(is.finite(res$variance_fraction)))
    expect_true(all(res$variance_fraction >= 0))
    expect_equal(sum(res$variance_fraction), 1)
})
