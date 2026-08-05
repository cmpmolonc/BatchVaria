test_that("evaluateCorrections returns structured outputs", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")
    evalResult <- evaluateCorrections(bv)

    ## top-level keys
    expect_true(all(c("variance_delta", "correlation_change", "pca_comparison") %in% names(evalResult)))

    ## varianceDelta structure
    vd <- evalResult$variance_delta
    expect_true(all(c("assay", "method", "term", "delta") %in% colnames(vd)))

    ## correlation structure
    cc <- evalResult$correlation_change
    expect_true(all(c("assay", "correlation_change") %in% colnames(cc)))

    ## pca structure
    pca <- evalResult$pca_comparison
    expect_true(all(c("assay", "component", "variance") %in% colnames(pca)))
})


test_that("compareCorrelations returns valid tidy output", {
    bv <- exampleBatchVaria()
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    df <- compareCorrelations(bv, baseline = "raw")

    ## structure
    expect_true(all(c("assay", "correlation_change") %in% colnames(df)))

    ## baseline excluded
    expect_false("raw" %in% df$assay)

    ## numeric + non-negative
    expect_true(is.numeric(df$correlation_change))
    expect_true(all(df$correlation_change >= 0))
})


test_that("compareCorrelations compares samples, not features", {
    set.seed(101)
    bv <- exampleBatchVaria(nGenes = 120, nSamples = 12)

    raw <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    noise <- as.matrix(SummarizedExperiment::assay(bv, "raw_noise"))

    ## value must equal the mean absolute change over unique sample pairs
    cor_raw <- stats::cor(raw)
    cor_noise <- stats::cor(noise)
    pairs <- upper.tri(cor_raw)

    df <- compareCorrelations(bv, assays = c("raw", "raw_noise"))

    expect_equal(
        df$correlation_change,
        mean(abs(cor_raw[pairs] - cor_noise[pairs]))
    )

    ## per-gene centring leaves gene-gene correlation untouched but does
    ## restructure sample-sample correlation, so it must not report zero
    expect_gt(
        compareCorrelations(bv, assays = c("raw", "raw_center"))$correlation_change,
        0.01
    )

    ## scaling the whole matrix preserves sample correlations exactly
    expect_equal(
        compareCorrelations(bv, assays = c("raw", "raw_scale"))$correlation_change,
        0
    )
})


test_that("compareCorrelations requires at least two samples", {
    mat <- matrix(
        stats::rnorm(20),
        ncol = 1,
        dimnames = list(paste0("g", seq_len(20)), "s1")
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat, corrected = mat),
        colData = S4Vectors::DataFrame(batch = "A", row.names = "s1")
    )

    expect_error(
        compareCorrelations(BatchVariaData(se)),
        "At least two samples"
    )
})


test_that("comparePCA returns tidy PCA summary", {
    bv <- exampleBatchVaria()
    res <- comparePCA(bv)

    ## structure checks
    expect_true(all(c("assay", "component", "variance") %in% colnames(res)))

    ## non-empty
    expect_true(nrow(res) > 0)

    ## expected components present
    expect_true(all(c("PC1", "PC2", "PC3") %in% res$component))

    ## variance is numeric and bounded
    expect_true(is.numeric(res$variance))
    expect_true(all(res$variance >= 0 & res$variance <= 1))

    ## each assay has multiple components
    comps_per_assay <- table(res$assay)
    expect_true(all(comps_per_assay >= 3))
})
