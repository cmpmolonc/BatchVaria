test_that("availableCorrectionMethods lists the correction engines", {
    m <- availableCorrectionMethods()
    expect_type(m, "character")
    expect_true("combat" %in% m)
    expect_false(anyDuplicated(m) > 0)
})


test_that("availableVarianceMethods lists the variance engines", {
    m <- availableVarianceMethods()
    expect_type(m, "character")
    expect_setequal(m, c("pca", "anova", "variancePartition"))
})


test_that("varianceHistory returns the ledger", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    ## empty before any profiling
    expect_null(varianceHistory(bv))

    bv <- profileVariance(bv, ~batch, assays = "raw", methods = "anova")
    vh <- varianceHistory(bv)

    expect_type(vh, "list")
    expect_length(vh, 1)
    expect_true(all(c("assay", "formula", "method", "result", "timestamp") %in%
        names(vh[[1]])))
    expect_equal(vh[[1]]$assay, "raw")
    expect_equal(vh[[1]]$method, "anova")

    ## grows with each profiling call
    bv <- profileVariance(bv, ~batch, assays = "raw_center", methods = "anova")
    expect_length(varianceHistory(bv), 2)
})


test_that("varianceResults filters by assay and method", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)
    bv <- profileVariance(
        bv, ~batch,
        assays = c("raw", "raw_center"), methods = c("anova", "pca")
    )

    all_res <- varianceResults(bv)
    expect_setequal(unique(all_res$assay), c("raw", "raw_center"))
    expect_setequal(unique(all_res$method), c("anova", "pca"))

    expect_setequal(unique(varianceResults(bv, assayName = "raw")$assay), "raw")
    expect_setequal(unique(varianceResults(bv, method = "pca")$method), "pca")

    both <- varianceResults(bv, assayName = "raw", method = "anova")
    expect_setequal(unique(both$assay), "raw")
    expect_setequal(unique(both$method), "anova")
})
