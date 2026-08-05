test_that("profileVariance records results correctly", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~ batch + group, methods = "anova")
    vh <- S4Vectors::metadata(bv)$variance_history
    expect_true(length(vh) >= 1)
    res <- varianceResults(bv)
    expect_true(length(res) >= 1)
})

test_that("varianceDelta returns tidy structure", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")
    vd <- varianceDelta(bv, method = "anova")

    ## structure
    expect_true(all(c("assay", "method", "term", "delta") %in% colnames(vd)))

    ## non-empty
    expect_true(nrow(vd) > 0)

    ## no baseline assay present
    expect_false("raw" %in% vd$assay)

    ## numeric delta
    expect_true(is.numeric(vd$delta))
})
