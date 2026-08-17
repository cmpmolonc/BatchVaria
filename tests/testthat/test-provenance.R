test_that("provenance returns structured output", {
    bv <- exampleBatchVaria()
    capture.output(p <- provenance(bv))

    expect_type(p, "list")
    expect_named(p, c("assays", "corrections", "variance_profiles"))
    expect_true(is.character(p$assays))
})


test_that("provenance prints expected content", {
    bv <- exampleBatchVaria()
    out <- capture.output(provenance(bv))

    expect_true(length(out) > 0)
    expect_true(any(grepl("assay", tolower(out))))
})


test_that("provenance reflects variance history after profiling", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")

    out <- capture.output(provenance(bv))
    expect_true(any(grepl("variance", tolower(out))))

    p <- provenance(bv)
    expect_s3_class(p$variance_profiles, "data.frame")
    expect_equal(nrow(p$variance_profiles), 1L)
    expect_equal(p$variance_profiles$assay, "raw")
})


test_that("provenance reports the correction ledger", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)
    bv <- runCorrection(bv, method = "combat", batch = "batch")

    p <- provenance(bv)

    expect_s3_class(p$corrections, "data.frame")
    expect_true(all(c("method", "input", "output", "batch", "preserve",
        "no_op") %in% colnames(p$corrections)))
    expect_equal(p$corrections$output, "raw_combat")
})


test_that("provenance reports empty ledgers rather than failing", {
    ## An object nothing has been done to is a legitimate input; the
    ## report should say so rather than error on the absent ledgers.
    set.seed(1)
    mat <- matrix(
        stats::rnorm(100), nrow = 10,
        dimnames = list(paste0("g", seq_len(10)), paste0("s", seq_len(10)))
    )
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = mat))

    out <- capture.output(p <- provenance(se))

    expect_true(any(grepl("Correction history: none", out)))
    expect_true(any(grepl("Variance profiling history: none", out)))
    expect_null(p$corrections)
    expect_null(p$variance_profiles)
})
