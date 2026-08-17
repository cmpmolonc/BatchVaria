## What BatchVaria requires of an object, and what it does not.
##
## There is no container class. The requirements are checked at the point
## of use by .check_se(), so these tests exercise them through public
## functions rather than through a constructor.

test_that("exampleBatchVaria returns an ordinary SummarizedExperiment", {
    bv <- exampleBatchVaria()

    expect_s4_class(bv, "SummarizedExperiment")

    expect_true("raw" %in% SummarizedExperiment::assayNames(bv))

    expect_true(is.null(S4Vectors::metadata(bv)$variance_history))
})


test_that("a plain SummarizedExperiment with no ledgers is accepted", {
    ## The ledgers are what runCorrection() and profileVariance() write.
    ## Requiring one to exist before either has run would make the
    ## package impossible to enter.
    set.seed(1)
    mat <- matrix(
        stats::rnorm(200), nrow = 20,
        dimnames = list(paste0("g", seq_len(20)), paste0("s", seq_len(10)))
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat),
        colData = S4Vectors::DataFrame(
            batch = factor(rep(c("A", "B"), 5)),
            row.names = paste0("s", seq_len(10))
        )
    )

    expect_no_error(assayVariance(se))
    expect_no_error(capture.output(provenance(se)))

    out <- profileVariance(se, ~batch, methods = "anova")
    expect_length(S4Vectors::metadata(out)$variance_history, 1L)

    corrected <- runCorrection(se, method = "combat", batch = "batch")
    expect_true(
        "raw_combat" %in% SummarizedExperiment::assayNames(corrected)
    )
})


test_that("anything that is not a SummarizedExperiment is rejected", {
    expect_error(assayVariance(list()), "must be a SummarizedExperiment")
    expect_error(provenance(data.frame()), "must be a SummarizedExperiment")
    expect_error(varianceResults(NULL), "must be a SummarizedExperiment")
})
