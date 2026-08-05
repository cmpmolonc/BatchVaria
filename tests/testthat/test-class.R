test_that("BatchVariaData initialises correctly", {
    bv <- exampleBatchVaria()

    expect_s4_class(bv, "BatchVariaData")

    expect_true("raw" %in% SummarizedExperiment::assayNames(bv))

    expect_true(is.null(S4Vectors::metadata(bv)$variance_history))
})
