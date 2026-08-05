test_that("summary returns structured output", {
    bv <- exampleBatchVaria()
    s <- summary(bv)
    expect_true(is.list(s) || is.data.frame(s))

    ## optional: check key elements if known
    if (is.list(s)) {
        expect_true(length(s) > 0)
    }
})


test_that("summary prints expected content", {
    bv <- exampleBatchVaria()
    out <- capture.output(summary(bv))
    expect_true(length(out) > 0)

    ## loose checks
    expect_true(any(grepl("assay", tolower(out))))
})

test_that("summary reflects variance history after profiling", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    out <- capture.output(summary(bv))

    expect_true(any(grepl("variance", tolower(out))))
})
