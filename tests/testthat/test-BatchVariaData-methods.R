test_that("show method runs without error", {
    bv <- exampleBatchVaria()
    expect_no_error(show(bv))
})


test_that("show method prints expected content", {
    bv <- exampleBatchVaria()
    out <- capture.output(show(bv))
    expect_true(length(out) > 0)

    ## light structural checks (not brittle)
    expect_true(any(grepl("BatchVariaData", out)))
    expect_true(any(grepl("assay", tolower(out))))
})

test_that("print method works", {
    bv <- exampleBatchVaria()
    expect_no_error(print(bv))
})
