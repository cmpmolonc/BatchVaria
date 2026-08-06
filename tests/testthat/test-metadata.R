test_that("metadata initialises without variance history", {
    bv <- exampleBatchVaria()
    md <- S4Vectors::metadata(bv)

    ## either NULL or empty list depending on implementation
    expect_true(
        is.null(md$variance_history) ||
            length(md$variance_history) == 0
    )
})

test_that("profileVariance records variance history", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    md <- S4Vectors::metadata(bv)

    expect_true("variance_history" %in% names(md))
    expect_true(length(md$variance_history) == 2)

    entry <- md$variance_history[[1]]

    ## structural checks (loose but meaningful)
    expect_true("assay" %in% names(entry))
    expect_true("formula" %in% names(entry))
    expect_true("method" %in% names(entry))
})

test_that("variance history accumulates across calls", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    bv <- profileVariance(bv, ~ batch + group, assayName = "raw")

    md <- S4Vectors::metadata(bv)

    expect_equal(length(md$variance_history), 4)
})


test_that("variance history tracks multiple assays", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat")
    md <- S4Vectors::metadata(bv)

    assays_recorded <- vapply(
        md$variance_history,
        function(x) x$assay,
        character(1)
    )

    expect_true("raw" %in% assays_recorded)
    expect_true(any(grepl("combat", assays_recorded)))
})

test_that("varianceResults accessor aligns with metadata", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    vr <- varianceResults(bv)
    md <- S4Vectors::metadata(bv)

    ## long format: one row per (ledger entry x model term)
    expected_rows <- sum(vapply(
        md$variance_history,
        function(x) nrow(x$result),
        integer(1)
    ))

    expect_s3_class(vr, "data.frame")
    expect_true(nrow(vr) == expected_rows)
    expect_true(
        all(
            c("source", "term", "variance_fraction", "assay", "method") %in%
                colnames(vr)
        )
    )
    expect_setequal(
        unique(vr$assay),
        unique(vapply(md$variance_history, function(x) x$assay, character(1)))
    )
})

test_that("variance history includes timestamps", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    entry <- S4Vectors::metadata(bv)$variance_history[[1]]

    expect_true("timestamp" %in% names(entry))
    expect_true(inherits(entry$timestamp, "POSIXct"))
})


test_that("variance history entries have consistent structure", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, ~batch, assayName = "raw")
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(bv, ~batch, assayName = "raw_combat")
    vh <- S4Vectors::metadata(bv)$variance_history

    keys <- lapply(vh, names)

    ## all entries share same structure
    expect_true(all(vapply(keys, function(k) identical(k, keys[[1]]), logical(1))))
})
