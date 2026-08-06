test_that("variance ledger initialises correctly", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, formula = ~ batch + group, methods = c("anova"))

    md <- S4Vectors::metadata(bv)

    expect_true("variance_history" %in% names(md))
    expect_true(is.list(md$variance_history))
})

test_that("recordVariance initialises ledger when NULL", {
    bv <- exampleBatchVaria()

    S4Vectors::metadata(bv)$variance_history <- NULL

    res <- data.frame(
        source = "test",
        term = "batch",
        variance_fraction = 1
    )

    bv <- recordVariance(bv, "raw", ~batch, res)

    vh <- S4Vectors::metadata(bv)$variance_history

    expect_true(length(vh) == 1)
})

test_that("recordVariance warns on duplicate entry", {
    bv <- exampleBatchVaria()

    res <- data.frame(
        source = "test",
        term = "batch",
        variance_fraction = 1
    )

    bv <- recordVariance(bv, "raw", ~batch, res, method = "anova")

    expect_warning(
        recordVariance(bv, "raw", ~batch, res, method = "anova"),
        "already recorded"
    )
})

test_that("varianceTable errors with empty ledger", {
    bv <- exampleBatchVaria()

    S4Vectors::metadata(bv)$variance_history <- list()

    expect_error(
        varianceTable(bv),
        "No variance history"
    )
})

test_that("varianceTable degrades to the percent table without a baseline", {
    bv <- exampleBatchVaria()

    res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(bv, "raw_center", ~batch, res)

    ## 'raw' exists on the object but was not profiled, so it is not among
    ## the assays being summarised and no delta can be formed. That limits
    ## the output rather than invalidating it.
    expect_warning(
        out <- varianceTable(bv, assays = "raw_center"),
        "No baseline assay"
    )
    expect_setequal(out$percent$component, c("batch", "Residuals"))
    expect_null(out$delta)

    ## naming a baseline that is present makes the same call produce one
    expect_type(
        varianceTable(bv, assays = "raw_center", baseline = "raw_center"),
        "list"
    )
})

test_that("varianceTable filters assays correctly", {
    bv <- exampleBatchVaria()

    res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(bv, "raw", ~batch, res)
    bv <- recordVariance(bv, "raw_center", ~batch, res)

    vt <- varianceTable(bv, assays = c("raw", "raw_center"))

    expect_true(all(c("raw", "raw_center") %in% colnames(vt$percent)))
})

test_that("varianceTable formats delta correctly", {
    bv <- exampleBatchVaria()

    res_raw <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    res_adj <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.1, 0.9)
    )

    bv <- recordVariance(bv, "raw", ~batch, res_raw)
    bv <- recordVariance(bv, "raw_center", ~batch, res_adj)

    vt <- varianceTable(bv, formatDelta = TRUE)

    expect_true(is.character(vt$delta[[2]][1]))
})

test_that("varianceChange returns correct structure", {
    bv <- exampleBatchVaria()

    res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(bv, "raw", ~batch, res)
    bv <- recordVariance(bv, "raw_center", ~batch, res)

    vc <- varianceChange(bv, term = "batch")

    expect_true(all(c("assay", "delta_variance") %in% colnames(vc)))
})

test_that("varianceChange errors on missing term", {
    bv <- exampleBatchVaria()

    res <- data.frame(
        source = "test",
        term = "batch",
        variance_fraction = 1
    )

    bv <- recordVariance(bv, "raw", ~batch, res)

    expect_error(
        varianceChange(bv, term = "group"),
        "Term not found"
    )
})


test_that("recordVariance appends to variance ledger", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, formula = ~ batch + group, methods = c("anova"))
    initial_len <- length(S4Vectors::metadata(bv)$variance_history)

    dummy_res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(
        bv,
        assayName = "raw",
        formula = ~batch,
        result = dummy_res,
        method = "anova"
    )

    new_len <- length(S4Vectors::metadata(bv)$variance_history)

    expect_equal(new_len, initial_len + 1)
})


test_that("variance ledger entries have correct structure", {
    bv <- exampleBatchVaria()

    dummy_res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(
        bv,
        assayName = "raw",
        formula = ~batch,
        result = dummy_res,
        method = "anova"
    )

    entry <- tail(S4Vectors::metadata(bv)$variance_history, 1)[[1]]

    expect_true(all(c("assay", "formula", "method", "result") %in% names(entry)))
    expect_equal(entry$assay, "raw")
    expect_equal(entry$method, "anova")
    expect_true(is.data.frame(entry$result))
})


test_that("variance ledger preserves result contract", {
    bv <- exampleBatchVaria()

    dummy_res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.3, 0.7)
    )

    bv <- recordVariance(
        bv,
        assayName = "raw",
        formula = ~batch,
        result = dummy_res,
        method = "anova"
    )

    entry <- tail(S4Vectors::metadata(bv)$variance_history, 1)[[1]]

    expect_true(all(c("source", "term", "variance_fraction") %in% colnames(entry$result)))
})


test_that("multiple variance entries accumulate correctly", {
    bv <- exampleBatchVaria()

    dummy_res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(bv, "raw", ~batch, dummy_res, "anova")
    bv <- recordVariance(bv, "raw", ~batch, dummy_res, "pca")

    vh <- S4Vectors::metadata(bv)$variance_history

    expect_equal(length(vh), 2)
    expect_equal(vh[[1]]$method, "anova")
    expect_equal(vh[[2]]$method, "pca")
})

test_that("variance entries correctly track assay", {
    bv <- exampleBatchVaria()

    dummy_res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(bv, "raw", ~batch, dummy_res, "anova")

    entry <- S4Vectors::metadata(bv)$variance_history[[1]]

    expect_equal(entry$assay, "raw")
})

test_that("variance results can be retrieved by assay and method", {
    bv <- exampleBatchVaria()

    dummy_res <- data.frame(
        source = "test",
        term = c("batch", "Residuals"),
        variance_fraction = c(0.2, 0.8)
    )

    bv <- recordVariance(bv, "raw", ~batch, dummy_res, "anova")

    res <- BatchVaria:::.getVarianceResult(
        bv,
        assayName = "raw",
        method = "anova"
    )

    expect_true(is.data.frame(res))
})


test_that("variance ledger initialises when missing", {
    bv <- exampleBatchVaria()
    bv <- profileVariance(bv, formula = ~ batch + group, methods = c("anova"))

    # manually remove ledger
    md <- S4Vectors::metadata(bv)
    md$variance_history <- NULL
    S4Vectors::metadata(bv) <- md

    dummy_res <- data.frame(
        source = "test",
        term = "batch",
        variance_fraction = 1
    )

    bv <- recordVariance(bv, "raw", ~batch, dummy_res)

    expect_true(!is.null(S4Vectors::metadata(bv)$variance_history))
})


test_that("recordVariance appends multiple entries", {
    bv <- exampleBatchVaria()

    res <- data.frame(
        source = "test",
        term = "batch",
        variance_fraction = 1
    )

    bv <- recordVariance(bv, "raw", ~batch, res)
    bv <- recordVariance(bv, "raw_center", ~batch, res)

    vh <- S4Vectors::metadata(bv)$variance_history

    expect_equal(length(vh), 2)
})


test_that(".getVarianceResult returns the most recent duplicate", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    older <- data.frame(
        source = "anova", term = c("model", "residual"),
        variance_fraction = c(0.2, 0.8)
    )
    newer <- data.frame(
        source = "anova", term = c("model", "residual"),
        variance_fraction = c(0.6, 0.4)
    )

    ## same assay + method twice (different formulas, so no duplicate warning)
    bv <- recordVariance(bv, "raw", ~batch, older, "anova")
    bv <- recordVariance(bv, "raw", ~group, newer, "anova")

    got <- BatchVaria:::.getVarianceResult(bv, assayName = "raw", method = "anova")
    expect_equal(got$variance_fraction, c(0.6, 0.4))

    ## selection is by timestamp, not ledger position: put the newest record
    ## first and it must still win
    vh <- S4Vectors::metadata(bv)$variance_history
    vh[[1]]$timestamp <- as.POSIXct("2026-01-02 00:00:00", tz = "UTC")
    vh[[2]]$timestamp <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
    S4Vectors::metadata(bv)$variance_history <- vh

    got <- BatchVaria:::.getVarianceResult(bv, assayName = "raw", method = "anova")
    expect_equal(got$variance_fraction, c(0.2, 0.8))
})


test_that(".getVarianceResult errors when nothing matches", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    expect_error(
        BatchVaria:::.getVarianceResult(bv, assayName = "raw", method = "anova"),
        "No variance history found"
    )

    bv <- profileVariance(bv, ~batch, assays = "raw", methods = "anova")
    expect_error(
        BatchVaria:::.getVarianceResult(bv, assayName = "raw", method = "pca"),
        "No variance result found for assay 'raw' with method 'pca'"
    )
})
