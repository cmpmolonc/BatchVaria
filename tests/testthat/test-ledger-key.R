## The ledger is keyed on (assay, method, formula, term). Summaries must
## resolve to exactly one decomposition rather than flattening the ledger.

test_that("formula keys ignore term order but not term content", {
    expect_equal(
        BatchVaria:::.formulaKey(~ batch + group),
        BatchVaria:::.formulaKey(~ group + batch)
    )

    expect_false(identical(
        BatchVaria:::.formulaKey(~batch),
        BatchVaria:::.formulaKey(~group)
    ))

    ## intercept and interactions are part of the model's identity
    expect_false(identical(
        BatchVaria:::.formulaKey(~batch),
        BatchVaria:::.formulaKey(~ 0 + batch)
    ))
    expect_false(identical(
        BatchVaria:::.formulaKey(~ batch + group),
        BatchVaria:::.formulaKey(~ batch * group)
    ))

    ## interaction terms sort deterministically
    expect_equal(
        BatchVaria:::.formulaKey(~ batch * group),
        BatchVaria:::.formulaKey(~ group * batch)
    )
})


test_that("recordVariance stores the canonical formula key", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(bv, ~ batch + group, methods = "anova")

    entry <- varianceHistory(bv)[[1]]
    expect_equal(entry$formula_key, "~batch + group")
})


test_that("reprofiling the same model warns; a different model does not", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    bv <- profileVariance(bv, ~batch, assays = "raw", methods = "anova")

    ## same model written with terms reordered is still the same model
    bv2 <- profileVariance(bv, ~ batch + group, assays = "raw", methods = "anova")
    expect_warning(
        profileVariance(bv2, ~ group + batch, assays = "raw", methods = "anova"),
        "already recorded"
    )

    ## a genuinely different model must not warn
    expect_no_warning(
        profileVariance(bv, ~group, assays = "raw", methods = "anova")
    )
})


test_that("varianceTable refuses to mix methods", {
    ## Previously this silently stacked PC1..PCn, model/residual and
    ## batch/Residuals in one 'component' column, producing assay columns
    ## that summed to roughly 300%.
    skip_if_not_installed("variancePartition")

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)
    bv <- suppressWarnings(suppressMessages(profileVariance(bv, ~batch)))

    expect_error(
        varianceTable(bv),
        "holds results for 2 methods"
    )

    tab <- varianceTable(bv, method = "anova")$percent
    expect_setequal(tab$component, c("batch", "shared", "residual"))

    ## every assay column is now a genuine composition
    for (col in setdiff(colnames(tab), "component")) {
        expect_equal(sum(tab[[col]]), 100, tolerance = 1e-6)
    }
})


test_that("varianceTable refuses to mix formulas", {
    ## Previously this produced 'Values are not uniquely identified',
    ## list-columns, and then 'non-numeric argument to binary operator'.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(bv, ~batch, methods = "anova")
    bv <- profileVariance(bv, ~group, methods = "anova")

    expect_error(
        varianceTable(bv),
        "has results for 2 formulas"
    )

    tab <- varianceTable(bv, formula = ~group)$percent
    expect_true(all(vapply(
        tab[, -1, drop = FALSE],
        is.numeric,
        logical(1)
    )))
})


test_that("repeat profiling of one key resolves to the most recent result", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(bv, ~batch, assays = "raw", methods = "anova")

    ## a second, deliberately distinguishable record under the same key
    stub <- data.frame(
        source = "anova",
        term = c("model", "residual"),
        variance_fraction = c(0.25, 0.75),
        stringsAsFactors = FALSE
    )
    bv <- suppressWarnings(
        recordVariance(bv, "raw", ~batch, stub, method = "anova")
    )

    expect_length(varianceHistory(bv), 2)

    tab <- varianceTable(bv, baseline = "raw")$percent
    expect_equal(tab$raw[tab$component == "model"], 25)
})


test_that("varianceTable accepts a baseline other than raw", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(
        bv, ~batch,
        assays = c("raw", "raw_noise"), methods = "anova"
    )

    res <- varianceTable(
        bv,
        assays = c("raw", "raw_noise"),
        baseline = "raw_noise"
    )

    ## the baseline column is dropped from the delta table, not 'raw'
    expect_true("raw" %in% colnames(res$delta))
    expect_false("raw_noise" %in% colnames(res$delta))
})


test_that("varianceDelta compares like with like across formulas", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(
        bv, ~batch,
        assays = c("raw", "raw_noise"), methods = "anova"
    )
    bv <- profileVariance(
        bv, ~group,
        assays = c("raw", "raw_noise"), methods = "anova"
    )

    expect_error(
        varianceDelta(bv, assays = c("raw", "raw_noise")),
        "has results for 2 formulas"
    )

    res <- varianceDelta(
        bv,
        assays = c("raw", "raw_noise"),
        formula = ~batch
    )
    expect_true(all(is.finite(res$delta)))
})


test_that("varianceTable works on objects with no assay named raw", {
    ## "raw" is a convention of exampleBatchVaria(), not a requirement of
    ## the container. An object built from a user's own SummarizedExperiment
    ## previously lost the percent table entirely.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(SummarizedExperiment::assay(bv, "raw"))),
        colData = SummarizedExperiment::colData(bv)
    )
    obj <- profileVariance(BatchVariaData(se), ~batch, methods = "anova")

    expect_warning(res <- varianceTable(obj), "No baseline assay")

    ## the composition does not depend on a baseline, so it is still returned
    expect_setequal(res$percent$component, c("batch", "shared", "residual"))
    expect_true("counts" %in% colnames(res$percent))
    expect_null(res$delta)
    expect_null(res$baseline)
})


test_that("the baseline is inferred from the correction lineage", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(SummarizedExperiment::assay(bv, "raw"))),
        colData = SummarizedExperiment::colData(bv)
    )
    obj <- BatchVariaData(se)
    obj <- runCorrection(obj, method = "combat", batch = "batch", assayName = "counts")
    obj <- profileVariance(obj, ~batch, methods = "anova")

    res <- varianceTable(obj)

    ## resolved from correction_history, not from the name "raw"
    expect_equal(res$baseline, "counts")
    expect_true("counts_combat" %in% colnames(res$delta))
    expect_false("counts" %in% colnames(res$delta))
})


test_that("an explicitly named absent baseline is an error, not an inference", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(bv, ~batch, assays = c("raw", "raw_noise"), methods = "anova")

    expect_error(
        varianceTable(bv, assays = c("raw", "raw_noise"), baseline = "nope"),
        "Baseline assay 'nope' is not among the assays"
    )
})


test_that("varianceChange refuses without a baseline", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(SummarizedExperiment::assay(bv, "raw"))),
        colData = SummarizedExperiment::colData(bv)
    )
    obj <- profileVariance(BatchVariaData(se), ~batch, methods = "anova")

    expect_error(
        suppressWarnings(varianceChange(obj, term = "residual")),
        "variance change cannot be computed"
    )
})


test_that("recordVariance enforces the result contract", {
    ## recordVariance() is exported, so it bypasses the check
    ## profileVariance() applies to engine output.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    expect_error(
        recordVariance(bv, "raw", ~batch, data.frame(nonsense = 1:2), method = "anova"),
        "missing required columns"
    )

    expect_error(
        recordVariance(
            bv, "raw", ~batch,
            data.frame(source = "x", term = "t", variance_fraction = -0.5),
            method = "anova"
        ),
        "must be non-negative"
    )

    ## a valid result is still accepted
    ok <- data.frame(
        source = "anova",
        term = c("model", "residual"),
        variance_fraction = c(0.3, 0.7)
    )
    expect_s4_class(
        recordVariance(bv, "raw", ~batch, ok, method = "anova"),
        "BatchVariaData"
    )
})


test_that("recency resolution survives a timestamp tie", {
    ## Sys.time() is sub-second, so real runs differ -- but two entries
    ## sharing a timestamp must still resolve to the appended-later one.
    ## which.max() returns the first maximum, which would be the stale one.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    bv <- profileVariance(bv, ~batch, assays = "raw", methods = "anova")
    bv <- suppressWarnings(
        profileVariance(bv, ~batch, assays = "raw", methods = "anova")
    )

    vh <- S4Vectors::metadata(bv)$variance_history
    expect_length(vh, 2)

    ## force the tie, and make the stale entry recognisable
    vh[[1]]$timestamp <- vh[[2]]$timestamp
    vh[[1]]$result$variance_fraction <- rep(0.99, nrow(vh[[1]]$result))
    S4Vectors::metadata(bv)$variance_history <- vh

    res <- BatchVaria:::.getVarianceResult(bv, "raw", "anova")
    expect_false(isTRUE(all.equal(res$variance_fraction[1], 0.99)))
})
