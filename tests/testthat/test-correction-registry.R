## A method registered through the public API only, indistinguishable from a
## built-in everywhere downstream. Mirrors test-engine-registry.R.

## shifts every batch onto the overall feature mean
centreByBatch <- function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
    b <- factor(sampleData[[batch]])
    out <- assayMatrix
    for (lv in levels(b)) {
        idx <- which(b == lv)
        out[, idx] <- out[, idx] - rowMeans(out[, idx, drop = FALSE])
    }
    out + rowMeans(assayMatrix)
}


test_that("built-in methods are registered through the public entry point", {
    expect_setequal(availableCorrectionMethods(), c("combat", "limma"))
})


test_that("a third-party method is indistinguishable from a built-in", {
    on.exit(unregisterCorrectionMethod("batchCentre"), add = TRUE)
    registerCorrectionMethod("batchCentre", centreByBatch)

    expect_true("batchCentre" %in% availableCorrectionMethods())

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 80, nSamples = 12)
    bv <- runCorrection(
        bv,
        method = "batchCentre", batch = "batch", preserve = "group",
        newAssayName = "centred"
    )

    ## the assay arrives like any other
    expect_true("centred" %in% SummarizedExperiment::assayNames(bv))

    ## the ledger records the method and the arguments it was called with
    ch <- S4Vectors::metadata(bv)$correction_history
    entry <- ch[[length(ch)]]
    expect_equal(entry$method, "batchCentre")
    expect_equal(entry$batch, "batch")
    expect_equal(entry$preserve, "group")
    expect_equal(entry$assay_in, "raw")
    expect_equal(entry$assay_out, "centred")
    expect_true(is.numeric(entry$total_variance_in))
    expect_true(is.numeric(entry$total_variance_out))

    ## lineage resolution finds the baseline without being told
    expect_equal(
        BatchVaria:::.resolveBaseline(bv, c("raw", "centred")),
        "raw"
    )

    ## and the delta view works off that lineage
    bv <- profileVariance(
        bv, ~ batch + group,
        assays = c("raw", "centred"), methods = "anova"
    )
    tab <- varianceTable(bv, assays = c("raw", "centred"), method = "anova")
    expect_false(is.null(tab$delta))
    expect_equal(tab$baseline, "raw")
})


test_that("registration rejects a non-conforming signature", {
    expect_error(
        registerCorrectionMethod("noDots", function(assayMatrix, batch,
                                                    sampleData, preserve) {
            assayMatrix
        }),
        "must accept '\\.\\.\\.'"
    )

    ## a method carrying only '...' would absorb 'preserve' silently
    expect_error(
        registerCorrectionMethod(
            "noPreserve",
            function(assayMatrix, batch, sampleData, ...) assayMatrix
        ),
        "missing preserve"
    )

    expect_error(
        registerCorrectionMethod("notAFunction", "combat"),
        "must be a function"
    )

    expect_error(
        registerCorrectionMethod("", function(assayMatrix, batch, sampleData,
                                              preserve = NULL, ...) {
            assayMatrix
        }),
        "single non-empty string"
    )

    expect_false(any(
        c("noDots", "noPreserve", "notAFunction") %in%
            availableCorrectionMethods()
    ))
})


test_that("registration refuses to clobber without overwrite", {
    ok <- function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
        assayMatrix
    }

    expect_error(registerCorrectionMethod("combat", ok), "already registered")

    ## and the built-in is untouched by the failed attempt
    expect_identical(
        BatchVaria:::.getCorrectionMethod("combat"),
        BatchVaria:::.computeCombat
    )
})


test_that("a non-conforming return is caught at use time", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    ## conforming signature, non-matrix return
    on.exit(unregisterCorrectionMethod("returnsList"), add = TRUE)
    registerCorrectionMethod(
        "returnsList",
        function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
            list(corrected = assayMatrix)
        }
    )
    expect_error(
        runCorrection(bv, method = "returnsList", batch = "batch"),
        "must return a numeric matrix"
    )

    ## right shape, reordered samples
    on.exit(unregisterCorrectionMethod("reorders"), add = TRUE)
    registerCorrectionMethod(
        "reorders",
        function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
            assayMatrix[, rev(seq_len(ncol(assayMatrix))), drop = FALSE]
        }
    )
    expect_error(
        runCorrection(bv, method = "reorders", batch = "batch"),
        "returned different dimnames"
    )

    ## right shape and order, dimnames stripped
    on.exit(unregisterCorrectionMethod("stripsNames"), add = TRUE)
    registerCorrectionMethod(
        "stripsNames",
        function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
            unname(assayMatrix)
        }
    )
    expect_error(
        runCorrection(bv, method = "stripsNames", batch = "batch"),
        "returned different dimnames"
    )

    ## wrong dimensions
    on.exit(unregisterCorrectionMethod("dropsFeatures"), add = TRUE)
    registerCorrectionMethod(
        "dropsFeatures",
        function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
            assayMatrix[seq_len(nrow(assayMatrix) - 1L), , drop = FALSE]
        }
    )
    expect_error(
        runCorrection(bv, method = "dropsFeatures", batch = "batch"),
        "returned a matrix of"
    )
})


test_that("unregistering removes the method and reports unknown names", {
    on.exit(
        {
            if (!"combat" %in% availableCorrectionMethods()) {
                registerCorrectionMethod("combat", BatchVaria:::.computeCombat)
            }
        },
        add = TRUE
    )

    unregisterCorrectionMethod("combat")
    expect_false("combat" %in% availableCorrectionMethods())

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 30, nSamples = 8)
    expect_error(
        runCorrection(bv, method = "combat", batch = "batch"),
        "Unknown correction method: combat"
    )

    expect_error(unregisterCorrectionMethod("combat"), "No correction method")

    registerCorrectionMethod("combat", BatchVaria:::.computeCombat)
    expect_true("combat" %in% availableCorrectionMethods())
})


test_that("a correction that changes nothing is flagged, not refused", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    on.exit(unregisterCorrectionMethod("passthrough"), add = TRUE)
    registerCorrectionMethod(
        "passthrough",
        function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
            assayMatrix
        }
    )

    expect_warning(
        out <- runCorrection(
            bv, method = "passthrough", batch = "batch",
            newAssayName = "unchanged"
        ),
        "indistinguishable from 'raw'"
    )

    ## the assay is still created -- flagged, not refused
    expect_true("unchanged" %in% SummarizedExperiment::assayNames(out))

    ## and the condition travels with the object rather than living in a
    ## console message that a later session never sees
    ch <- S4Vectors::metadata(out)$correction_history
    expect_true(ch[[length(ch)]]$no_op)

    ## a real correction is not flagged (ComBat reports progress through
    ## message(), so this asserts on warnings rather than on silence)
    expect_no_warning(
        good <- suppressMessages(runCorrection(
            bv, method = "combat", batch = "batch", newAssayName = "corrected"
        ))
    )
    ch2 <- S4Vectors::metadata(good)$correction_history
    expect_false(ch2[[length(ch2)]]$no_op)
})


test_that("summary() shows what distinguishes two runs of one method", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 12)
    bv <- runCorrection(
        bv, method = "combat", batch = "batch", newAssayName = "naive"
    )
    bv <- runCorrection(
        bv, method = "combat", batch = "batch", preserve = "group",
        newAssayName = "protected"
    )

    out <- paste(capture.output(summary(bv)), collapse = "\n")

    ## without a preserve column the two rows are identical but for names
    expect_match(out, "preserve")
    expect_match(out, "group")
})


test_that("batch is limited to a single variable", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 30, nSamples = 8)

    expect_error(
        runCorrection(bv, method = "combat", batch = c("batch", "group")),
        "single column name"
    )
})
