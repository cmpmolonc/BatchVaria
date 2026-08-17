## What 'verbose' governs, and what it must not.
##
## The distinction is the point of these tests: engine progress reporting
## is suppressible, everything BatchVaria derives about the data is not.
## Kept in one file so the boundary stays visible.

test_that("engine progress is suppressed by default and shown on request", {
    chatty <- function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
        message("engine progress: standardizing")
        assayMatrix + 1
    }

    registerCorrectionMethod("chatty", chatty)
    on.exit(unregisterCorrectionMethod("chatty"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_no_message(
        runCorrection(bv, method = "chatty", batch = "batch")
    )

    expect_message(
        runCorrection(bv, method = "chatty", batch = "batch", verbose = TRUE),
        "engine progress"
    )
})


test_that("suppression does not reach warnings the engine raises", {
    ## An engine reporting a finding about the data rather than its own
    ## progress uses warning(), which verbose = FALSE must not touch.
    finder <- function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
        message("progress chatter")
        warning("3 features had uniform expression", call. = FALSE)
        assayMatrix + 1
    }

    registerCorrectionMethod("finder", finder)
    on.exit(unregisterCorrectionMethod("finder"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_warning(
        runCorrection(bv, method = "finder", batch = "batch"),
        "uniform expression"
    )
})


test_that("no-op detection survives suppression", {
    ## The detector compares matrices rather than reading what the method
    ## reported, so silencing a method cannot blind it. This is the one
    ## place suppression could have cost a safety signal.
    silentNoOp <- function(assayMatrix, batch, sampleData, preserve = NULL,
                           ...) {
        message("Coefficients not estimable: batch")
        assayMatrix
    }

    registerCorrectionMethod("silentNoOp", silentNoOp)
    on.exit(unregisterCorrectionMethod("silentNoOp"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_warning(
        runCorrection(bv, method = "silentNoOp", batch = "batch"),
        "indistinguishable"
    )

    bv2 <- suppressWarnings(
        runCorrection(bv, method = "silentNoOp", batch = "batch")
    )

    ch <- metadata(bv2)$correction_history
    expect_true(isTRUE(ch[[length(ch)]]$no_op))
})


test_that("profileVariance gates engine progress the same way", {
    chattyEngine <- function(assayMatrix, formula, sampleData, ...) {
        message("engine progress: fitting")
        newVarianceSummary(
            source = "chattyEngine",
            term = c("batch", "residual"),
            varianceFraction = c(0.4, 0.6),
            nFeatures = nrow(assayMatrix)
        )
    }

    registerVarianceEngine("chattyEngine", chattyEngine)
    on.exit(unregisterVarianceEngine("chattyEngine"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_no_message(
        profileVariance(bv, ~batch, assays = "raw", methods = "chattyEngine")
    )

    expect_message(
        profileVariance(
            bv, ~batch,
            assays = "raw", methods = "chattyEngine", verbose = TRUE
        ),
        "engine progress"
    )
})


test_that("engine failures are still reported when quiet", {
    ## A failed engine warns and is skipped. verbose governs progress, not
    ## whether a failure is announced.
    broken <- function(assayMatrix, formula, sampleData, ...) {
        message("engine progress: starting")
        stop("engine exploded")
    }

    registerVarianceEngine("broken", broken)
    on.exit(unregisterVarianceEngine("broken"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_warning(
        try(
            profileVariance(bv, ~batch, assays = "raw", methods = "broken"),
            silent = TRUE
        ),
        "engine exploded"
    )
})


test_that("verbose does not collide with arguments forwarded to engines", {
    ## 'verbose' sits after ... in both signatures, so an engine with its
    ## own 'verb'-prefixed argument still receives it rather than having it
    ## partially matched away by the wrapper's parameter.
    seen <- NULL

    picky <- function(assayMatrix, batch, sampleData, preserve = NULL,
                      verbosity = "none", ...) {
        seen <<- verbosity
        assayMatrix + 1
    }

    registerCorrectionMethod("picky", picky)
    on.exit(unregisterCorrectionMethod("picky"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    runCorrection(bv, method = "picky", batch = "batch", verbosity = "high")
    expect_equal(seen, "high")
})
