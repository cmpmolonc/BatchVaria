## The harmonisation claim is that engines are extensible, not that two of
## them happen to exist. These tests register an engine defined entirely
## outside the package and assert it is indistinguishable from a built-in
## everywhere downstream.

## A third-party engine: no BatchVaria internals, only the public
## constructor and the documented signature.
dummyEngine <- function(assayMatrix, formula, sampleData, ...) {
    newVarianceSummary(
        source = "dummy",
        term = c("batch", "residual"),
        varianceFraction = c(0.25, 0.75),
        nFeatures = nrow(assayMatrix)
    )
}

withDummyEngine <- function(code) {
    registerVarianceEngine("dummy", dummyEngine)
    on.exit(unregisterVarianceEngine("dummy"), add = TRUE)
    force(code)
}


test_that("a registered engine appears alongside the built-ins", {
    builtins <- availableVarianceMethods()
    expect_false("dummy" %in% builtins)

    withDummyEngine({
        expect_setequal(
            availableVarianceMethods(),
            c(builtins, "dummy")
        )
    })

    ## and is gone again afterwards
    expect_setequal(availableVarianceMethods(), builtins)
})


test_that("a registered engine flows through the whole pipeline", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    withDummyEngine({
        res <- profileVariance(
            bv, ~batch,
            assays = c("raw", "raw_noise"), methods = "dummy"
        )

        ## recordVariance(): keyed like any other engine
        entry <- varianceHistory(res)[[1]]
        expect_equal(entry$method, "dummy")
        expect_equal(entry$formula_key, "~batch")

        ## varianceResults(): same schema
        vr <- varianceResults(res, method = "dummy")
        expect_true(all(
            c("source", "term", "variance_fraction", "n_features",
                "assay", "method", "formula") %in% colnames(vr)
        ))
        expect_setequal(unique(vr$assay), c("raw", "raw_noise"))

        ## varianceTable(): summarised without special-casing
        vt <- varianceTable(
            res,
            assays = c("raw", "raw_noise"), method = "dummy"
        )
        expect_setequal(vt$percent$component, c("batch", "residual"))
        expect_equal(sum(vt$percent$raw), 100)
        expect_equal(vt$baseline, "raw")

        ## varianceDelta() and the plotting layer accept it too
        vd <- varianceDelta(
            res,
            assays = c("raw", "raw_noise"), method = "dummy"
        )
        expect_true(all(is.finite(vd$delta)))
        expect_s3_class(
            plotVarianceComposition(res, method = "dummy"),
            "ggplot"
        )
    })
})


test_that("a registered engine is subject to the same ledger rules", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    withDummyEngine({
        res <- profileVariance(bv, ~batch, assays = "raw", methods = "dummy")

        ## mixing it with a built-in is refused, exactly as two built-ins are
        res <- profileVariance(res, ~batch, assays = "raw", methods = "anova")
        expect_error(varianceTable(res, assays = "raw"), "holds results for 2 methods")

        ## re-profiling the same key warns
        expect_warning(
            profileVariance(res, ~batch, assays = "raw", methods = "dummy"),
            "already recorded"
        )
    })
})


test_that("registration enforces the engine contract", {
    ## missing '...': profileVariance() forwards engine-specific arguments
    expect_error(
        registerVarianceEngine(
            "bad", function(assayMatrix, formula, sampleData) NULL
        ),
        "must accept '\\.\\.\\.'"
    )

    ## missing a required argument
    expect_error(
        registerVarianceEngine("bad", function(assayMatrix, ...) NULL),
        "is missing formula, sampleData"
    )

    expect_error(registerVarianceEngine("bad", "not a function"), "must be a function")
    expect_error(registerVarianceEngine(character(0), dummyEngine), "single non-empty string")

    expect_false("bad" %in% availableVarianceMethods())
})


test_that("registration refuses to clobber silently", {
    withDummyEngine({
        expect_error(
            registerVarianceEngine("dummy", dummyEngine),
            "already registered"
        )

        other <- function(assayMatrix, formula, sampleData, ...) {
            newVarianceSummary(
                source = "dummy", term = "all",
                varianceFraction = 1, nFeatures = nrow(assayMatrix)
            )
        }
        expect_silent(
            registerVarianceEngine("dummy", other, overwrite = TRUE)
        )
    })

    expect_error(unregisterVarianceEngine("dummy"), "No variance engine named")
})


test_that("an engine breaking the result contract is rejected, not recorded", {
    ## The contract is enforced at use, not only at registration: a
    ## conforming signature does not guarantee conforming output.
    badOutput <- function(assayMatrix, formula, sampleData, ...) {
        data.frame(nonsense = 1:2)
    }

    registerVarianceEngine("badOutput", badOutput)
    on.exit(unregisterVarianceEngine("badOutput"), add = TRUE)

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_error(
        suppressWarnings(
            profileVariance(bv, ~batch, assays = "raw", methods = "badOutput")
        ),
        "All 1 variance profiling attempt failed"
    )
})


test_that("newVarianceSummary requires the feature count", {
    expect_error(
        newVarianceSummary(
            source = "x", term = "a", varianceFraction = 1
        ),
        "'nFeatures' must be supplied"
    )

    ok <- newVarianceSummary(
        source = "x", term = c("a", "b"),
        varianceFraction = c(0.4, 0.6), nFeatures = 12
    )
    expect_equal(ok$n_features, c(12, 12))
})
