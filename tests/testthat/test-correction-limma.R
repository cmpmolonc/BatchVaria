## The two built-in methods estimate different things -- ComBat shrinks each
## batch's location and scale, removeBatchEffect removes a location shift --
## so these assert direction and never magnitude. A tolerance loose enough to
## admit both would test nothing, and a tight one would break on any sva or
## limma update.

skip_if_no_limma <- function() {
    testthat::skip_if_not_installed("limma")
}


test_that("limma is registered and returns a conforming matrix", {
    skip_if_no_limma()

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100, nSamples = 12)
    raw <- as.matrix(SummarizedExperiment::assay(bv, "raw"))

    out <- runCorrection(
        bv,
        method = "limma", batch = "batch", preserve = "group",
        newAssayName = "limma_corrected"
    )
    corrected <- as.matrix(
        SummarizedExperiment::assay(out, "limma_corrected")
    )

    expect_true(is.numeric(corrected))
    expect_identical(dimnames(corrected), dimnames(raw))
    expect_false(identical(corrected, raw))
})


test_that("both methods reduce batch variance and retain group", {
    skip_if_no_limma()

    ## balanced design: batch and group are orthogonal, so both are
    ## separately estimable and neither method is asked to do the impossible
    set.seed(42)
    bv <- exampleBatchVaria(nGenes = 400, nSamples = 20, confounding = 0)

    bv <- runCorrection(
        bv, method = "combat", batch = "batch", preserve = "group",
        newAssayName = "combat"
    )
    bv <- runCorrection(
        bv, method = "limma", batch = "batch", preserve = "group",
        newAssayName = "limma"
    )

    bv <- profileVariance(
        bv, ~ batch + group,
        assays = c("raw", "combat", "limma"), methods = "anova"
    )

    ## absolute variance, not fractions: fractions are compositional, so a
    ## rising group fraction would not by itself show group was retained
    res <- varianceResults(bv, method = "anova")
    tot <- assayVariance(bv, assays = c("raw", "combat", "limma"))

    absVar <- function(assayName, term) {
        frac <- res$variance_fraction[
            res$assay == assayName & res$term == term
        ]
        frac * tot$total_variance[tot$assay == assayName]
    }

    for (m in c("combat", "limma")) {
        expect_lt(absVar(m, "batch"), absVar("raw", "batch"))
        expect_gt(
            absVar(m, "group"),
            0.5 * absVar("raw", "group")
        )
    }

    ## the direction is the same for both, which is what a backwards
    ## 'preserve' mapping would break: wiring preserve to limma's
    ## 'covariates' would strip group instead of protecting it
    expect_gt(absVar("limma", "group"), absVar("limma", "batch"))
    expect_gt(absVar("combat", "group"), absVar("combat", "batch"))
})


test_that("preserving group changes what limma returns", {
    skip_if_no_limma()

    set.seed(7)
    bv <- exampleBatchVaria(nGenes = 200, nSamples = 16, confounding = 0.4)

    naive <- runCorrection(
        bv, method = "limma", batch = "batch", newAssayName = "out"
    )
    protectedRun <- runCorrection(
        bv, method = "limma", batch = "batch", preserve = "group",
        newAssayName = "out"
    )

    expect_false(identical(
        as.matrix(SummarizedExperiment::assay(naive, "out")),
        as.matrix(SummarizedExperiment::assay(protectedRun, "out"))
    ))
})


test_that("removeBatchEffect stays on the reference basis", {
    skip_if_no_limma()

    ## removeBatchEffect is unambiguously a per-feature linear shift, so it
    ## should retain the original axes as closely as ComBat does. This is an
    ## independent check on basisRetention() separating structure-preserving
    ## corrections from distortions, using a method whose linearity is not
    ## in question.
    set.seed(3)
    bv <- exampleBatchVaria(nGenes = 400, nSamples = 20, confounding = 0.6)

    bv <- runCorrection(
        bv, method = "combat", batch = "batch", newAssayName = "combat"
    )
    bv <- runCorrection(
        bv, method = "limma", batch = "batch", newAssayName = "limma"
    )

    ret <- basisRetention(
        bv,
        assays = c("combat", "limma", "raw_noise"), reference = "raw"
    )

    getRet <- function(a) ret$retention[ret$assay == a]

    expect_gt(getRet("limma"), 0.95)
    expect_gt(getRet("combat"), 0.95)

    ## and the distorted control is clearly separated from both
    expect_lt(getRet("raw_noise"), getRet("limma"))
    expect_lt(getRet("raw_noise"), getRet("combat"))
})
