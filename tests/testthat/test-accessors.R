test_that("availableCorrectionMethods lists the correction engines", {
    m <- availableCorrectionMethods()
    expect_type(m, "character")
    expect_true("combat" %in% m)
    expect_false(anyDuplicated(m) > 0)
})


test_that("availableVarianceMethods lists the variance engines", {
    m <- availableVarianceMethods()
    expect_type(m, "character")
    expect_setequal(m, c("anova", "variancePartition"))

    ## PCA is an embedding, not a variance engine
    expect_false("pca" %in% m)
})


test_that("varianceHistory returns the ledger", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    ## empty before any profiling
    expect_null(varianceHistory(bv))

    bv <- profileVariance(bv, ~batch, assays = "raw", methods = "anova")
    vh <- varianceHistory(bv)

    expect_type(vh, "list")
    expect_length(vh, 1)
    expect_true(all(c("assay", "formula", "method", "result", "timestamp") %in%
        names(vh[[1]])))
    expect_equal(vh[[1]]$assay, "raw")
    expect_equal(vh[[1]]$method, "anova")

    ## grows with each profiling call
    bv <- profileVariance(bv, ~batch, assays = "raw_center", methods = "anova")
    expect_length(varianceHistory(bv), 2)
})


test_that("varianceResults filters by assay and method", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)
    bv <- profileVariance(
        bv, ~batch,
        assays = c("raw", "raw_center"),
        methods = c("anova", "variancePartition")
    )

    all_res <- varianceResults(bv)
    expect_setequal(unique(all_res$assay), c("raw", "raw_center"))
    expect_setequal(unique(all_res$method), c("anova", "variancePartition"))

    expect_setequal(unique(varianceResults(bv, assayName = "raw")$assay), "raw")
    expect_setequal(
        unique(varianceResults(bv, method = "variancePartition")$method),
        "variancePartition"
    )

    both <- varianceResults(bv, assayName = "raw", method = "anova")
    expect_setequal(unique(both$assay), "raw")
    expect_setequal(unique(both$method), "anova")
})


test_that("assayVariance returns one row per assay with the expected schema", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)

    av <- assayVariance(bv)

    expect_true(all(
        c("assay", "total_variance", "n_features", "n_samples") %in%
            colnames(av)
    ))
    expect_setequal(av$assay, SummarizedExperiment::assayNames(bv))
    expect_true(is.numeric(av$total_variance))
    expect_true(all(av$total_variance > 0))
    expect_true(all(av$n_features == 100))
    expect_true(all(av$n_samples == ncol(bv)))
})


test_that("assayVariance equals the sum of per-feature variances", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)

    mat <- as.matrix(SummarizedExperiment::assay(bv, "raw"))

    expect_equal(
        assayVariance(bv, assays = "raw")$total_variance,
        sum(apply(mat, 1, stats::var))
    )
})


test_that("assayVariance tracks known variance-preserving transforms", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)

    av <- assayVariance(bv)
    tv <- stats::setNames(av$total_variance, av$assay)

    ## per-feature centring shifts means only, so variance is unchanged
    expect_equal(tv[["raw_center"]], tv[["raw"]])

    ## raw_scale is raw * 0.9, so variance scales by 0.9^2
    expect_equal(tv[["raw_scale"]], tv[["raw"]] * 0.81)

    ## added noise increases total variance
    expect_gt(tv[["raw_noise"]], tv[["raw"]])
})


test_that("assayVariance validates assay names", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)

    expect_error(assayVariance(bv, assays = "nope"), "Assays not found")
    expect_error(assayVariance(list()), "is\\(object, \"BatchVariaData\"\\)")
})


test_that("runCorrection records total variance in the correction ledger", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)
    bv <- runCorrection(bv, method = "combat", batch = "batch")

    entry <- S4Vectors::metadata(bv)$correction_history[[1]]

    expect_true(all(
        c("total_variance_in", "total_variance_out") %in% names(entry)
    ))

    av <- assayVariance(bv, assays = c("raw", "raw_combat"))
    expect_equal(entry$total_variance_in, av$total_variance[av$assay == "raw"])
    expect_equal(
        entry$total_variance_out,
        av$total_variance[av$assay == "raw_combat"]
    )
})


test_that("a rising variance fraction can mask falling absolute variance", {
    ## This is the compositional artefact assayVariance() exists to expose:
    ## ComBat raises the group *fraction* while lowering the group
    ## *absolute* variance. If this ever stops holding, the vignette's
    ## interpretation section needs revisiting.
    skip_if_not_installed("variancePartition")

    set.seed(1)
    bv <- exampleBatchVaria()
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    bv <- profileVariance(
        bv,
        formula = ~ batch + group,
        assays = c("raw", "raw_combat"),
        methods = "variancePartition"
    )

    vr <- varianceResults(bv, method = "variancePartition")
    av <- assayVariance(bv, assays = c("raw", "raw_combat"))

    frac <- function(a) {
        vr$variance_fraction[vr$assay == a & vr$term == "group"]
    }
    total <- function(a) av$total_variance[av$assay == a]

    ## fraction rises
    expect_gt(frac("raw_combat"), frac("raw"))

    ## total variance falls
    expect_lt(total("raw_combat"), total("raw"))

    ## and the absolute group variance falls despite the rising fraction
    expect_lt(
        frac("raw_combat") * total("raw_combat"),
        frac("raw") * total("raw")
    )
})
