test_that("runCorrection creates new assay", {
    bv <- exampleBatchVaria()

    bv <- runCorrection(bv, method = "combat", batch = "batch")

    assays <- SummarizedExperiment::assayNames(bv)

    expect_true(any(grepl("combat", assays)))
})


test_that("runCorrection accepts covariates", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)

    ## exercises the model.matrix() path passed to ComBat as 'mod'
    out <- runCorrection(
        bv,
        method = "combat", batch = "batch", covariates = "group"
    )

    expect_s4_class(out, "BatchVariaData")
    expect_true("raw_combat" %in% SummarizedExperiment::assayNames(out))
    expect_equal(dim(SummarizedExperiment::assay(out, "raw_combat")), c(60L, 8L))

    ## the covariate is recorded in the correction ledger
    ch <- S4Vectors::metadata(out)$correction_history
    expect_equal(ch[[length(ch)]]$covariates, "group")

    ## protecting a covariate should change the result
    plain <- runCorrection(bv, method = "combat", batch = "batch")
    expect_false(identical(
        SummarizedExperiment::assay(out, "raw_combat"),
        SummarizedExperiment::assay(plain, "raw_combat")
    ))
})


test_that("runCorrection rejects unknown covariates", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    expect_error(
        runCorrection(bv, method = "combat", batch = "batch", covariates = "nope"),
        "Covariates not found in colData: nope"
    )
    expect_error(
        runCorrection(
            bv,
            method = "combat", batch = "batch",
            covariates = c("group", "nope")
        ),
        "Covariates not found in colData: nope"
    )
})
