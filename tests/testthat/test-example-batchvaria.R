test_that("exampleBatchVaria returns valid object", {
    bv <- exampleBatchVaria()
    expect_s4_class(bv, "BatchVariaData")
    expect_true("raw" %in% SummarizedExperiment::assayNames(bv))
    expect_true(ncol(SummarizedExperiment::assay(bv, "raw")) > 0)
})


test_that("exampleBatchVaria builds a balanced 2 x 2 design", {
    for (n in c(4, 8, 12, 20)) {
        bv <- exampleBatchVaria(nGenes = 20, nSamples = n)
        cd <- SummarizedExperiment::colData(bv)

        expect_equal(ncol(bv), n)
        ## both factors split evenly, and the two are crossed not nested
        expect_equal(unname(table(cd$batch)), unname(table(cd$group)))
        expect_true(all(table(cd$batch) == n / 2))
        expect_true(all(table(cd$batch, cd$group) == n / 4))
    }
})


test_that("exampleBatchVaria rejects sample counts it cannot balance", {
    ## previously failed deep inside DataFrame() with
    ## "different row counts implied by arguments"
    for (n in c(1, 2, 3, 5, 6, 7, 10, 14, 21)) {
        expect_error(
            exampleBatchVaria(nGenes = 20, nSamples = n),
            "must be a multiple of 4"
        )
    }

    ## message points at the nearest usable values
    expect_error(
        exampleBatchVaria(nGenes = 20, nSamples = 10),
        "8 or 12"
    )

    expect_error(exampleBatchVaria(nSamples = "8"), "single number")
    expect_error(exampleBatchVaria(nSamples = c(4, 8)), "single number")
})


test_that("exampleBatchVaria is reproducible under a caller-set seed", {
    set.seed(1)
    a <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw")
    set.seed(1)
    b <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw")
    expect_identical(a, b)

    ## the random 'raw_noise' assay must be reproducible too
    set.seed(2)
    n1 <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw_noise")
    set.seed(2)
    n2 <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw_noise")
    expect_identical(n1, n2)

    ## a different seed gives different data
    set.seed(99)
    c1 <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw")
    expect_false(identical(a, c1))
})


test_that("exampleBatchVaria does not seed the RNG itself", {
    ## seeding is the caller's job, so two unseeded calls must differ and the
    ## stream must keep advancing rather than being reset inside the function
    x <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw")
    y <- SummarizedExperiment::assay(exampleBatchVaria(nGenes = 40, nSamples = 8), "raw")
    expect_false(identical(x, y))
})


test_that("exampleBatchVaria rejects gene counts with empty effect blocks", {
    ## round(0.10 * nGenes) is 0 below 10 genes, leaving no 'both' block
    expect_error(exampleBatchVaria(nGenes = 5, nSamples = 8), "at least 10")
    expect_error(exampleBatchVaria(nGenes = NA, nSamples = 8), "single number")
})


test_that("confounding slides the design from orthogonal to aliased", {
    crossTab <- function(cf) {
        cd <- SummarizedExperiment::colData(
            exampleBatchVaria(nGenes = 40, nSamples = 24, confounding = cf)
        )
        unname(as.vector(t(table(cd$batch, cd$group))))
    }

    expect_equal(crossTab(0), c(6, 6, 6, 6))
    expect_equal(crossTab(0.5), c(9, 3, 3, 9))
    expect_equal(crossTab(1), c(12, 0, 0, 12))
})


test_that("the default design is unchanged by the new argument", {
    set.seed(1)
    a <- exampleBatchVaria(nGenes = 30)
    set.seed(1)
    b <- exampleBatchVaria(nGenes = 30, confounding = 0)

    expect_identical(
        SummarizedExperiment::colData(a)$group,
        SummarizedExperiment::colData(b)$group
    )
    expect_equal(
        as.matrix(SummarizedExperiment::assay(a, "raw")),
        as.matrix(SummarizedExperiment::assay(b, "raw"))
    )
})


test_that("confounding is validated", {
    expect_error(
        exampleBatchVaria(nGenes = 20, confounding = 1.5),
        "between 0 and 1"
    )
    expect_error(
        exampleBatchVaria(nGenes = 20, confounding = -0.1),
        "between 0 and 1"
    )
    expect_error(
        exampleBatchVaria(nGenes = 20, confounding = "lots"),
        "between 0 and 1"
    )
})


test_that("confounded data makes the shared term non-zero", {
    ## The whole point of the argument: on an orthogonal design every unit
    ## of variance belongs to one term or the other, so 'shared' is zero
    ## and the confounding story cannot be told with this data.
    sharedOf <- function(cf) {
        set.seed(1)
        bv <- exampleBatchVaria(nGenes = 300, confounding = cf)
        bv <- suppressWarnings(suppressMessages(
            profileVariance(bv, ~ batch + group, assays = "raw", methods = "anova")
        ))
        res <- varianceResults(bv, assayName = "raw")
        res$variance_fraction[res$term == "shared"]
    }

    expect_equal(sharedOf(0), 0, tolerance = 1e-8)
    expect_gt(sharedOf(0.6), 0.05)
    expect_gt(sharedOf(0.9), sharedOf(0.6))
})


test_that("the engines diverge on confounded data and agree without it", {
    ## variancePartition has no way to express variance belonging to
    ## neither term, so it attributes the ambiguous share to both.
    skip_if_not_installed("variancePartition")

    terms <- function(cf) {
        set.seed(1)
        bv <- exampleBatchVaria(nGenes = 300, confounding = cf)
        bv <- suppressWarnings(suppressMessages(profileVariance(
            bv, ~ batch + group,
            assays = "raw", methods = c("anova", "variancePartition")
        )))
        vr <- varianceResults(bv, assayName = "raw")
        c(
            anova = vr$variance_fraction[
                vr$method == "anova" & vr$term == "batch"
            ],
            vp = vr$variance_fraction[
                vr$method == "variancePartition" & vr$term == "batch"
            ]
        )
    }

    balanced <- terms(0)
    expect_equal(balanced[["anova"]], balanced[["vp"]], tolerance = 1e-3)

    confounded <- terms(0.8)
    expect_gt(confounded[["vp"]], confounded[["anova"]] * 1.2)
})
