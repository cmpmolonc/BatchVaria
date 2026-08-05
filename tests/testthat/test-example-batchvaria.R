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
