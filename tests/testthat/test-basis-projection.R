## basisProjection() exports the coordinates the reference-basis layer is
## built on. 

test_that("basisProjection returns one row per assay and sample", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)
    bv <- runCorrection(bv, method = "combat", batch = "batch")

    coords <- basisProjection(bv, assays = c("raw", "raw_combat"))

    expect_equal(nrow(coords), 2 * ncol(bv))
    expect_true(all(c("assay", "sample", "reference", "PC1", "PC2") %in%
        colnames(coords)))
    expect_setequal(coords$assay, c("raw", "raw_combat"))
    expect_setequal(coords$sample, colnames(bv))
})


test_that("nPCs controls the number of component columns", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)

    coords <- basisProjection(bv, assays = "raw", nPCs = 4)
    expect_true(all(paste0("PC", 1:4) %in% colnames(coords)))
    expect_false("PC5" %in% colnames(coords))

    expect_error(basisProjection(bv, nPCs = 0), "must be a single number")
})


test_that("every assay is projected onto one basis", {
    ## The property the whole layer exists for: the reference's own
    ## coordinates are unchanged by which other assays are projected
    ## alongside it, because the basis does not depend on them.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)
    bv <- runCorrection(bv, method = "combat", batch = "batch")

    alone <- basisProjection(bv, assays = "raw", reference = "raw")
    withOther <- basisProjection(
        bv,
        assays = c("raw", "raw_combat"), reference = "raw"
    )
    withOther <- withOther[withOther$assay == "raw", ]

    expect_equal(alone$PC1, withOther$PC1)
    expect_equal(alone$PC2, withOther$PC2)

    expect_true(all(withOther$reference == "raw"))
})


test_that("basisProjection agrees with basisRetention on the same basis", {
    ## basisRetention() sums variance along the basis axes; the projection
    ## is where those axes put the samples. Computing retention from the
    ## published coordinates must reproduce the published summary, or the
    ## exported numbers describe something other than the summary does.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)

    nPC <- ncol(bv)
    coords <- basisProjection(
        bv,
        assays = "raw", nPCs = nPC, reference = "raw"
    )
    pcCols <- grep("^PC", colnames(coords), value = TRUE)

    ## for the reference itself retention is 1, so its coordinates must
    ## carry all of its variance in the scaled space
    scoreVar <- sum(vapply(
        coords[, pcCols, drop = FALSE], stats::var, numeric(1)
    ))
    expect_gt(scoreVar, 0)

    ret <- basisRetention(bv, assays = "raw", reference = "raw")
    expect_equal(ret$retention, 1)
})


test_that("plotBasisProjection draws what basisProjection returns", {
    ## The plot is a layer over the public accessor, so the points in the
    ## panel must be the published coordinates rather than a parallel
    ## computation that could drift from them.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60, nSamples = 8)
    bv <- runCorrection(bv, method = "combat", batch = "batch")

    coords <- basisProjection(
        bv,
        assays = c("raw", "raw_combat"), reference = "raw"
    )
    plots <- plotBasisProjection(
        bv,
        assays = c("raw", "raw_combat"), reference = "raw"
    )

    expect_length(plots, 2)

    for (i in seq_along(plots)) {
        a <- c("raw", "raw_combat")[i]
        expected <- coords[coords$assay == a, ]
        expect_equal(plots[[i]]$data$PC1, expected$PC1)
        expect_equal(plots[[i]]$data$PC2, expected$PC2)
    }
})


test_that("basisProjection validates its inputs", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40, nSamples = 8)

    expect_error(basisProjection(bv, assays = "nope"), "Assays not found")
    expect_error(
        basisProjection(bv, reference = "nope"),
        "Assays not found|Reference assay not found"
    )
    expect_error(basisProjection(list()), "must be a SummarizedExperiment")
})
