## Every cross-assay PCA view shares one basis, fitted on a reference
## assay. Independent per-assay fits make "PC1" a different direction in
## each row; a joint fit over the pair being plotted makes the same samples
## move depending on what they are compared against.

pcaFixture <- function(nGenes = 300) {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = nGenes)
    suppressMessages(runCorrection(bv, method = "combat", batch = "batch"))
}


test_that("the reference assay's own decomposition is unchanged", {
    bv <- pcaFixture()

    res <- comparePCA(bv, assays = c("raw", "raw_combat"), nPCs = 3)
    ours <- res$variance[res$assay == "raw"]

    ## projecting the reference onto its own basis is ordinary PCA
    p <- stats::prcomp(
        t(as.matrix(SummarizedExperiment::assay(bv, "raw"))),
        scale. = TRUE
    )
    theirs <- (p$sdev^2 / sum(p$sdev^2))[seq_len(3)]

    expect_equal(ours, theirs, tolerance = 1e-8)
})


test_that("fractions sum to one for the reference and less for others", {
    bv <- pcaFixture()
    basis <- BatchVaria:::.referenceBasis(bv, "raw")

    totalOf <- function(a) {
        sum(BatchVaria:::.basisVarianceFractions(bv, a, basis))
    }

    expect_equal(totalOf("raw"), 1, tolerance = 1e-8)

    ## noise adds structure the raw basis does not span, so much of its
    ## variance falls outside the projection
    expect_lt(totalOf("raw_noise"), 0.9)

    ## ComBat reshapes variance within the existing structure
    expect_gt(totalOf("raw_combat"), 0.9)
    expect_lte(totalOf("raw_combat"), 1 + 1e-8)
})


test_that("comparePCA reports the reference it used", {
    bv <- pcaFixture()

    res <- comparePCA(bv, assays = c("raw", "raw_combat"))

    ## resolved from the correction lineage, not hard-coded
    expect_setequal(unique(res$reference), "raw")

    explicit <- comparePCA(
        bv,
        assays = c("raw", "raw_combat"),
        reference = "raw_combat"
    )
    expect_setequal(unique(explicit$reference), "raw_combat")

    ## the nominated reference is the one that sums to one
    expect_equal(
        sum(BatchVaria:::.basisVarianceFractions(
            bv, "raw_combat", BatchVaria:::.referenceBasis(bv, "raw_combat")
        )),
        1,
        tolerance = 1e-8
    )
})


test_that("trajectory before-positions do not depend on the comparator", {
    ## The joint-fit version placed the same raw samples differently in a
    ## raw -> combat plot and a raw -> noise plot, so two trajectories from
    ## one object could not be overlaid.
    bv <- pcaFixture()

    coords <- function(after) {
        p <- plotPCATrajectory(
            bv,
            assayBefore = "raw", assayAfter = after, colourBy = "batch"
        )
        p$data[, c("x1", "y1")]
    }

    expect_equal(coords("raw_combat"), coords("raw_noise"))
})


test_that("plotBasisProjection panels share a basis across assays", {
    bv <- pcaFixture()

    plots <- plotBasisProjection(
        bv,
        assays = c("raw", "raw_combat"), colourBy = "batch"
    )

    expect_length(plots, 2)

    basis <- BatchVaria:::.referenceBasis(bv, "raw")
    expected <- BatchVaria:::.projectOntoBasis(bv, "raw_combat", basis)

    expect_equal(plots[[2]]$data$PC1, unname(expected[, 1]))
})


test_that("the basis tolerates constant features and rejects empty ones", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)
    mat <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    mat[3, ] <- 5

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat),
        colData = SummarizedExperiment::colData(bv)
    )
    obj <- se

    ## A constant feature is kept in the basis on its original scale
    ## rather than dropped. Dropping it would remove it from the coordinate
    ## system, so variance a correction later introduced there would vanish
    ## from the measurement instead of counting as off-basis.
    basis <- BatchVaria:::.referenceBasis(obj, "raw")
    expect_length(basis$features, 60)
    expect_equal(basis$flat, rownames(mat)[3])
    expect_equal(basis$scale[[3]], 1)

    dead <- mat
    dead[] <- 3
    seDead <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = dead),
        colData = SummarizedExperiment::colData(bv)
    )
    expect_error(
        BatchVaria:::.referenceBasis(seDead, "raw"),
        "fewer than two features"
    )
})


test_that("variance introduced in a reference-flat feature is not lost", {
    ## A correction that puts signal into a previously flat feature is
    ## exactly what a distorting correction does. Excluding such features
    ## from the basis made that variance vanish from the measurement -
    ## retention came back as 1.0 - rather than counting as off-basis.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 100)
    mat <- as.matrix(SummarizedExperiment::assay(bv, "raw"))
    flat <- rownames(mat)[c(3, 7, 11)]
    mat[flat, ] <- 5

    revived <- mat
    revived[flat, ] <- stats::rnorm(length(flat) * ncol(mat), sd = 6)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat, revived = revived),
        colData = SummarizedExperiment::colData(bv)
    )
    obj <- se

    res <- basisRetention(obj, reference = "raw")

    expect_equal(res$retention[res$assay == "raw"], 1)
    expect_lt(res$retention[res$assay == "revived"], 0.9)
    expect_gt(res$off_basis[res$assay == "revived"], 0.1)
})


test_that("basisRetention separates reshaped from newly introduced variance", {
    bv <- pcaFixture(nGenes = 200)

    res <- basisRetention(
        bv,
        assays = c("raw", "raw_combat", "raw_noise")
    )

    expect_true(all(
        c("assay", "reference", "retention", "off_basis", "n_features") %in%
            colnames(res)
    ))
    expect_setequal(unique(res$reference), "raw")

    ret <- stats::setNames(res$retention, res$assay)

    ## the reference is exactly one, not one give or take rounding
    expect_identical(ret[["raw"]], 1)

    ## ComBat redistributes variance within the existing structure
    expect_gt(ret[["raw_combat"]], 0.95)

    ## added noise occupies directions the raw data did not
    expect_lt(ret[["raw_noise"]], 0.6)

    expect_equal(res$retention + res$off_basis, rep(1, nrow(res)))
})
