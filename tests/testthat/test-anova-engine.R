## The anova engine decomposes variance per model term using Type II sums
## of squares, with an explicit 'shared' row for variance an unbalanced
## design cannot attribute to any single term.

makeDesign <- function(batch, group, nGenes = 300, seed = 1) {
    set.seed(seed)
    cd <- S4Vectors::DataFrame(
        batch = batch, group = group,
        row.names = paste0("s", seq_along(batch))
    )
    mat <- matrix(stats::rnorm(nGenes * length(batch)), nGenes)
    mat <- mat +
        outer(stats::rnorm(nGenes, 0, 1.5), as.numeric(batch == "B")) +
        outer(stats::rnorm(nGenes, 0, 1.5), as.numeric(group == "Y"))
    dimnames(mat) <- list(paste0("g", seq_len(nGenes)), rownames(cd))
    list(mat = mat, colData = cd)
}

fractions <- function(res) {
    stats::setNames(res$variance_fraction, res$term)
}


test_that("anova reports one row per model term plus shared and residual", {
    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))

    res <- BatchVaria:::.computeAnovaVariance(d$mat, ~ batch + group, d$colData)

    expect_equal(res$term, c("batch", "group", "shared", "residual"))
    expect_equal(sum(res$variance_fraction), 1, tolerance = 1e-8)
})


test_that("the composition is a partition for one-term and interaction models", {
    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))

    one <- BatchVaria:::.computeAnovaVariance(d$mat, ~batch, d$colData)
    expect_equal(one$term, c("batch", "shared", "residual"))
    expect_equal(sum(one$variance_fraction), 1, tolerance = 1e-8)

    inter <- BatchVaria:::.computeAnovaVariance(
        d$mat, ~ batch * group, d$colData
    )
    expect_true("batch:group" %in% inter$term)
    expect_equal(sum(inter$variance_fraction), 1, tolerance = 1e-8)
})


test_that("shared is zero when balanced and grows with confounding", {
    balanced <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))
    partial <- makeDesign(
        rep(c("A", "B"), each = 12),
        c(rep("X", 11), "Y", "X", rep("Y", 11))
    )
    complete <- makeDesign(
        rep(c("A", "B"), each = 12),
        rep(c("X", "Y"), each = 12)
    )

    sharedOf <- function(d) {
        fractions(
            BatchVaria:::.computeAnovaVariance(d$mat, ~ batch + group, d$colData)
        )[["shared"]]
    }

    ## an orthogonal design leaves nothing unattributable
    expect_equal(sharedOf(balanced), 0, tolerance = 1e-8)

    expect_gt(sharedOf(partial), 0.05)
    expect_gt(sharedOf(complete), sharedOf(partial))
})


test_that("complete confounding attributes nothing to either term", {
    ## When batch and group are perfectly aliased neither is identifiable
    ## given the other, so Type II assigns each zero and the whole of the
    ## explained variance falls to 'shared'. The old lumped model/residual
    ## split reported this as a large 'model' fraction with no indication
    ## that none of it was attributable.
    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), each = 12))

    f <- fractions(
        BatchVaria:::.computeAnovaVariance(d$mat, ~ batch + group, d$colData)
    )

    expect_equal(f[["batch"]], 0, tolerance = 1e-8)
    expect_equal(f[["group"]], 0, tolerance = 1e-8)
    expect_gt(f[["shared"]], 0.1)
})


test_that("anova agrees with variancePartition on a balanced design", {
    ## Independent implementations of the same question should converge.
    skip_if_not_installed("variancePartition")

    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))

    ours <- fractions(
        BatchVaria:::.computeAnovaVariance(d$mat, ~ batch + group, d$colData)
    )
    theirs <- colMeans(suppressMessages(
        variancePartition::fitExtractVarPartModel(
            exprObj = d$mat, formula = ~ batch + group,
            data = as.data.frame(d$colData)
        )
    ))

    expect_equal(ours[["batch"]], theirs[["batch"]], tolerance = 1e-3)
    expect_equal(ours[["group"]], theirs[["group"]], tolerance = 1e-3)
})


test_that("term order does not change the decomposition", {
    ## Type II sums of squares are order independent, unlike Type I.
    d <- makeDesign(
        rep(c("A", "B"), each = 12),
        c(rep("X", 8), rep("Y", 4), rep("X", 4), rep("Y", 8))
    )

    a <- fractions(
        BatchVaria:::.computeAnovaVariance(d$mat, ~ batch + group, d$colData)
    )
    b <- fractions(
        BatchVaria:::.computeAnovaVariance(d$mat, ~ group + batch, d$colData)
    )

    expect_equal(a[["batch"]], b[["batch"]])
    expect_equal(a[["group"]], b[["group"]])
    expect_equal(a[["shared"]], b[["shared"]])
})


test_that("the weighting choice is exposed and recorded", {
    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))

    byFeature <- BatchVaria:::.computeAnovaVariance(
        d$mat, ~ batch + group, d$colData,
        weighting = "feature"
    )
    pooled <- BatchVaria:::.computeAnovaVariance(
        d$mat, ~ batch + group, d$colData,
        weighting = "pooled"
    )

    expect_equal(unique(byFeature$weighting), "feature")
    expect_equal(unique(pooled$weighting), "pooled")
    expect_equal(unique(byFeature$ss_type), "II")

    ## both remain partitions, but weight features differently
    expect_equal(sum(pooled$variance_fraction), 1, tolerance = 1e-8)
    expect_false(isTRUE(all.equal(
        byFeature$variance_fraction, pooled$variance_fraction
    )))
})


test_that("Type III is refused rather than silently approximated", {
    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))

    expect_error(
        BatchVaria:::.computeAnovaVariance(
            d$mat, ~ batch + group, d$colData,
            ssType = "III"
        ),
        "Only Type II sums of squares are supported"
    )
})


test_that("no-intercept designs still yield a partition", {
    d <- makeDesign(rep(c("A", "B"), each = 12), rep(c("X", "Y"), 12))

    res <- BatchVaria:::.computeAnovaVariance(
        d$mat, ~ 0 + batch + group, d$colData
    )

    expect_equal(sum(res$variance_fraction), 1, tolerance = 1e-8)
    expect_true(all(res$variance_fraction[res$term != "shared"] >= 0))
})


test_that("varianceResults tolerates engine-specific columns", {
    ## The anova engine carries ss_type and weighting; the others do not.
    skip_if_not_installed("variancePartition")

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 80)
    bv <- suppressWarnings(suppressMessages(
        profileVariance(bv, ~batch, assays = "raw")
    ))

    res <- varianceResults(bv)

    expect_true(all(c("ss_type", "weighting") %in% colnames(res)))
    expect_true(any(is.na(res$ss_type)))
    expect_setequal(
        stats::na.omit(unique(res$ss_type)),
        "II"
    )
})


test_that("the two covariate engines agree on a balanced design", {
    ## The vignette's concordance section rests on this. If it stops
    ## holding, that text needs revisiting rather than quietly going stale.
    skip_if_not_installed("variancePartition")

    set.seed(1)
    bv <- exampleBatchVaria()
    bv <- runCorrection(bv, method = "combat", batch = "batch")
    assaysUsed <- c("raw", "raw_combat", "raw_noise")

    bv <- suppressWarnings(suppressMessages(profileVariance(
        bv, ~ batch + group,
        assays = assaysUsed, methods = c("anova", "variancePartition")
    )))

    vr <- varianceResults(bv)

    for (a in assaysUsed) {
        pick <- function(m, term) {
            vr$variance_fraction[
                vr$assay == a & vr$method == m & vr$term == term
            ]
        }
        expect_equal(pick("anova", "batch"), pick("variancePartition", "batch"),
            tolerance = 1e-3, info = a)
        expect_equal(pick("anova", "group"), pick("variancePartition", "group"),
            tolerance = 1e-3, info = a)

        ## orthogonal design leaves nothing unattributable
        expect_equal(pick("anova", "shared"), 0, tolerance = 1e-6, info = a)
    }
})


test_that("plotVarianceComposition warns on negative segments", {
    ## A stacked bar cannot honestly render a negative share, so the chart
    ## must say so rather than silently drawing below the axis.
    set.seed(3)
    n <- 40
    g <- 150
    x1 <- stats::rnorm(n)
    x2 <- 0.9 * x1 + sqrt(1 - 0.81) * stats::rnorm(n)
    cd <- S4Vectors::DataFrame(a = x1, b = x2, row.names = paste0("s", seq_len(n)))

    mat <- matrix(stats::rnorm(g * n, sd = 0.4), g, n) +
        outer(stats::rnorm(g, 0, 2), x1 - x2)
    dimnames(mat) <- list(paste0("g", seq_len(g)), rownames(cd))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = mat), colData = cd
    )
    obj <- suppressWarnings(
        profileVariance(BatchVariaData(se), ~ a + b, methods = "anova")
    )

    res <- varianceResults(obj)
    expect_lt(res$variance_fraction[res$term == "shared"], 0)

    ## the partition still holds despite the negative component
    expect_equal(sum(res$variance_fraction), 1, tolerance = 1e-8)

    expect_warning(
        plotVarianceComposition(obj, method = "anova"),
        "Negative variance percentages"
    )
})
