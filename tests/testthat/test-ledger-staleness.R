## Staleness detection and S4 validity.
##
## These two mechanisms cover disjoint failures, and the division is not
## arbitrary: `[`, `assay<-`, `colData<-` and `metadata<-` do not invoke
## validity, so a validity method cannot see any mutation of an existing
## object. Validity therefore covers construction, and the fingerprint
## check covers everything afterwards. Tests for both are kept together so
## that the boundary stays visible.

makeProfiled <- function(nGenes = 40, nSamples = 20) {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = nGenes, nSamples = nSamples)
    bv <- suppressMessages(
        runCorrection(bv, method = "combat", batch = "batch")
    )
    profileVariance(bv, ~batch)
}

test_that("an unmodified object reports no staleness", {
    bv <- makeProfiled()

    expect_no_warning(varianceTable(bv))
    expect_no_warning(varianceResults(bv))
})

test_that("subsetting samples is detected at read time", {
    bv <- makeProfiled(nSamples = 20)

    ## The case that motivated the check: bv[, 1:6] is ordinary documented
    ## SummarizedExperiment behaviour, it does not invoke validity, and
    ## before this check varianceTable() returned percentages describing
    ## samples the object no longer held.
    expect_warning(
        varianceTable(bv[, 1:6]),
        "computed on 20 samples, object now has 6"
    )
})

test_that("subsetting features is detected at read time", {
    bv <- makeProfiled(nGenes = 40)

    expect_warning(
        varianceResults(bv[1:10, ]),
        "computed on 40 features, object now has 10"
    )
})

test_that("an assay replaced in place is detected despite equal dims", {
    bv <- makeProfiled()

    mat <- SummarizedExperiment::assay(bv, "raw")
    mat[1, 1] <- mat[1, 1] + 1e-9
    SummarizedExperiment::assay(bv, "raw") <- mat

    ## Dimensions are unchanged, so only a content-sensitive fingerprint
    ## can see this. The message says so rather than reporting a dimension
    ## mismatch that did not occur.
    expect_warning(
        varianceResults(bv, assayName = "raw"),
        "contents have changed.*dimensions have not"
    )
})

test_that("a stale table is still returned, not withheld", {
    bv <- makeProfiled()

    tab <- suppressWarnings(varianceTable(bv[, 1:6]))

    ## A stale result remains evidence of what was computed. Refusing to
    ## return it would also make bv[, 1:6] effectively illegal.
    expect_true(is.data.frame(tab$percent))
    expect_true("batch" %in% tab$percent$component)
})

test_that("staleness is reported per assay, not per ledger entry", {
    bv <- makeProfiled()

    w <- tryCatch(
        varianceTable(bv[, 1:6]),
        warning = function(w) conditionMessage(w)
    )

    ## Every assay shares one reason under subsetting, so the fact is
    ## stated once with the assays grouped under it.
    expect_equal(
        lengths(regmatches(w, gregexpr("computed on 20 samples", w)))[[1]],
        1L
    )
    expect_match(w, "'raw'")
    expect_match(w, "'raw_combat'")
})

test_that("varianceResults() warns only about the assays it reports", {
    bv <- makeProfiled()

    ## "raw" is untouched; "raw_combat" is replaced. Filtering to the
    ## former must not warn about the latter.
    mat <- SummarizedExperiment::assay(bv, "raw_combat")
    mat[1, 1] <- mat[1, 1] + 1
    SummarizedExperiment::assay(bv, "raw_combat") <- mat

    expect_no_warning(varianceResults(bv, assayName = "raw"))
    expect_warning(varianceResults(bv, assayName = "raw_combat"))
})

test_that("entries without a fingerprint do not warn", {
    bv <- makeProfiled()

    ## Objects built before fingerprints were recorded have none. An absent
    ## fingerprint is missing evidence, not evidence of a mismatch.
    vh <- S4Vectors::metadata(bv)$variance_history
    vh <- lapply(vh, function(e) {
        e$fingerprint <- NULL
        e
    })
    S4Vectors::metadata(bv)$variance_history <- vh

    expect_no_warning(varianceTable(bv[, 1:6]))
})

## ---- requirements checking ---------------------------------------------
##
## These rules used to be an S4 validity method, which could only run at
## construction and on an explicit validObject(). They now run at the
## point of use, so the mutations below are caught by the next public
## call rather than by nothing at all.

test_that("a clean object and a subset of one are accepted", {
    bv <- makeProfiled()

    expect_no_error(assayVariance(bv))

    ## Subsetting preserves every assay name and colData column, so the
    ## referential rules must not fire on it. Staleness is a separate
    ## question answered elsewhere in this file.
    expect_no_error(assayVariance(bv[, 1:6]))
})

test_that("a malformed ledger is rejected", {
    bv <- makeProfiled()

    b <- bv
    S4Vectors::metadata(b)$variance_history <- "not a ledger"
    expect_error(assayVariance(b), "variance_history is not a list")

    b <- bv
    S4Vectors::metadata(b)$correction_history <- list(list(method = 42))
    expect_error(assayVariance(b), "missing required field")

    b <- bv
    S4Vectors::metadata(b)$correction_history[[1]]$batch <- c("a", "b")
    expect_error(assayVariance(b), "single character string")
})

test_that("references to absent assays and columns are rejected", {
    ## The case the old validity method could not see. Neither `assay<-`
    ## nor `colData<-` invokes S4 validity, so before the check moved to
    ## the point of use these objects passed through every function in
    ## the package with a ledger describing things they no longer had.
    bv <- makeProfiled()

    b <- bv
    SummarizedExperiment::assay(b, "raw_combat") <- NULL
    expect_error(assayVariance(b), "refers to assay")

    b <- bv
    SummarizedExperiment::colData(b)$batch <- NULL
    expect_error(assayVariance(b), "refers to colData column")
})

test_that("a dangling reference is caught at the first public call", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = matrix(rnorm(40), 10, 4)),
        colData = S4Vectors::DataFrame(batch = rep(c("A", "B"), each = 2))
    )
    S4Vectors::metadata(se)$correction_history <- list(
        list(
            method = "combat", assay_in = "raw",
            assay_out = "absent", batch = "batch"
        )
    )

    expect_error(assayVariance(se), "refers to assay")
    expect_error(profileVariance(se, ~batch), "refers to assay")
})
