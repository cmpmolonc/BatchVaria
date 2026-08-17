## Design matrices are built per engine. Building one centrally meant a
## random-effects formula written for variancePartition was forced through
## stats::model.matrix(), where '|' is logical OR rather than random-effects
## notation.

withBatchType <- function(bv, type) {
    cd <- SummarizedExperiment::colData(bv)
    cd$batch <- switch(type,
        character = as.character(cd$batch),
        factor = factor(cd$batch),
        numeric = as.numeric(factor(cd$batch))
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = as.matrix(SummarizedExperiment::assay(bv, "raw"))),
        colData = cd
    )
    se
}


test_that("random-effects formulas reach the engines that support them", {
    skip_if_not_installed("variancePartition")

    ## ~(1 | batch) is the formula variancePartition's own error message
    ## recommends. It previously failed in profileVariance() before any
    ## engine ran, because the design matrix was built centrally.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 80)

    res <- suppressWarnings(suppressMessages(
        profileVariance(
            bv, ~(1 | batch),
            assays = "raw",
            methods = c("anova", "variancePartition")
        )
    ))

    methods <- vapply(varianceHistory(res), function(x) x$method, character(1))
    expect_true("variancePartition" %in% methods)
    expect_false("anova" %in% methods)

    vp <- Filter(function(x) x$method == "variancePartition", varianceHistory(res))
    expect_setequal(vp[[1]]$result$term, c("batch", "residual"))
})


test_that("the anova engine refuses random effects rather than mis-fitting", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    expect_warning(
        try(
            suppressMessages(
                profileVariance(bv, ~(1 | batch), assays = "raw", methods = "anova")
            ),
            silent = TRUE
        ),
        "models fixed effects only"
    )
})


test_that("random-effects refusal does not depend on covariate type", {
    ## Previously stats::model.matrix(~(1 | batch)) errored for character
    ## covariates, produced a zero-row design for factors, and produced a
    ## constant column for numerics -- and anova silently fitted the last
    ## two. All three must now behave identically.
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    for (type in c("character", "factor", "numeric")) {
        reasons <- character()

        tryCatch(
            withCallingHandlers(
                suppressMessages(profileVariance(
                    withBatchType(bv, type), ~(1 | batch), methods = "anova"
                )),
                warning = function(w) {
                    reasons <<- c(reasons, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }
            ),
            error = function(e) NULL
        )

        expect_match(
            paste(reasons, collapse = " "),
            "models fixed effects only",
            info = paste("covariate type:", type)
        )
    }
})


test_that("fixed-effects formulas are unaffected", {
    skip_if_not_installed("variancePartition")

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    res <- suppressWarnings(suppressMessages(
        profileVariance(
            bv, ~ batch + group,
            assays = "raw",
            methods = c("anova", "variancePartition")
        )
    ))

    methods <- vapply(varianceHistory(res), function(x) x$method, character(1))
    expect_setequal(methods, c("anova", "variancePartition"))

    vp <- Filter(function(x) x$method == "variancePartition", varianceHistory(res))
    expect_setequal(vp[[1]]$result$term, c("batch", "group", "residual"))
})


test_that("the categorical mixing guard fires on character colData", {
    ## exampleBatchVaria() builds colData from character vectors, which
    ## as.data.frame() does not convert to factors under R >= 4.0. Testing
    ## is.factor() alone left this guard unreachable on the package's own
    ## demonstration data.
    skip_if_not_installed("variancePartition")

    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 60)

    expect_true(is.character(
        SummarizedExperiment::colData(bv)$group
    ))

    expect_warning(
        try(
            suppressMessages(profileVariance(
                bv, ~ (1 | batch) + group,
                assays = "raw", methods = "variancePartition"
            )),
            silent = TRUE
        ),
        "does not allow mixing categorical fixed and random effects"
    )
})


test_that("profileVariance validates formula variables against colData", {
    set.seed(1)
    bv <- exampleBatchVaria(nGenes = 40)

    expect_error(
        profileVariance(bv, ~nosuchvariable, methods = "anova"),
        "Formula variables not found in colData: nosuchvariable"
    )

    expect_error(
        profileVariance(bv, "batch", methods = "anova"),
        "'formula' must be a formula"
    )
})


test_that(".buildFixedModelMatrix names the fixed-effect rewrite", {
    d <- data.frame(batch = c("A", "A", "B", "B"))

    expect_error(
        BatchVaria:::.buildFixedModelMatrix(~(1 | batch), d, "anova"),
        "Rewrite .* with batch as fixed effects"
    )

    mm <- BatchVaria:::.buildFixedModelMatrix(~batch, d, "anova")
    expect_equal(colnames(mm), c("(Intercept)", "batchB"))
    expect_equal(nrow(mm), 4)
})
