#' Available variance profiling methods
#'
#' Returns the names of supported variance decomposition engines.
#'
#' @return Character vector of method identifiers
#' @examples
#' availableVarianceMethods()
#' @export
availableVarianceMethods <- function() {
    c("pca", "anova", "variancePartition")
}

## Features with zero or non-finite variance carry no information for a
## variance decomposition, and each engine fails differently on them:
## prcomp(scale. = TRUE) cannot rescale a constant column, the anova
## engine divides by a zero total sum of squares, and variancePartition
## drops them silently inside colMeans(na.rm = TRUE). Excluding them once,
## up front, makes the three engines agree on which features were used.
##
## All-zero features are routine in unfiltered count data, so this is the
## common case rather than an edge case.
.dropConstantFeatures <- function(assayMatrix, assayName = NULL) {
    featureVar <- .rowVars(assayMatrix)
    keep <- is.finite(featureVar) & featureVar > 0

    where <- if (is.null(assayName)) "" else paste0(" in assay '", assayName, "'")

    if (!any(keep)) {
        stop(
            "All ", nrow(assayMatrix), " features", where,
            " have zero or non-finite variance; there is no variance to ",
            "decompose"
        )
    }

    nDropped <- sum(!keep)
    if (nDropped > 0) {
        warning(
            nDropped, " of ", nrow(assayMatrix), " features", where,
            " have zero or non-finite variance and were excluded from ",
            "variance profiling",
            call. = FALSE
        )
    }

    assayMatrix[keep, , drop = FALSE]
}

## Build a fixed-effects design matrix for engines that need one.
##
## This is deliberately per-engine rather than hoisted into
## profileVariance(): of the three engines only anova consumes a design
## matrix at all, and building one centrally meant a formula written for
## variancePartition was silently forced through stats::model.matrix().
## There, '|' is parsed as logical OR rather than as random-effects
## notation, which errors for character covariates but for factor and
## numeric covariates yields a nonsense design that anova fits without
## complaint.
.buildFixedModelMatrix <- function(formula, sampleData, method) {
    randomTerms <- reformulas::findbars(formula)

    if (length(randomTerms) > 0) {
        randomVars <- vapply(randomTerms, function(x) deparse1(x[[3L]]), character(1))
        stop(
            "the '", method, "' engine models fixed effects only and cannot ",
            "consume random-effects notation. Rewrite ",
            deparse1(formula), " with ",
            paste0(randomVars, collapse = " + "),
            " as fixed effects, or restrict this formula to methods that ",
            "support random effects (see availableVarianceMethods())"
        )
    }

    stats::model.matrix(formula, data = as.data.frame(sampleData))
}

.computePCAVariance <- function(assayMatrix, formula = NULL, sampleData = NULL, ...) {
    stopifnot(
        is.matrix(assayMatrix)
    )

    pcs <- stats::prcomp(
        t(assayMatrix),
        scale. = TRUE
    )

    percent_var <- pcs$sdev^2 / sum(pcs$sdev^2)

    .newVarianceSummary(
        source = "pca",
        term = paste0("PC", seq_along(percent_var)),
        varianceFraction = percent_var,
        n = ncol(assayMatrix)
    )
}

.computeAnovaVariance <- function(assayMatrix, formula, sampleData, ...) {
    stopifnot(is.matrix(assayMatrix))

    modelMatrix <- .buildFixedModelMatrix(formula, sampleData, "anova")

    stopifnot(
        is.matrix(modelMatrix),
        ncol(assayMatrix) == nrow(modelMatrix)
    )

    # Fit linear models per feature
    fits <- apply(
        assayMatrix,
        1,
        function(y) stats::lm.fit(modelMatrix, y)
    )

    # Residual sum of squares
    rss <- vapply(
        fits,
        function(f) sum(f$residuals^2),
        numeric(1)
    )

    # Total sum of squares
    tss <- apply(
        assayMatrix,
        1,
        function(y) sum((y - mean(y))^2)
    )

    # Fractions
    #
    # A feature with zero total sum of squares gives rss / tss == Inf, and
    # na.rm = TRUE does not remove Inf -- one such feature would drive the
    # mean to Inf and the model fraction to -Inf. profileVariance() screens
    # these out before calling any engine; this guard keeps the engine
    # correct when it is called directly.
    ratio <- rss / tss
    usable <- is.finite(ratio)

    if (!any(usable)) {
        stop(
            "No features with non-zero total variance; ",
            "cannot compute an anova variance decomposition"
        )
    }

    residual_fraction <- mean(ratio[usable])
    model_fraction <- 1 - residual_fraction

    # Canonical variance-summary output
    .newVarianceSummary(
        source = "anova",
        term = c("model", "residual"),
        varianceFraction = c(model_fraction, residual_fraction),
        n = nrow(assayMatrix)
    )
}

.computeVariancePartitionVariance <- function(
    assayMatrix,
    formula,
    sampleData,
    REML = TRUE,
    ...
) {
    # Check for variancePartition package
    if (!requireNamespace("variancePartition", quietly = TRUE)) {
        stop(
            "The 'variancePartition' package is required for this method. ",
            "Please install it from Bioconductor."
        )
    }

    # Validate inputs
    stopifnot(is.matrix(assayMatrix))
    stopifnot(inherits(formula, "formula"))

    # Convert sampleData to data.frame if needed
    sampleData <- as.data.frame(sampleData)

    # Ensure rownames are present
    if (is.null(rownames(sampleData))) {
        stop("sampleData must have rownames corresponding to sample names (columns of assayMatrix)")
    }

    # Align samples: reorder metadata to match assay columns
    if (!identical(colnames(assayMatrix), rownames(sampleData))) {
        if (all(colnames(assayMatrix) %in% rownames(sampleData))) {
            sampleData <- sampleData[colnames(assayMatrix), , drop = FALSE]
        } else {
            stop(
                "Sample names mismatch: columns of assayMatrix and rows of colData do not match"
            )
        }
    }

    # Detect categorical variables in sampleData.
    #
    # Character columns count as categorical: as.data.frame() on a DataFrame
    # does not convert them to factors under R >= 4.0, so testing is.factor()
    # alone left this guard unreachable for the common case of colData built
    # from character vectors.
    categorical_cols <- names(Filter(
        function(x) is.factor(x) || is.character(x),
        sampleData
    ))

    # Extract random effect variable names
    random_vars <- vapply(
        reformulas::findbars(formula),
        function(x) deparse1(x[[3L]]),
        character(1)
    )

    # Extract fixed-effect variables
    fixed_vars <- all.vars(stats::delete.response(stats::terms(formula)))

    # Identify fixed categorical variables that are NOT random
    fixed_categorical <- setdiff(
        intersect(fixed_vars, categorical_cols),
        random_vars
    )

    # Defensive check: cannot mix categorical fixed + random
    if (length(random_vars) > 0 && length(fixed_categorical) > 0) {
        stop(
            "variancePartition does not allow mixing categorical fixed and random effects.\n",
            "Use a formula like: ~ (1 | batch) + (1 | condition)\n",
            "or: ~ batch + condition (all fixed)"
        )
    }

    # Run variancePartition
    vp <- variancePartition::fitExtractVarPartModel(
        exprObj = assayMatrix,
        formula = formula,
        data    = sampleData,
        REML    = REML,
        ...
    )

    # Compute mean variance fractions across features
    variance_fraction <- colMeans(vp, na.rm = TRUE)

    # Return tidy variance summary contract
    .newVarianceSummary(
        source = "variancePartition",
        term = names(variance_fraction),
        varianceFraction = as.numeric(variance_fraction),
        n = nrow(assayMatrix)
    )
}
