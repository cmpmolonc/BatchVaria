#' Available variance profiling methods
#'
#' Returns the names of supported variance decomposition engines.
#'
#' @return Character vector of method identifiers
#'
#' @details
#' The list is the contents of the engine registry, so it grows when
#' \code{\link{registerVarianceEngine}} is called. Every engine listed here answers the same question -- how much of the
#' variance is attributable to each term of a model -- so their results
#' share a schema and can be compared directly.
#'
#' PCA is deliberately not among them. It is unsupervised, so its output is
#' not attributable to any covariate, and its terms are latent axes rather
#' than model terms: putting both in one \code{variance_fraction} column
#' under one \code{term} key would assert a commensurability that does not
#' hold. PCA remains available as an embedding through
#' \code{\link{comparePCA}}, \code{\link{basisRetention}} and the
#' plotting layer, where a shared reference basis makes cross-assay
#' comparison well defined.
#'
#' @seealso \code{\link{registerVarianceEngine}} to add one;
#'   \code{\link{comparePCA}} for the PCA-based view.
#'
#' @examples
#' availableVarianceMethods()
#' @export
availableVarianceMethods <- function() {
    sort(ls(.varianceEngines))
}

## Features with zero or non-finite variance carry no information for a
## variance decomposition, and each engine fails differently on them: the
## anova engine divides by a zero total sum of squares, and
## variancePartition drops them silently inside colMeans(na.rm = TRUE).
## Excluding them once, up front, makes both engines agree on which
## features were used.
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
## profileVariance(): only the anova engine consumes a design matrix at
## all, and building one centrally meant a formula written for
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

## Residual sum of squares for every column of Yt against a design.
##
## One QR decomposition serves all features, rather than one lm.fit object
## per feature: the old approach held residuals, effects and a QR for every
## gene before reducing each to a single number.
.rssPerFeature <- function(Yt, X) {
    if (ncol(X) == 0L) {
        return(colSums(Yt^2))
    }
    colSums(qr.resid(qr(X), Yt)^2)
}

## Which model terms contain term j, in the marginality sense used by
## Type II sums of squares: term j is contained in term k when j's
## variables are a strict subset of k's.
.termsContaining <- function(factors, j) {
    inJ <- factors[, j] > 0
    vapply(
        seq_len(ncol(factors)),
        function(k) {
            inK <- factors[, k] > 0
            all(inJ <= inK) && !all(inJ == inK)
        },
        logical(1)
    )
}

.computeAnovaVariance <- function(
    assayMatrix,
    formula,
    sampleData,
    ssType = "II",
    weighting = c("feature", "pooled"),
    ...
) {
    stopifnot(is.matrix(assayMatrix))
    weighting <- match.arg(weighting)

    if (!identical(ssType, "II")) {
        stop(
            "Only Type II sums of squares are supported. Type III requires ",
            "sum-to-zero contrasts to be meaningful and silently produces ",
            "misleading results under R's default treatment contrasts"
        )
    }

    modelMatrix <- .buildFixedModelMatrix(formula, sampleData, "anova")

    stopifnot(
        is.matrix(modelMatrix),
        ncol(assayMatrix) == nrow(modelMatrix)
    )

    tt <- stats::terms(formula)
    labels <- attr(tt, "term.labels")
    factors <- attr(tt, "factors")
    assign <- attr(modelMatrix, "assign")

    if (length(labels) == 0L) {
        stop("The formula has no terms to decompose variance across")
    }

    ## Features are columns from here on: one QR serves all of them.
    Yt <- t(assayMatrix)

    ## Total sum of squares is the residual of the null model, which makes
    ## intercept and no-intercept designs fall out of the same code. With
    ## an intercept this is the usual sum of squares about the mean; with
    ## ~0 + x it is the uncorrected sum of squares, so that the fractions
    ## remain a partition of what the model is actually fitting.
    tss <- .rssPerFeature(Yt, modelMatrix[, assign == 0L, drop = FALSE])
    rssFull <- .rssPerFeature(Yt, modelMatrix)

    usable <- is.finite(tss) & tss > 0

    if (!any(usable)) {
        stop(
            "No features with non-zero total variance; ",
            "cannot compute an anova variance decomposition"
        )
    }

    ## Type II: each term is tested against a model holding every term
    ## that does not contain it. Under an unbalanced design these sums of
    ## squares do not add up to the model sum of squares -- variance that
    ## the design cannot attribute to any single term is left over, and is
    ## reported below as 'shared' rather than silently dropped.
    ssTerms <- vapply(seq_along(labels), function(j) {
        containing <- .termsContaining(factors, j)
        baseTerms <- setdiff(which(!containing), j)

        baseCols <- assign == 0L | assign %in% baseTerms
        fullCols <- baseCols | assign == j

        .rssPerFeature(Yt, modelMatrix[, baseCols, drop = FALSE]) -
            .rssPerFeature(Yt, modelMatrix[, fullCols, drop = FALSE])
    }, numeric(ncol(Yt)))

    ssTerms <- matrix(ssTerms, nrow = ncol(Yt), dimnames = list(NULL, labels))

    ## Whatever the terms and the residual do not account for
    ssShared <- tss - rowSums(ssTerms) - rssFull

    ssAll <- cbind(ssTerms, shared = ssShared, residual = rssFull)

    ## 'feature' averages per-feature fractions, matching how
    ## variancePartition summarises across features, so the two engines
    ## describe the typical feature and remain comparable. 'pooled' sums
    ## the sums of squares first, weighting each feature by how much
    ## variance it actually carries.
    fractions <- if (weighting == "feature") {
        colMeans(ssAll[usable, , drop = FALSE] / tss[usable])
    } else {
        colSums(ssAll[usable, , drop = FALSE]) / sum(tss[usable])
    }

    ## Rounding leaves fractions a hair from zero; a materially negative
    ## shared component is real, and means the terms jointly explain more
    ## than they do apart (a suppression effect).
    tol <- sqrt(.Machine$double.eps)
    fractions[abs(fractions) < tol] <- 0

    if (fractions[["shared"]] < -tol) {
        warning(
            "The shared variance component is negative (",
            format(fractions[["shared"]], digits = 3),
            "), meaning the model terms explain more jointly than ",
            "separately. This indicates a suppression effect in the design",
            call. = FALSE
        )
    }

    out <- newVarianceSummary(
        source = "anova",
        term = names(fractions),
        varianceFraction = as.numeric(fractions),
        nFeatures = sum(usable)
    )

    ## The analysis choices travel with the estimate rather than being
    ## resolved out of sight.
    out$ss_type <- ssType
    out$weighting <- weighting

    out
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
    newVarianceSummary(
        source = "variancePartition",
        term = names(variance_fraction),
        varianceFraction = as.numeric(variance_fraction),
        nFeatures = nrow(assayMatrix)
    )
}
