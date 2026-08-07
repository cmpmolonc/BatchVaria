# Internal contract for variance decomposition results
# This file defines the canonical structure expected by
# plotting, storage, and downstream analysis layers.
# Returns data frame with each row corresponding to model term
# validated via .validateVarianceSummary
#
# Columns:
#   source            engine that produced the row
#   term              model term the variance is attributed to
#   variance_fraction share of total variance, in [0, 1] except for the
#                     reserved 'shared' term (see below)
#   n_features        number of features the decomposition was computed
#                     over, after zero-variance features were excluded
#
# n_features is required rather than optional: it appears in
# varianceResults() output beside engine-specific columns, and a bare 'n'
# there was ambiguous -- it meant samples for one engine and features for
# the others. Engines may add their own columns alongside these.
#
# Reserved term names. Two terms mean the same thing in every engine and
# are therefore spelled the same way in every engine:
#
#   residual  variance no model term accounts for
#   shared    variance an unbalanced design cannot attribute to any
#             single term; may be negative, unlike any other term
#
# An engine whose underlying implementation names these differently --
# variancePartition returns "Residuals" -- normalises before returning.
# Harmonising the columns but not the vocabulary would leave results that
# look comparable and cannot be joined, which is worse than not
# harmonising at all. Engines are free to name every other term as their
# model does.

.validateVarianceSummary <- function(summaryDf) {
    required <- c("source", "term", "variance_fraction")

    if (!is.data.frame(summaryDf)) {
        stop("variance summary must be a data.frame")
    }

    missing <- setdiff(required, colnames(summaryDf))
    if (length(missing) > 0) {
        stop(
            "variance summary missing required columns: ",
            paste(missing, collapse = ", ")
        )
    }

    if (!is.numeric(summaryDf$variance_fraction)) {
        stop("variance_fraction must be numeric")
    }

    ## 'shared' is reserved for variance an unbalanced design cannot
    ## attribute to any single term. It is normally non-negative, but goes
    ## negative under suppression, where terms explain more jointly than
    ## separately. Reporting that is more honest than clamping it, so the
    ## reserved term is exempt from this check.
    attributed <- summaryDf[summaryDf$term != "shared", , drop = FALSE]

    if (any(attributed$variance_fraction < 0, na.rm = TRUE)) {
        stop("variance_fraction must be non-negative")
    }

    invisible(TRUE)
}

#' Build a conforming variance summary
#'
#' Constructs the data.frame that \code{\link{profileVariance}} expects
#' from a variance engine, and validates it. Use this rather than
#' assembling the frame by hand when writing an engine.
#'
#' @param source Character. Name of the engine producing the result.
#' @param term Character vector of model terms.
#' @param varianceFraction Numeric vector of variance shares, one per term.
#' @param nFeatures Number of features the decomposition was computed over.
#' @param level Optional character vector of factor levels per term.
#'
#' @return A validated data.frame with columns \code{source},
#'   \code{term}, \code{variance_fraction} and \code{n_features}.
#'
#' @seealso \code{\link{registerVarianceEngine}} for the full contract.
#'
#' @examples
#' newVarianceSummary(
#'     source = "example",
#'     term = c("batch", "residual"),
#'     varianceFraction = c(0.3, 0.7),
#'     nFeatures = 100
#' )
#'
#' @export
newVarianceSummary <- function(
    source,
    term,
    varianceFraction,
    nFeatures,
    level = NULL
) {
    if (missing(nFeatures) || !is.numeric(nFeatures) || anyNA(nFeatures)) {
        stop("'nFeatures' must be supplied as a number of features")
    }

    ## 'variance_fraction' is the documented column name in the result
    ## contract, so it stays snake_case even though the argument does not
    df <- data.frame(
        source = source,
        term = term,
        variance_fraction = varianceFraction,
        n_features = nFeatures,
        stringsAsFactors = FALSE
    )

    if (!is.null(level)) df$level <- level

    .validateVarianceSummary(df)
    df
}
