#' Access variance profiling history
#'
#' @param object BatchVariaData object
#'
#' @return list of variance analysis blocks
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' length(varianceHistory(bv))
#'
#' @export
varianceHistory <- function(object) {
    md <- S4Vectors::metadata(object)
    md$variance_history
}

## Per-feature variance, vectorised. Equivalent to apply(m, 1, stats::var)
## but avoids materialising one closure result per feature.
.rowVars <- function(m) {
    n <- ncol(m)
    if (n < 2L) {
        return(rep(NA_real_, nrow(m)))
    }
    mu <- rowMeans(m, na.rm = TRUE)
    rowSums((m - mu)^2, na.rm = TRUE) / (n - 1L)
}

## Total variance of an assay: the sum of per-feature variances, i.e. the
## trace of the feature covariance matrix. This is the denominator that
## variance *fractions* are expressed relative to.
.assayTotalVariance <- function(object, assayName) {
    m <- as.matrix(SummarizedExperiment::assay(object, assayName))
    sum(.rowVars(m), na.rm = TRUE)
}

#' Total variance of each assay
#'
#' Reports the total variance carried by each assay, independent of any
#' variance decomposition. Total variance is the sum of per-feature
#' variances, i.e. the trace of the feature covariance matrix.
#'
#' @param object A \code{BatchVariaData} object.
#' @param assays Character vector of assay names (default: all assays).
#'
#' @return A data.frame with columns \code{assay}, \code{total_variance},
#'   \code{n_features} and \code{n_samples}, one row per assay.
#'
#' @details
#' Variance fractions reported by \code{\link{profileVariance}} are
#' compositional: they sum to one within an assay. A correction that
#' removes variance associated with one term therefore *increases* the
#' fraction attributed to every remaining term, even when the absolute
#' variance of those terms is unchanged or reduced.
#'
#' Because of this, a rising variance fraction cannot on its own be read
#' as preserved or enhanced signal. \code{assayVariance()} supplies the
#' denominator, so that a fraction can be converted to an absolute
#' quantity:
#'
#' \deqn{absolute variance = variance fraction \times total variance}
#'
#' Comparing absolute variances across assays is immune to the
#' compositional effect and is the appropriate basis for judging whether
#' a correction preserved biological signal.
#'
#' Assays with fewer than two samples yield \code{NA} per-feature
#' variances, which are excluded from the total.
#'
#' @seealso \code{\link{varianceResults}}, \code{\link{profileVariance}}
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' assayVariance(bv, assays = c("raw", "raw_combat"))
#'
#' @export
assayVariance <- function(object, assays = NULL) {
    stopifnot(is(object, "BatchVariaData"))

    available <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available
    }

    missingAssays <- setdiff(assays, available)
    if (length(missingAssays) > 0) {
        stop("Assays not found: ", paste(missingAssays, collapse = ", "))
    }

    rows <- lapply(assays, function(a) {
        m <- as.matrix(SummarizedExperiment::assay(object, a))
        data.frame(
            assay = a,
            total_variance = sum(.rowVars(m), na.rm = TRUE),
            n_features = nrow(m),
            n_samples = ncol(m),
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#'
#' Retrieve variance decomposition results
#'
#' Extracts variance decomposition results stored in the
#' \code{BatchVariaData} object metadata.
#'
#' @param object A \code{BatchVariaData} object.
#' @param assayName Optional assay name to filter results.
#' @param method Optional method name to filter results.
#'
#' @return A data.frame containing variance decomposition results
#' in long format.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' varianceResults(bv, assayName = "raw")
#'
#' @export
varianceResults <- function(object, assayName = NULL, method = NULL) {
    stopifnot(is(object, "BatchVariaData"))

    vh <- S4Vectors::metadata(object)$variance_history

    if (is.null(vh) || length(vh) == 0) {
        return(data.frame())
    }

    out <- list()

    for (entry in vh) {
        # each entry corresponds to one recordVariance() call
        res <- entry$result

        if (is.null(res)) {
            next
        }

        df <- as.data.frame(res)

        df$assay <- entry$assay
        df$method <- entry$method
        df$formula <- .entryFormulaKey(entry)
        df$timestamp <- entry$timestamp

        out[[length(out) + 1]] <- df
    }

    ## Engines share the required columns but may carry their own
    ## alongside them -- the anova engine records the sums-of-squares type
    ## and weighting it used, for instance. bind_rows() fills those with NA
    ## for engines to which they do not apply, where rbind() would fail on
    ## the differing widths.
    out <- dplyr::bind_rows(out)

    # optional filtering
    if (!is.null(assayName)) {
        out <- out[out$assay %in% assayName, , drop = FALSE]
    }

    if (!is.null(method)) {
        out <- out[out$method %in% method, , drop = FALSE]
    }

    rownames(out) <- NULL

    out
}
