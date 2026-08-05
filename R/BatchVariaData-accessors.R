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
        df$timestamp <- entry$timestamp

        out[[length(out) + 1]] <- df
    }

    out <- do.call(rbind, out)

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
