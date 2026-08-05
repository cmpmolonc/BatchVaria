#'
#' Show method for BatchVariaData
#'
#' Provides a concise summary of the container including assays,
#' sample annotations, and analysis history.
#'
#' @import methods
#' @param object A BatchVariaData object
#' @return prints a summary of \code{object}. Returns \code{object} invisibly.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv
#'
#' @export
setMethod("show", "BatchVariaData", function(object) {
    cat("class:", class(object), "\n")

    cat(
        "dimensions:",
        nrow(object), "features x",
        ncol(object), "samples\n"
    )

    ## ---- assays ------------------------------------------------------
    assays <- SummarizedExperiment::assayNames(object)

    cat("\nassays(", length(assays), "):\n", sep = "")
    print(utils::head(assays, 6))

    if (length(assays) > 6) {
        cat("...\n")
    }

    ## ---- colData -----------------------------------------------------
    cd <- colnames(SummarizedExperiment::colData(object))

    cat("\ncolData variables(", length(cd), "):\n", sep = "")
    print(utils::head(cd, 6))

    if (length(cd) > 6) {
        cat("...\n")
    }

    ## ---- history -----------------------------------------------------
    ch <- metadata(object)$correction_history
    vh <- metadata(object)$variance_history

    n_corr <- if (is.null(ch)) 0 else length(ch)
    n_var <- if (is.null(vh)) 0 else length(vh)

    last_method <- NULL
    if (n_corr > 0) {
        last_method <- ch[[n_corr]]$method
    }

    cat("\nanalysis history:\n")

    if (n_corr > 0) {
        cat(
            "  corrections:", n_corr,
            "(last:", last_method, ")\n"
        )
    } else {
        cat("  corrections:", n_corr, "\n")
    }

    cat("  variance profiles:", n_var, "\n")
    invisible(object)
})
