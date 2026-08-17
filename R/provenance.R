#' Report the analysis history recorded on an object
#'
#' Prints what has been done to an object: its dimensions and assays, the
#' corrections that produced each assay, and the variance decompositions
#' recorded against them. Returns the two tables invisibly, as
#' data.frames.
#'
#' This is the readable view of the two ledgers BatchVaria writes into
#' \code{metadata()}. \code{\link{varianceHistory}} returns the variance
#' ledger unprocessed; the correction ledger is
#' \code{metadata(object)$correction_history}.
#'
#' @param object A \code{SummarizedExperiment} meeting BatchVaria's
#'   requirements. See \link{BatchVaria-requirements}.
#'
#' @return Invisibly, a list with elements \code{assays} (character),
#'   \code{corrections} and \code{variance_profiles} (data.frames, or
#'   \code{NULL} where the ledger is empty).
#'
#' @details
#' The correction table reports \code{preserve} and \code{no_op} for each
#' step. \code{preserve} is what distinguishes two runs of one method on
#' one assay, and \code{no_op} marks a correction that returned an assay
#' indistinguishable from its input - a condition that is easy to miss
#' on the console and would otherwise be invisible once the session that
#' produced it has ended.
#'
#' @seealso \code{\link{varianceHistory}}, \code{\link{varianceResults}}
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' bv <- profileVariance(bv, ~batch, assays = c("raw", "raw_combat"))
#' provenance(bv)
#'
#' @export
provenance <- function(object) {
    .check_se(object)

    cat("BatchVaria provenance\n")
    cat("---------------------\n")

    ## ---- basic structure ------------------------------------------------
    cat("\nDimensions:\n")
    cat("  features:", nrow(object), "\n")
    cat("  samples :", ncol(object), "\n")

    assays <- SummarizedExperiment::assayNames(object)

    cat("\nAssays (", length(assays), "):\n", sep = "")
    print(assays)

    ## ---- correction history ---------------------------------------------
    ch <- S4Vectors::metadata(object)$correction_history

    correctionTable <- NULL

    if (!is.null(ch) && length(ch) > 0) {
        correctionTable <- do.call(
            rbind,
            lapply(seq_along(ch), function(i) {
                entry <- ch[[i]]

                data.frame(
                    step = i,
                    method = entry$method,
                    input = entry$assay_in,
                    output = entry$assay_out,
                    batch = entry$batch,
                    ## the field that distinguishes two runs of the
                    ## same method on the same assay
                    preserve = if (is.null(entry$preserve)) {
                        "-"
                    } else {
                        paste(entry$preserve, collapse = ", ")
                    },
                    no_op = isTRUE(entry$no_op),
                    timestamp = as.character(entry$timestamp),
                    stringsAsFactors = FALSE
                )
            })
        )

        cat("\nCorrection history:\n")
        print(correctionTable, row.names = FALSE)
    } else {
        cat("\nCorrection history: none\n")
    }

    ## ---- variance history -----------------------------------------------
    vh <- S4Vectors::metadata(object)$variance_history

    varianceProfileTable <- NULL

    if (!is.null(vh) && length(vh) > 0) {
        varianceProfileTable <- do.call(
            rbind,
            lapply(seq_along(vh), function(i) {
                entry <- vh[[i]]

                data.frame(
                    step = i,
                    assay = entry$assay,
                    formula = deparse(entry$formula),
                    methods = entry$method,
                    timestamp = as.character(entry$timestamp),
                    stringsAsFactors = FALSE
                )
            })
        )

        cat("\nVariance profiling history:\n")
        print(varianceProfileTable, row.names = FALSE)
    } else {
        cat("\nVariance profiling history: none\n")
    }

    invisible(
        list(
            assays = assays,
            corrections = correctionTable,
            variance_profiles = varianceProfileTable
        )
    )
}
