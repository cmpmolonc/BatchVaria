#'
#' Summary method for BatchVariaData
#'
#' Provides a structured overview of assays, correction history,
#' and variance profiling analyses.
#'
#' @param object A BatchVariaData object
#' @param ... unused
#'
#' @return Invisibly returns a list containing summary tables
#'
#' @importFrom S4Vectors summary
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
#' summary(bv)
#' @export
setMethod(
    "summary",
    "BatchVariaData",
    function(object, ...) {
        se <- object

        cat("BatchVariaData summary\n")
        cat("----------------------\n")

        ## ---- basic structure --------------------------------------------
        cat("\nDimensions:\n")
        cat("  features:", nrow(se), "\n")
        cat("  samples :", ncol(se), "\n")

        assays <- SummarizedExperiment::assayNames(se)

        cat("\nAssays (", length(assays), "):\n", sep = "")
        print(assays)

        ## ---- correction history -----------------------------------------
        ch <- metadata(se)$correction_history

        correction_table <- NULL

        if (!is.null(ch) && length(ch) > 0) {
            correction_table <- do.call(
                rbind,
                lapply(seq_along(ch), function(i) {
                    entry <- ch[[i]]

                    data.frame(
                        step = i,
                        method = entry$method,
                        input = entry$assay_in,
                        output = entry$assay_out,
                        batch = entry$batch,
                        timestamp = as.character(entry$timestamp),
                        stringsAsFactors = FALSE
                    )
                })
            )

            cat("\nCorrection history:\n")
            print(correction_table, row.names = FALSE)
        } else {
            cat("\nCorrection history: none\n")
        }

        ## ---- variance history -------------------------------------------
        vh <- metadata(se)$variance_history

        variance_table <- NULL

        if (!is.null(vh) && length(vh) > 0) {
            variance_table <- do.call(
                rbind,
                lapply(seq_along(vh), function(i) {
                    entry <- vh[[i]]

                    methods <- entry$method

                    data.frame(
                        step = i,
                        assay = entry$assay,
                        formula = deparse(entry$formula),
                        methods = methods,
                        timestamp = as.character(entry$timestamp),
                        stringsAsFactors = FALSE
                    )
                })
            )

            cat("\nVariance profiling history:\n")
            print(variance_table, row.names = FALSE)
        } else {
            cat("\nVariance profiling history: none\n")
        }

        invisible(
            list(
                assays = assays,
                corrections = correction_table,
                variance_profiles = variance_table
            )
        )
    }
)
