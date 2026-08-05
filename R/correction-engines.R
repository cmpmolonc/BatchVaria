#'
#' Available correction methods
#'
#' Returns the names of supported correction engines.
#'
#' @return Character vector of method identifiers
#' @examples
#' availableCorrectionMethods()
#' @export

availableCorrectionMethods <- function() {
    c("combat")
}


.getCorrectionEngine <- function(method) {
    switch(method,
        combat = .computeCombat,
        stop("Unknown correction method: ", method)
    )
}

.computeCombat <- function(
    assayMatrix,
    batch,
    modelMatrix = NULL,
    parPrior = TRUE,
    priorPlots = FALSE
) {
    if (!requireNamespace("sva", quietly = TRUE)) {
        stop(
            "The 'sva' package is required for ComBat. ",
            "Please install it from Bioconductor."
        )
    }

    stopifnot(
        is.matrix(assayMatrix),
        length(batch) == ncol(assayMatrix)
    )

    corrected <- sva::ComBat(
        dat = assayMatrix,
        batch = batch,
        mod = modelMatrix,
        par.prior = parPrior,
        prior.plots = priorPlots
    )

    list(
        corrected_assay = corrected,
        engine_output = list()
    )
}
