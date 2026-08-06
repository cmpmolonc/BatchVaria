#'
#' Available correction methods
#'
#' Returns the names of the registered batch correction methods.
#'
#' @return Character vector of method identifiers
#' @seealso \code{\link{registerCorrectionMethod}}
#' @examples
#' availableCorrectionMethods()
#' @export

availableCorrectionMethods <- function() {
    sort(ls(.correctionMethods))
}


.computeCombat <- function(
    assayMatrix,
    batch,
    sampleData,
    preserve = NULL,
    parPrior = TRUE,
    priorPlots = FALSE,
    ...
) {
    if (!requireNamespace("sva", quietly = TRUE)) {
        stop(
            "The 'sva' package is required for ComBat. ",
            "Please install it from Bioconductor."
        )
    }

    sva::ComBat(
        dat = assayMatrix,
        batch = sampleData[[batch]],
        mod = .preserveDesign(preserve, sampleData),
        par.prior = parPrior,
        prior.plots = priorPlots
    )
}


.computeLimma <- function(
    assayMatrix,
    batch,
    sampleData,
    preserve = NULL,
    ...
) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop(
            "The 'limma' package is required for the 'limma' correction ",
            "method. Please install it from Bioconductor."
        )
    }

    ## 'preserve' maps to 'design' -- the conditions to keep -- and never to
    ## removeBatchEffect()'s 'covariates', which means the opposite.
    design <- .preserveDesign(preserve, sampleData)

    if (is.null(design)) {
        ## removeBatchEffect() defaults to exactly this but messages about
        ## it; building it here keeps an unprotected correction quiet
        design <- matrix(1, ncol(assayMatrix), 1)
    }

    limma::removeBatchEffect(
        x = assayMatrix,
        batch = sampleData[[batch]],
        design = design,
        ...
    )
}
