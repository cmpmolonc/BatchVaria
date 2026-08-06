#'
#' Apply batch correction to expression data
#'
#' Runs a specified batch correction method on an assay and stores the
#' corrected matrix as a new assay within the \code{BatchVariaData} object.
#'
#' @param object A \code{BatchVariaData} object.
#' @param method Character string specifying the correction method
#' @param batch Column name in \code{colData(object)} defining batch structure.
#' @param covariates Optional character vector of additional covariates.
#' @param assayName Name of the input assay (default \code{"raw"}).
#' @param newAssayName Optional corrected assay name.
#' @param ... Additional arguments passed to the correction engine.
#'
#' @return Updated \code{BatchVariaData} object with corrected assay added.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' SummarizedExperiment::assayNames(bv)
#'
#' @export
runCorrection <- function(
    object,
    method = "combat",
    batch,
    covariates = NULL,
    assayName = "raw",
    newAssayName = NULL,
    ...
) {
    stopifnot(is(object, "BatchVariaData"))

    if (!assayName %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay not found: ", assayName)
    }

    sampleData <- as.data.frame(SummarizedExperiment::colData(object))

    if (!batch %in% colnames(sampleData)) {
        stop("Batch variable '", batch, "' not found in colData")
    }

    if (!is.null(covariates)) {
        missing <- setdiff(covariates, colnames(sampleData))

        if (length(missing) > 0) {
            stop(
                "Covariates not found in colData: ",
                paste(missing, collapse = ", ")
            )
        }
    }

    assayMatrix <- as.matrix(
        SummarizedExperiment::assay(object, assayName)
    )

    ## ---- Model matrix ---------------------------------------------------

    modelMatrix <- NULL

    if (!is.null(covariates)) {
        modelMatrix <- stats::model.matrix(
            stats::as.formula(
                paste("~", paste(covariates, collapse = " + "))
            ),
            data = sampleData
        )
    }

    ## ---- Select correction engine --------------------------------------

    engine_fun <- .getCorrectionEngine(method)

    engine_res <- engine_fun(
        assayMatrix = assayMatrix,
        batch = sampleData[[batch]],
        modelMatrix = modelMatrix,
        ...
    )

    corrected <- engine_res$corrected_assay

    ## ---- Name corrected assay ------------------------------------------

    if (is.null(newAssayName)) {
        newAssayName <- paste0(assayName, "_", method)
    }

    ## ---- Store corrected assay -----------------------------------------

    SummarizedExperiment::assay(object, newAssayName) <- corrected

    ## ---- Update correction ledger --------------------------------------

    ## Total variance before and after correction. Recorded here because
    ## variance fractions are compositional: without the denominator a
    ## rising fraction cannot be distinguished from preserved signal.
    entry <- list(
        method             = method,
        assay_in           = assayName,
        assay_out          = newAssayName,
        batch              = batch,
        covariates         = covariates,
        params             = list(...),
        total_variance_in  = .assayTotalVariance(object, assayName),
        total_variance_out = .assayTotalVariance(object, newAssayName),
        timestamp          = Sys.time()
    )

    if (is.null(metadata(object)$correction_history)) {
        metadata(object)$correction_history <- list()
    }

    metadata(object)$correction_history <-
        c(metadata(object)$correction_history, list(entry))

    object
}
