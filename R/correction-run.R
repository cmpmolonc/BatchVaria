#'
#' Apply batch correction to expression data
#'
#' Runs a specified batch correction method on an assay and stores the
#' corrected matrix as a new assay on the object.
#'
#' @param object A \code{SummarizedExperiment} meeting BatchVaria's
#'   requirements. See \link{BatchVaria-requirements}.
#' @param method Character string specifying the correction method
#' @param batch Column name in \code{colData(object)} defining batch structure.
#' @param preserve Optional character vector of \code{colData} column names
#'   whose variation should be \emph{kept}. The batch effect is estimated
#'   adjusted for these, so biological signal they carry is protected rather
#'   than removed.
#' @param assayName Name of the input assay (default \code{"raw"}).
#' @param newAssayName Optional corrected assay name.
#' @param ... Additional arguments passed to the correction engine.
#' @param verbose Show progress output from the correction engine
#'   (default \code{FALSE}). This governs only what the engine reports
#'   about its progress. Diagnostics BatchVaria derives itself are
#'   warnings and are always raised, whatever \code{verbose} is set to.
#'
#' @return Updated \code{SummarizedExperiment} with corrected assay added.
#'
#' @details
#' Two methods are registered on load. \code{"combat"} calls
#' \code{sva::ComBat()}, which applies empirical Bayes shrinkage to each
#' batch's location \emph{and} scale. \code{"limma"} calls
#' \code{limma::removeBatchEffect()}, which removes a per-feature location
#' shift only. The two therefore differ in magnitude even where they agree
#' in direction, and neither should be expected to reproduce the other.
#'
#' Both expect log-scale data. Neither is appropriate on raw counts, and
#' \code{removeBatchEffect()} in particular is intended for exploratory work
#' - visualisation, clustering, distance calculations - rather than as
#' input to differential expression testing, where the batch term belongs in
#' the model instead. Correcting first and testing afterwards understates
#' the residual variance and overstates significance.
#'
#' \code{removeBatchEffect()} fits \code{batch} and the \code{preserve}
#' design simultaneously, so a design in which the two are close to aliased
#' is rank deficient and the batch effect is not separately estimable. That
#' is a property of the design rather than of the method.
#'
#' Further methods can be added with
#' \code{\link{registerCorrectionMethod}}.
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
    preserve = NULL,
    assayName = "raw",
    newAssayName = NULL,
    ...,
    verbose = FALSE
) {
    .check_se(object)

    if (!assayName %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay not found: ", assayName)
    }

    ## One batch variable only. ComBat cannot take more, so a multi-batch
    ## contract would describe one method rather than the layer; methods
    ## that accept a second batch expose it through '...'.
    if (!is.character(batch) || length(batch) != 1L || is.na(batch)) {
        stop("'batch' must be a single column name in colData")
    }

    sampleData <- as.data.frame(SummarizedExperiment::colData(object))

    if (!batch %in% colnames(sampleData)) {
        stop("Batch variable '", batch, "' not found in colData")
    }

    if (!is.null(preserve)) {
        missing <- setdiff(preserve, colnames(sampleData))

        if (length(missing) > 0) {
            stop(
                "Variables to preserve not found in colData: ",
                paste(missing, collapse = ", ")
            )
        }
    }

    assayMatrix <- as.matrix(
        SummarizedExperiment::assay(object, assayName)
    )

    ## ---- Apply the correction method ------------------------------------

    ## Methods receive column names and build their own design matrix, so a
    ## design assembled for one is never forced through another.
    correctionFun <- .getCorrectionMethod(method)

    corrected <- .withEngineOutput(verbose, correctionFun(
        assayMatrix = assayMatrix,
        batch = batch,
        sampleData = sampleData,
        preserve = preserve,
        ...
    ))

    corrected <- .validateCorrectionResult(corrected, assayMatrix, method)

    ## ---- No-op detection -------------------------------------------------

    ## A method can succeed and change nothing: limma::removeBatchEffect()
    ## does exactly this when batch is aliased with the preserved design,
    ## reporting it through message() rather than a condition a caller can
    ## trap. Left undetected the result flows into the ledger, the delta
    ## tables and the PCA basis as a correction, and every comparison
    ## against the input assay is then a comparison of an assay with
    ## itself. Flagged rather than refused, since a correction that finds
    ## nothing to remove is legitimate on some data.
    ##
    ## all.equal() rather than identical(), deliberately, and the opposite
    ## choice from .validateCorrectionResult()'s dimnames check. That check
    ## asks whether the method altered something it must not, where any
    ## drift is a defect and a tolerance would hide it. This asks whether
    ## the method did anything at all, where a sub-tolerance change is
    ## functionally nothing. Different questions, so do not align them.
    noOp <- isTRUE(all.equal(corrected, assayMatrix))

    if (noOp) {
        warning(
            "Correction method '", method, "' returned an assay ",
            "indistinguishable from '", assayName, "'. Comparisons between ",
            "them are comparisons of an assay with itself. This happens ",
            "when the design leaves the batch effect inestimable, typically ",
            "because batch is aliased with the variables being preserved",
            call. = FALSE
        )
    }

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
        preserve           = preserve,
        no_op              = noOp,
        params             = list(...),
        total_variance_in  = .assayTotalVariance(object, assayName),
        total_variance_out = .assayTotalVariance(object, newAssayName),
        ## The recorded totals are the quantities the variance fractions
        ## are relative to, so they go stale exactly as a variance result
        ## does. Fingerprint both assays for the same reason.
        fingerprint        = .assayFingerprint(object, assayName),
        fingerprint_out    = .assayFingerprint(object, newAssayName),
        timestamp          = Sys.time()
    )

    if (is.null(metadata(object)$correction_history)) {
        metadata(object)$correction_history <- list()
    }

    metadata(object)$correction_history <-
        c(metadata(object)$correction_history, list(entry))

    object
}
