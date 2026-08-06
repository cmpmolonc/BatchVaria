## Registry of batch correction methods.
##
## Mirrors .varianceEngines: held in an environment rather than a switch()
## so that the set of methods is data rather than code, and another package
## can add one without editing BatchVaria. The built-ins are registered in
## .onLoad().
.correctionMethods <- new.env(parent = emptyenv())

## Arguments runCorrection() supplies to every method. 'preserve' is
## required rather than optional: a method carrying only '...' would absorb
## it silently, and an unprotected correction that was asked to protect
## something returns a plausible-looking matrix rather than an error.
.correctionMethodArgs <- c("assayMatrix", "batch", "sampleData", "preserve")

#' Register a batch correction method
#'
#' Adds a correction method to the set that \code{\link{runCorrection}} can
#' dispatch to, so that the methods are extensible rather than fixed.
#'
#' @param name Character. Name callers will pass as \code{method}.
#' @param method A function implementing the correction contract below.
#' @param overwrite Logical. Replace a method of the same name.
#'
#' @return The registered name, invisibly.
#'
#' @details
#' A correction method is a function with the signature
#'
#' \preformatted{function(assayMatrix, batch, sampleData, preserve = NULL, ...)}
#'
#' where \code{assayMatrix} is features by samples, \code{batch} is the name
#' of a single \code{sampleData} column holding the batch structure,
#' \code{sampleData} is the \code{colData} of the object, and \code{preserve}
#' is an optional character vector of \code{sampleData} column names. The
#' \code{...} is required: \code{runCorrection()} forwards method-specific
#' arguments through it, so each method must tolerate arguments intended for
#' another.
#'
#' \code{preserve} is required in the signature rather than optional, which
#' is where this contract diverges from
#' \code{\link{registerVarianceEngine}}. A method declaring only \code{...}
#' would absorb it silently and return a plausible matrix from a correction
#' that protected nothing. The variance engines have no argument whose
#' silent absorption changes the answer in that way.
#'
#' Methods build their own design matrices from \code{preserve} rather than
#' receiving one, so a design assembled for one method is never forced
#' through another.
#'
#' \strong{Direction of preserve.} \code{preserve} names variables whose
#' variation should be \emph{kept}. A method protects them, typically by
#' estimating the batch effect adjusted for them. It maps to \code{mod} in
#' \code{sva::ComBat()} and to \code{design} in
#' \code{limma::removeBatchEffect()}, which its own documentation describes
#' as the conditions "to be preserved".
#'
#' Note that \code{removeBatchEffect()} \emph{also} has an argument spelled
#' \code{covariates}, and it means the opposite: numeric covariates to be
#' adjusted for, that is, removed. An adapter that wires \code{preserve} to
#' \code{covariates} will strip exactly the signal it was asked to protect,
#' and the result will still look like a successful correction -- batch
#' variance falls and the table reads as expected. Check the direction.
#'
#' \strong{Batch.} \code{batch} is a single column name, validated by
#' \code{\link{runCorrection}}. Methods accepting more than one batch
#' variable must expose the others through \code{...}.
#'
#' \strong{Return value.} A numeric matrix whose \code{dim} and
#' \code{dimnames} are \code{identical()} to those of \code{assayMatrix}.
#' This makes sample order part of the contract rather than a convention:
#' delta tables, the PCA basis and correlation structure are all keyed on
#' sample identity, so a reordered return would misalign every downstream
#' comparison without raising an error. The check is applied to whatever a
#' method returns, because a conforming signature does not guarantee a
#' conforming result.
#'
#' Methods take a matrix and return a matrix. The \code{BatchVariaData}
#' container is never passed in, and recording the correction in the ledger
#' is \code{\link{runCorrection}}'s responsibility, not the method's.
#'
#' @seealso \code{\link{availableCorrectionMethods}},
#'   \code{\link{unregisterCorrectionMethod}}, \code{\link{runCorrection}}
#'
#' @examples
#' ## a method that centres each batch on the overall feature mean
#' centreByBatch <- function(assayMatrix, batch, sampleData,
#'                           preserve = NULL, ...) {
#'     b <- factor(sampleData[[batch]])
#'     out <- assayMatrix
#'     for (lv in levels(b)) {
#'         idx <- which(b == lv)
#'         out[, idx] <- out[, idx] - rowMeans(out[, idx, drop = FALSE])
#'     }
#'     out + rowMeans(assayMatrix)
#' }
#'
#' registerCorrectionMethod("batchCentre", centreByBatch)
#' availableCorrectionMethods()
#'
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 50)
#' bv <- runCorrection(bv, method = "batchCentre", batch = "batch")
#' SummarizedExperiment::assayNames(bv)
#'
#' unregisterCorrectionMethod("batchCentre")
#'
#' @export
registerCorrectionMethod <- function(name, method, overwrite = FALSE) {
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !nzchar(name)) {
        stop("'name' must be a single non-empty string")
    }

    if (!is.function(method)) {
        stop("'method' must be a function")
    }

    formalNames <- names(formals(method))

    if (!"..." %in% formalNames) {
        stop(
            "A correction method must accept '...': runCorrection() ",
            "forwards method-specific arguments to the method, so each ",
            "must tolerate arguments intended for another"
        )
    }

    missingArgs <- setdiff(.correctionMethodArgs, formalNames)
    if (length(missingArgs) > 0) {
        stop(
            "A correction method must accept ",
            paste(.correctionMethodArgs, collapse = ", "),
            "; '", name, "' is missing ",
            paste(missingArgs, collapse = ", ")
        )
    }

    if (!overwrite &&
        exists(name, envir = .correctionMethods, inherits = FALSE)) {
        stop(
            "A correction method named '", name, "' is already registered. ",
            "Pass overwrite = TRUE to replace it"
        )
    }

    assign(name, method, envir = .correctionMethods)

    invisible(name)
}

#' Remove a registered batch correction method
#'
#' @param name Character. Name of the method to remove.
#'
#' @return The removed name, invisibly.
#'
#' @seealso \code{\link{registerCorrectionMethod}}
#'
#' @examples
#' dummy <- function(assayMatrix, batch, sampleData, preserve = NULL, ...) {
#'     assayMatrix
#' }
#' registerCorrectionMethod("dummy", dummy)
#' unregisterCorrectionMethod("dummy")
#' availableCorrectionMethods()
#'
#' @export
unregisterCorrectionMethod <- function(name) {
    if (!exists(name, envir = .correctionMethods, inherits = FALSE)) {
        stop(
            "No correction method named '", name, "' is registered. ",
            "Registered methods: ",
            paste(availableCorrectionMethods(), collapse = ", ")
        )
    }

    rm(list = name, envir = .correctionMethods)

    invisible(name)
}

.getCorrectionMethod <- function(method) {
    if (!exists(method, envir = .correctionMethods, inherits = FALSE)) {
        stop(
            "Unknown correction method: ", method,
            ". Available methods are: ",
            paste(availableCorrectionMethods(), collapse = ", ")
        )
    }

    get(method, envir = .correctionMethods, inherits = FALSE)
}

## Use-time half of the two-stage check. Registration validates the
## signature; this validates what the method actually returned, because a
## conforming signature does not imply a conforming result.
.validateCorrectionResult <- function(corrected, assayMatrix, method) {
    if (!is.matrix(corrected) || !is.numeric(corrected)) {
        stop(
            "Correction method '", method, "' must return a numeric matrix, ",
            "not ", paste(class(corrected), collapse = "/")
        )
    }

    if (!identical(dim(corrected), dim(assayMatrix))) {
        stop(
            "Correction method '", method, "' returned a matrix of ",
            paste(dim(corrected), collapse = " x "),
            "; the input was ", paste(dim(assayMatrix), collapse = " x ")
        )
    }

    ## identical() rather than setequal(): sample order is part of the
    ## contract, and a reordered return would misalign every downstream
    ## comparison rather than fail
    if (!identical(dimnames(corrected), dimnames(assayMatrix))) {
        stop(
            "Correction method '", method, "' returned different dimnames ",
            "from the input. Feature and sample names must be preserved in ",
            "order, since downstream comparisons are keyed on them"
        )
    }

    corrected
}

## Both built-in methods want the same thing from 'preserve': a model
## matrix of the variables to protect. Shared so the two cannot drift.
.preserveDesign <- function(preserve, sampleData) {
    if (is.null(preserve)) {
        return(NULL)
    }

    stats::model.matrix(
        stats::as.formula(paste("~", paste(preserve, collapse = " + "))),
        data = sampleData
    )
}
