## Structural and referential invariants of the two ledgers, checked at
## the point of use.

## Fields present, and each of the naming fields a usable name.
.validateLedgerEntry <- function(entry, i, ledger, required) {
    prefix <- paste0(ledger, "[[", i, "]]")

    if (!is.list(entry)) {
        return(paste0(prefix, " is not a list"))
    }

    missing <- setdiff(required, names(entry))

    if (length(missing) > 0) {
        return(paste0(
            prefix, " is missing required field(s): ",
            paste(missing, collapse = ", ")
        ))
    }

    ## Every required field except 'result' names something: an assay, a
    ## method, a colData column. A vector of length != 1 cannot be a name.
    scalars <- setdiff(required, "result")

    bad <- scalars[!vapply(
        entry[scalars],
        function(x) is.character(x) && length(x) == 1L && !is.na(x),
        logical(1)
    )]

    if (length(bad) > 0) {
        return(paste0(
            prefix, " field(s) must each be a single character string: ",
            paste(bad, collapse = ", ")
        ))
    }

    NULL
}

## Structural and referential rules for the variance ledger.
##
## An absent ledger is valid: an object that has not been profiled is a
## perfectly good input to profileVariance().
.validateVarianceLedger <- function(vh, assayNamesPresent) {
    if (is.null(vh)) {
        return(character(0))
    }

    if (!is.list(vh)) {
        return("metadata()$variance_history is not a list")
    }

    problems <- character(0)

    for (i in seq_along(vh)) {
        msg <- .validateLedgerEntry(
            vh[[i]], i, "variance_history", c("assay", "method", "result")
        )

        if (!is.null(msg)) {
            problems <- c(problems, msg)
            next
        }

        if (!is.data.frame(vh[[i]]$result)) {
            problems <- c(problems, paste0(
                "variance_history[[", i, "]]$result is not a data.frame"
            ))
        } else if (!all(c("term", "variance_fraction") %in%
            colnames(vh[[i]]$result))) {
            problems <- c(problems, paste0(
                "variance_history[[", i, "]]$result lacks the required ",
                "columns term and variance_fraction"
            ))
        }

        if (!vh[[i]]$assay %in% assayNamesPresent) {
            problems <- c(problems, paste0(
                "variance_history[[", i, "]] refers to assay '",
                vh[[i]]$assay, "', which the object does not have"
            ))
        }
    }

    problems
}

## Structural and referential rules for the correction ledger.
.validateCorrectionLedger <- function(ch, assayNamesPresent, colDataNames) {
    if (is.null(ch)) {
        return(character(0))
    }

    if (!is.list(ch)) {
        return("metadata()$correction_history is not a list")
    }

    problems <- character(0)

    for (i in seq_along(ch)) {
        msg <- .validateLedgerEntry(
            ch[[i]], i, "correction_history",
            c("method", "assay_in", "assay_out", "batch")
        )

        if (!is.null(msg)) {
            problems <- c(problems, msg)
            next
        }

        dangling <- setdiff(
            c(ch[[i]]$assay_in, ch[[i]]$assay_out), assayNamesPresent
        )

        if (length(dangling) > 0) {
            problems <- c(problems, paste0(
                "correction_history[[", i, "]] refers to assay(s) ",
                paste(dangling, collapse = ", "),
                ", which the object does not have"
            ))
        }

        ## batch and preserve name the colData columns the correction was
        ## defined by. Without them the entry cannot be interpreted.
        missingCols <- setdiff(
            c(ch[[i]]$batch, ch[[i]]$preserve), colDataNames
        )

        if (length(missingCols) > 0) {
            problems <- c(problems, paste0(
                "correction_history[[", i, "]] refers to colData column(s) ",
                paste(missingCols, collapse = ", "),
                ", which the object does not have"
            ))
        }
    }

    problems
}


#' Requirements BatchVaria places on a SummarizedExperiment
#'
#' Every public function in the package opens by calling the internal
#' \code{.check_se()}, which is the whole of what BatchVaria requires of
#' an object. There is no container class: an object is acceptable if it
#' is a \code{SummarizedExperiment} whose ledgers, if it has any, are
#' well formed and refer to things it still has.
#'
#' The requirements are:
#'
#' \itemize{
#'   \item The object is a \code{SummarizedExperiment}. Expression
#'     matrices are assays; sample annotation is \code{colData}.
#'   \item \code{metadata(object)$correction_history}, if present, is a
#'     list of entries each naming a \code{method}, an \code{assay_in},
#'     an \code{assay_out} and a \code{batch}. The assays it names must
#'     exist, and \code{batch} and \code{preserve} must name
#'     \code{colData} columns.
#'   \item \code{metadata(object)$variance_history}, if present, is a
#'     list of entries each naming an \code{assay} and a \code{method}
#'     and carrying a \code{result} data.frame with \code{term} and
#'     \code{variance_fraction} columns. The assay it names must exist.
#' }
#'
#' Both ledgers are optional. A plain \code{SummarizedExperiment} with
#' neither is a valid input to \code{\link{runCorrection}} and
#' \code{\link{profileVariance}}, which is how one is normally built:
#' the ledgers are what those functions write.
#'
#' @section Why not a class:
#' These rules were previously a \code{setValidity()} method on a
#' \code{BatchVariaData} class. S4 validity runs at construction and on
#' an explicit \code{validObject()}, and at no other point -- notably
#' not on \code{[}, \code{assay<-}, \code{colData<-} or
#' \code{metadata<-}. Checking at the point of use sees the object as it
#' is when a function is asked to work with it, which is the only moment
#' the answer matters.
#'
#' @return No return value. This page documents the requirements
#'   BatchVaria places on a \code{SummarizedExperiment} rather than a
#'   function; \code{\link{exampleBatchVaria}} returns an object that
#'   meets them.
#'
#' @name BatchVaria-requirements
#' @seealso \code{\link{BatchVaria}} for an overview of the package,
#'   \code{\link{provenance}} for reading the ledgers, and
#'   \code{\link{exampleBatchVaria}} for an object that satisfies them.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' class(bv)
#' SummarizedExperiment::assayNames(bv)
#'
#' ## an ordinary SummarizedExperiment is accepted; the ledgers are
#' ## written by the functions that need them
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(raw = SummarizedExperiment::assay(bv, "raw")),
#'     colData = SummarizedExperiment::colData(bv)
#' )
#' assayVariance(se)
#'
NULL


## The single entry-point guard. Returns the object invisibly so it can
## be used for its error alone or in a pipeline.
.check_se <- function(object, arg = "object") {
    if (!methods::is(object, "SummarizedExperiment")) {
        stop(arg, " must be a SummarizedExperiment", call. = FALSE)
    }

    md <- S4Vectors::metadata(object)

    problems <- c(
        .validateVarianceLedger(
            md$variance_history, SummarizedExperiment::assayNames(object)
        ),
        .validateCorrectionLedger(
            md$correction_history,
            SummarizedExperiment::assayNames(object),
            colnames(SummarizedExperiment::colData(object))
        )
    )

    if (length(problems) > 0) {
        stop(
            arg, " does not meet BatchVaria's requirements:\n",
            paste0("  - ", problems, collapse = "\n"),
            call. = FALSE
        )
    }

    invisible(object)
}
