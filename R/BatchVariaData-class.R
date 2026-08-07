#'
#' BatchVariaData container
#'
#' An extension of SummarizedExperiment for batch-effect exploration,
#' correction, and variance decomposition analysis.
#'
#' Expression matrices (raw and corrected) are stored as assays.
#' Analysis provenance, including correction parameters and variance
#' decomposition history, is stored in metadata().
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @return A \code{BatchVariaData} object, an S4 class extending
#'   \code{SummarizedExperiment}.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' class(bv)
#' SummarizedExperiment::assayNames(bv)
#'
#' @export
setClass("BatchVariaData", contains = "SummarizedExperiment")

## Structural and referential invariants of the two ledgers.
##
## Scope, stated plainly because it is easy to over-read: S4 validity runs
## when an object is constructed and when validObject() is called on it. It
## does not run on `[`, `assay<-`, `colData<-` or `metadata<-`, none of
## which validate. So this method catches a malformed or dangling ledger
## handed to the constructor, and catches nothing that happens to an object
## afterwards. Staleness introduced by later mutation is detected at read
## time instead -- see .warnIfStale() -- because that is the only place it
## can be.
##
## Referential rules are safe under subsetting: bv[, 1:6] preserves every
## assay name and colData column, so an entry valid before a subset is
## still valid after one. What it stops describing is the data, which is a
## different question and a different mechanism.
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

setValidity("BatchVariaData", function(object) {
    md <- S4Vectors::metadata(object)
    assayNamesPresent <- SummarizedExperiment::assayNames(object)
    colDataNames <- colnames(SummarizedExperiment::colData(object))

    problems <- c(
        .validateVarianceLedger(md$variance_history, assayNamesPresent),
        .validateCorrectionLedger(
            md$correction_history, assayNamesPresent, colDataNames
        )
    )

    if (length(problems) == 0) TRUE else problems
})
