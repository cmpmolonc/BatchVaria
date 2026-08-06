## Canonical key for a model formula.
##
## Term order is not part of a model's identity: ~ batch + group and
## ~ group + batch describe the same model. The formula objects are not
## identical(), though, so a ledger keyed on the formula object alone
## records them as two distinct entries and then cannot tell that the
## results describe the same fit.
.formulaKey <- function(formula) {
    if (is.null(formula)) {
        return(NA_character_)
    }

    if (!inherits(formula, "formula")) {
        return(as.character(formula)[1L])
    }

    tt <- stats::terms(formula)

    ## Sorting the term labels is not enough: an interaction label keeps
    ## the order it was written in, so batch * group yields "batch:group"
    ## while group * batch yields "group:batch". Canonicalise within each
    ## label before sorting across them.
    labels <- vapply(
        strsplit(attr(tt, "term.labels"), ":", fixed = TRUE),
        function(parts) paste(sort(parts), collapse = ":"),
        character(1)
    )
    labels <- sort(labels)

    prefix <- if (attr(tt, "intercept") == 0) "~0 + " else "~"

    paste0(
        prefix,
        if (length(labels) > 0) paste(labels, collapse = " + ") else "1"
    )
}

## Formula key of a ledger entry, tolerating entries written before the
## key was stored alongside the formula.
.entryFormulaKey <- function(entry) {
    if (!is.null(entry$formula_key)) {
        return(entry$formula_key)
    }
    .formulaKey(entry$formula)
}

## Ledger entries matching an (assay, method, formula) key. NULL means
## "any" for that component.
.matchLedger <- function(vh, assayName = NULL, method = NULL, formulaKey = NULL) {
    if (is.null(vh) || length(vh) == 0) {
        return(list())
    }

    keep <- vapply(
        vh,
        function(x) {
            (is.null(assayName) || identical(x$assay, assayName)) &&
                (is.null(method) || identical(x$method, method)) &&
                (is.null(formulaKey) ||
                    identical(.entryFormulaKey(x), formulaKey))
        },
        logical(1)
    )

    vh[keep]
}

## Resolve which method a summary should describe. A ledger holding more
## than one method cannot be flattened into a single table: the engines
## answer different questions and their terms are not commensurable.
.resolveMethod <- function(vh, method) {
    available <- unique(vapply(vh, function(x) x$method, character(1)))

    if (is.null(method)) {
        if (length(available) == 1L) {
            return(available)
        }
        stop(
            "The variance ledger holds results for ", length(available),
            " methods (", paste(available, collapse = ", "),
            "). Their terms are not comparable, so specify method = ",
            "to choose one"
        )
    }

    if (!method %in% available) {
        stop(
            "No variance results recorded for method '", method,
            "'. Available: ", paste(available, collapse = ", ")
        )
    }

    method
}

## Resolve which formula a summary should describe, within a method.
.resolveFormulaKey <- function(vh, method, formula) {
    entries <- .matchLedger(vh, method = method)

    ## Report the missing method rather than "0 formulas" when a caller
    ## names a method the ledger has no results for.
    if (length(entries) == 0) {
        stop("No variance results recorded for method '", method, "'")
    }

    available <- unique(vapply(entries, .entryFormulaKey, character(1)))

    if (is.null(formula)) {
        if (length(available) == 1L) {
            return(available)
        }
        stop(
            "Method '", method, "' has results for ", length(available),
            " formulas (", paste(available, collapse = "; "),
            "). Specify formula = to choose one"
        )
    }

    key <- .formulaKey(formula)

    if (!key %in% available) {
        stop(
            "No variance results recorded for method '", method,
            "' with formula ", key,
            ". Available: ", paste(available, collapse = "; ")
        )
    }

    key
}

## Choose the assay that changes are measured against.
##
## "raw" is a convention of exampleBatchVaria(), not a requirement of the
## container: an object built from a user's own SummarizedExperiment may
## have no assay by that name. The correction ledger already records which
## assay each corrected assay was derived from, so prefer that lineage and
## fall back to the convention only when it is present.
##
## Returns NULL when no baseline can be determined, which callers treat as
## "report composition, but no delta" rather than as an error.
.resolveBaseline <- function(object, assays, baseline = NULL) {
    if (!is.null(baseline)) {
        return(baseline)
    }

    ch <- S4Vectors::metadata(object)$correction_history

    if (!is.null(ch) && length(ch) > 0) {
        parents <- vapply(ch, function(x) x$assay_in, character(1))
        derived <- vapply(ch, function(x) x$assay_out, character(1))

        ## parents of the assays being summarised, that are themselves
        ## among them -- a single common ancestor is unambiguous
        candidates <- unique(parents[derived %in% assays])
        candidates <- intersect(candidates, assays)

        if (length(candidates) == 1L) {
            return(candidates)
        }
    }

    if ("raw" %in% assays) {
        return("raw")
    }

    NULL
}

.getVarianceResult <- function(object, assayName, method, formulaKey = NULL) {
    vh <- S4Vectors::metadata(object)$variance_history

    if (is.null(vh) || length(vh) == 0) {
        stop("No variance history found in object.")
    }

    matches <- .matchLedger(
        vh,
        assayName = assayName,
        method = method,
        formulaKey = formulaKey
    )

    if (length(matches) == 0) {
        formulaNote <- if (is.null(formulaKey)) {
            "."
        } else {
            paste0(" and formula ", formulaKey, ".")
        }

        stop(
            "No variance result found for assay '", assayName,
            "' with method '", method, "'", formulaNote
        )
    }

    ## The ledger is append-only, so re-profiling the same assay, method
    ## and formula leaves more than one match. Use the most recent.
    ##
    ## which.max() alone is not enough: it returns the *first* maximum, so
    ## two entries sharing a timestamp would resolve to the stale one.
    ## Since entries are appended in order, the last of the tied maxima is
    ## the newest.
    if (length(matches) > 1) {
        timestamps <- vapply(
            matches,
            function(x) as.numeric(x$timestamp),
            numeric(1)
        )
        newest <- max(which(timestamps == max(timestamps)))
        res <- matches[[newest]]
    } else {
        res <- matches[[1]]
    }

    res$result
}
