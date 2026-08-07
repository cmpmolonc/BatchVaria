## Staleness detection for ledger entries.
##
## A ledger entry describes a computation performed on an assay at a
## moment in time. The object it is attached to remains mutable
## afterwards: bv[, 1:6] subsets it, assay(bv, "x") <- m replaces a
## matrix, colData(bv)$batch <- NULL removes a variable. None of these
## invoke S4 validity -- neither `[`, `assay<-`, `colData<-` nor
## `metadata<-` calls validObject() -- so a validity method cannot see
## them, and the entry silently becomes a description of data that is no
## longer present.
##
## The failure is quiet and confident rather than loud: varianceTable()
## on a subset object returns percentages that look entirely reasonable
## and describe samples the object no longer contains.
##
## So the check is made at read time rather than write time. Recording a
## fingerprint alongside each entry and comparing it when the entry is
## read catches every mutation, including ones not anticipated here,
## because it asks whether the data still matches rather than whether a
## known operation occurred.

## Fingerprint of the assay a ledger entry describes.
##
## The hash covers content and dimnames both, so it detects an assay
## replaced in place with a matrix of identical dimensions -- which
## dimensions alone cannot. Dimensions are stored beside it anyway, not
## for detection but so that a mismatch can be reported in terms a reader
## can act on: "computed on 20 samples, object now has 6" locates the
## problem, where a pair of differing hashes does not.
##
## NULL for an assay that is not present: recordVariance() is exported
## and may be handed an assay name the object does not have, which is a
## separate problem and not one this function should raise.
.assayFingerprint <- function(object, assayName) {
    if (!assayName %in% SummarizedExperiment::assayNames(object)) {
        return(NULL)
    }

    mat <- SummarizedExperiment::assay(object, assayName)

    list(
        n_features = nrow(mat),
        n_samples = ncol(mat),
        hash = rlang::hash(as.matrix(mat))
    )
}

## Compare a stored fingerprint against the assay as it stands now,
## returning NULL when they agree and a human-readable reason when they
## do not.
##
## Entries written before fingerprints were recorded have none. Those are
## reported as agreeing rather than as stale: an absent fingerprint is
## missing evidence, not evidence of a mismatch, and warning about it
## would fire on every object built by an earlier version.
.fingerprintMismatch <- function(object, assayName, stored) {
    if (is.null(stored) || is.null(stored$hash)) {
        return(NULL)
    }

    current <- .assayFingerprint(object, assayName)

    if (is.null(current)) {
        return("no longer present in the object")
    }

    if (identical(current$hash, stored$hash)) {
        return(NULL)
    }

    ## Dimension changes are the common case and the informative one, so
    ## report them specifically before falling back to "content changed".
    dimNote <- character(0)

    if (!identical(current$n_samples, stored$n_samples)) {
        dimNote <- c(dimNote, paste0(
            "computed on ", stored$n_samples, " samples, object now has ",
            current$n_samples
        ))
    }

    if (!identical(current$n_features, stored$n_features)) {
        dimNote <- c(dimNote, paste0(
            "computed on ", stored$n_features, " features, object now has ",
            current$n_features
        ))
    }

    if (length(dimNote) == 0) {
        dimNote <- paste0(
            "contents have changed since the result was recorded, ",
            "though the dimensions have not"
        )
    }

    paste(dimNote, collapse = " and ")
}

## Warn once for a set of ledger entries whose assays no longer match the
## fingerprints recorded with them.
##
## A warning rather than an error, and the entries are still returned. A
## stale result remains evidence of what was computed; it has simply
## stopped being a description of the current object. Refusing to return
## it would also make bv[, 1:6] illegal, which it should not be.
.warnIfStale <- function(object, entries, context) {
    if (length(entries) == 0) {
        return(invisible(NULL))
    }

    assays <- vapply(entries, function(e) e$assay, character(1))

    reasons <- vapply(
        entries,
        function(e) {
            m <- .fingerprintMismatch(object, e$assay, e$fingerprint)
            if (is.null(m)) NA_character_ else m
        },
        character(1)
    )

    keep <- !is.na(reasons)

    if (!any(keep)) {
        return(invisible(NULL))
    }

    ## Subsetting changes every assay at once, so one reason typically
    ## covers all of them. Grouping assays under a shared reason states
    ## the fact once instead of once per assay.
    assays <- assays[keep]
    reasons <- reasons[keep]

    notes <- vapply(
        unique(reasons),
        function(r) {
            paste0(
                paste0("'", unique(assays[reasons == r]), "'", collapse = ", "),
                ": ", r
            )
        },
        character(1),
        USE.NAMES = FALSE
    )

    warning(
        context, " is reporting results recorded before the object was ",
        "modified, so they do not describe its current contents. ",
        paste(notes, collapse = "; "),
        ". Re-run profileVariance() to bring the ledger up to date",
        call. = FALSE
    )

    invisible(NULL)
}
