.getVarianceResult <- function(object, assayName, method) {
    vh <- S4Vectors::metadata(object)$variance_history

    if (is.null(vh) || length(vh) == 0) {
        stop("No variance history found in object.")
    }

    matches <- vapply(
        vh,
        function(x) {
            identical(x$assay, assayName) && identical(x$method, method)
        },
        logical(1)
    )

    if (!any(matches)) {
        stop(
            "No variance result found for assay '",
            assayName, "' with method '", method, "'."
        )
    }

    res <- vh[matches]

    ## return most recent result if multiple matches
    if (length(res) > 1) {
        timestamps <- vapply(res, function(x) x$timestamp, numeric(1))
        res <- res[[which.max(timestamps)]]
    } else {
        res <- res[[1]]
    }

    res$result
}
