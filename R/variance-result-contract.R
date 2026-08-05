# Internal contract for variance decomposition results
# This file defines the canonical structure expected by
# plotting, storage, and downstream analysis layers.
# Returns data frame with each row corresponding to model term
# validated via .validateVarianceSummary

.validateVarianceSummary <- function(summaryDf) {
    required <- c("source", "term", "variance_fraction")

    if (!is.data.frame(summaryDf)) {
        stop("variance summary must be a data.frame")
    }

    missing <- setdiff(required, colnames(summaryDf))
    if (length(missing) > 0) {
        stop(
            "variance summary missing required columns: ",
            paste(missing, collapse = ", ")
        )
    }

    if (!is.numeric(summaryDf$variance_fraction)) {
        stop("variance_fraction must be numeric")
    }

    if (any(summaryDf$variance_fraction < 0, na.rm = TRUE)) {
        stop("variance_fraction must be non-negative")
    }

    invisible(TRUE)
}

# variance result constructor
.newVarianceSummary <- function(
    source,
    term,
    varianceFraction,
    level = NULL,
    n = NULL
) {
    ## 'variance_fraction' is the documented column name in the result
    ## contract, so it stays snake_case even though the argument does not
    df <- data.frame(
        source = source,
        term = term,
        variance_fraction = varianceFraction,
        stringsAsFactors = FALSE
    )

    if (!is.null(level)) df$level <- level
    if (!is.null(n)) df$n <- n

    .validateVarianceSummary(df)
    df
}
