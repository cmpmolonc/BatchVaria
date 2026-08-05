# Internal: validate plotting data frame
.validatePlotDf <- function(
    df,
    requiredCols,
    fnName = "plot function",
    allowNA = FALSE
) {
    if (!is.data.frame(df)) {
        stop(fnName, ": input must be a data.frame")
    }

    ## check required columns
    missing_cols <- setdiff(requiredCols, colnames(df))

    if (length(missing_cols) > 0) {
        stop(
            fnName, ": missing required columns: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    ## optional NA check
    if (!allowNA) {
        na_cols <- requiredCols[
            vapply(df[requiredCols], function(x) any(is.na(x)), logical(1))
        ]

        if (length(na_cols) > 0) {
            stop(
                fnName, ": NA values found in columns: ",
                paste(na_cols, collapse = ", ")
            )
        }
    }

    ## zero-row guard
    if (nrow(df) == 0) {
        stop(fnName, ": no data available for plotting")
    }

    invisible(df)
}
