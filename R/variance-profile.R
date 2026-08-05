#' Profile sources of variance in expression data
#'
#' @param object BatchVariaData object
#' @param formula Model formula specifying covariates
#' @param methods Character vector of variance engines to apply
#' @param assayName Name of a single assay to profile. A convenience
#'   alternative to \code{assays} when only one assay is of interest;
#'   \code{assays} takes precedence if both are supplied.
#' @param assays Character vector of assays to profile (default: all assays)
#' @param ... Additional arguments passed to variance engines
#'
#' @return Updated BatchVariaData object
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' length(varianceHistory(bv))
#'
#' @export

profileVariance <- function(
    object,
    formula,
    methods = availableVarianceMethods(),
    assayName = NULL,
    assays = NULL,
    ...
) {
    stopifnot(is(object, "BatchVariaData"))

    available_assays <- SummarizedExperiment::assayNames(object)

    # --------------------------------
    # Resolve assays to run
    # --------------------------------
    if (!is.null(assays) && !is.null(assayName)) {
        warning("Both 'assays' and 'assayName' provided; using 'assays'")
    }

    if (!is.null(assays)) {
        selected_assays <- assays
    } else if (!is.null(assayName)) {
        selected_assays <- assayName
    } else {
        selected_assays <- available_assays
    }

    # Validate assays
    missing_assays <- setdiff(selected_assays, available_assays)
    if (length(missing_assays) > 0) {
        stop("Assays not found: ", paste(missing_assays, collapse = ", "))
    }

    # --------------------------------
    # Prepare model matrix once
    # --------------------------------
    sampleData <- SummarizedExperiment::colData(object)
    df <- as.data.frame(sampleData)
    modelMatrix <- stats::model.matrix(formula, data = df)

    # --------------------------------
    # Loop over assays
    # --------------------------------
    for (assayNameI in selected_assays) {
        assayMatrix <- SummarizedExperiment::assay(object, assayNameI)
        assayMatrix <- as.matrix(assayMatrix)

        for (method in methods) {
            engine_fun <- .getVarianceEngine(method)

            summaryDf <- engine_fun(
                assayMatrix        = assayMatrix,
                modelMatrix = modelMatrix,
                formula      = formula,
                sampleData     = sampleData,
                ...
            )

            # --------------------------------
            # Enforce variance result contract
            # --------------------------------
            requiredCols <- c("source", "term", "variance_fraction")

            if (!is.data.frame(summaryDf)) {
                stop("Variance engine '", method, "' must return a data.frame")
            }

            if (!all(requiredCols %in% colnames(summaryDf))) {
                stop(
                    "Variance engine '", method,
                    "' returned invalid output. Required columns: ",
                    paste(requiredCols, collapse = ", ")
                )
            }

            # --------------------------------
            # Record in variance ledger
            # --------------------------------
            object <- recordVariance(
                object,
                assayName = assayNameI,
                formula = formula,
                result = summaryDf,
                method = method
            )
        }
    }

    invisible(object)
}


.getVarianceEngine <- function(method) {
    if (!method %in% availableVarianceMethods()) {
        stop(
            "Unknown variance method: ", method,
            ". Available methods are: ",
            paste(availableVarianceMethods(), collapse = ", ")
        )
    }

    switch(method,
        pca               = .computePCAVariance,
        anova             = .computeAnovaVariance,
        variancePartition = .computeVariancePartitionVariance
    )
}
