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
#' @details
#' Features with zero or non-finite variance are excluded from every
#' engine before profiling, and the number excluded is reported once per
#' assay. Such features are common in unfiltered count data and carry no
#' information for a variance decomposition; retaining them causes
#' engine-specific failures rather than engine-specific answers.
#'
#' Each assay/method combination is attempted independently. If one engine
#' fails, a warning names the assay, the method and the reason, and
#' profiling continues with the remaining combinations. An error is raised
#' only when every attempt fails, so a partial run always returns the
#' results that succeeded rather than discarding them.
#'
#' @seealso \code{\link{assayVariance}}, \code{\link{varianceHistory}}
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
    #
    # Each assay/method combination is attempted independently: one engine
    # failing must not discard the results of the engines that succeeded,
    # nor abandon the remaining assays.
    # --------------------------------
    nAttempted <- 0L
    nSucceeded <- 0L

    for (assayNameI in selected_assays) {
        assayMatrix <- SummarizedExperiment::assay(object, assayNameI)
        assayMatrix <- as.matrix(assayMatrix)

        ## Screen zero-variance features once per assay, so that every
        ## engine decomposes the same feature set and the exclusion is
        ## reported once rather than once per method.
        assayMatrix <- tryCatch(
            .dropConstantFeatures(assayMatrix, assayName = assayNameI),
            error = function(e) {
                warning(
                    "Skipping assay '", assayNameI, "': ",
                    conditionMessage(e),
                    call. = FALSE
                )
                NULL
            }
        )

        if (is.null(assayMatrix)) {
            nAttempted <- nAttempted + length(methods)
            next
        }

        for (method in methods) {
            nAttempted <- nAttempted + 1L
            engine_fun <- .getVarianceEngine(method)

            summaryDf <- tryCatch(
                {
                    out <- engine_fun(
                        assayMatrix = assayMatrix,
                        modelMatrix = modelMatrix,
                        formula = formula,
                        sampleData = sampleData,
                        ...
                    )

                    # ------------------------------------
                    # Enforce variance result contract
                    # ------------------------------------
                    requiredCols <- c("source", "term", "variance_fraction")

                    if (!is.data.frame(out)) {
                        stop("engine must return a data.frame")
                    }

                    if (!all(requiredCols %in% colnames(out))) {
                        stop(
                            "engine returned invalid output. Required ",
                            "columns: ", paste(requiredCols, collapse = ", ")
                        )
                    }

                    out
                },
                error = function(e) {
                    warning(
                        "Variance method '", method, "' failed for assay '",
                        assayNameI, "': ", conditionMessage(e),
                        call. = FALSE
                    )
                    NULL
                }
            )

            if (is.null(summaryDf)) {
                next
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

            nSucceeded <- nSucceeded + 1L
        }
    }

    if (nSucceeded == 0L) {
        stop(
            "All ", nAttempted, " variance profiling attempt",
            if (nAttempted == 1L) "" else "s", " failed; ",
            "the reason for each was reported above"
        )
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
