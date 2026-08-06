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
#' Engines differ in the formulas they accept, and each builds its own
#' design matrix. \code{variancePartition} accepts random-effects notation
#' such as \code{~ (1 | batch)}; \code{anova} models fixed effects only and
#' declines such formulas by name. Supplying a random-effects formula to
#' the default set of methods is therefore expected to profile with
#' \code{variancePartition} and report why \code{anova} did not
#' participate.
#'
#' The \code{anova} engine decomposes variance per model term using Type II
#' sums of squares, emitting one row per term plus \code{residual} and
#' \code{shared}. Type II sums of squares are order independent but do not
#' add up to the model sum of squares under an unbalanced design: variance
#' that the design cannot attribute to any single term is reported as
#' \code{shared} rather than left out, so the fractions remain a partition.
#' Under an orthogonal design \code{shared} is zero; under complete
#' confounding of two terms it absorbs everything they explain, since
#' neither is identifiable given the other. It can be negative where terms
#' explain more jointly than separately, which is reported with a warning.
#'
#' Fractions are averaged across features by default, matching how
#' \code{variancePartition} summarises, so the two engines describe the
#' typical feature and stay comparable. Pass \code{weighting = "pooled"} to
#' weight each feature by the variance it carries instead. The choice and
#' the sums-of-squares type are recorded alongside each result.
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
    # Validate the formula against colData
    #
    # Design matrices are built inside the engines that need them, not
    # here: only the anova engine consumes one, and constructing it
    # centrally forced random-effects formulas written for
    # variancePartition through stats::model.matrix(), which does not
    # understand them. Checking that the variables exist is engine
    # independent, so it stays.
    # --------------------------------
    if (!inherits(formula, "formula")) {
        stop("'formula' must be a formula")
    }

    sampleData <- SummarizedExperiment::colData(object)

    missing_vars <- setdiff(all.vars(formula), colnames(sampleData))
    if (length(missing_vars) > 0) {
        stop(
            "Formula variables not found in colData: ",
            paste(missing_vars, collapse = ", ")
        )
    }

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

