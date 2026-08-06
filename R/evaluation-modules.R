#'
#' Evaluate batch correction performance across assays
#'
#' Runs a suite of diagnostic comparisons between a baseline assay and
#' one or more corrected assays, including variance decomposition,
#' PCA structure, and sample correlation analysis.
#'
#' @param object A \code{BatchVariaData} object.
#' @param baseline Character string specifying the reference assay
#'   (default \code{"raw"}).
#' @param diagnostics Character vector of diagnostics to compute.
#'   Supported options include \code{"variance"}, \code{"pca"},
#'   and \code{"correlation"}.
#' @param assays Optional character vector of assay names to evaluate.
#'   If \code{NULL}, all assays except the baseline are used.
#' @param varianceMethod Character string specifying the variance
#'   decomposition method to use (default \code{"anova"}).
#'
#' @return A named list containing diagnostic results for each assay
#'   comparison.
#'
#' @details
#' This function acts as a high-level orchestration layer, combining
#' multiple evaluation strategies into a unified interface. By default,
#' all supported diagnostics are run across all non-baseline assays.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' res <- evaluateCorrections(bv, assays = c("raw", "raw_center"))
#' names(res)
#'
#' @export

evaluateCorrections <- function(
    object,
    baseline = "raw",
    diagnostics = c("variance", "pca", "correlation"),
    assays = NULL,
    varianceMethod = "anova"
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    available_assays <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available_assays
    }

    ## validate assays
    missing_assays <- setdiff(assays, available_assays)
    if (length(missing_assays) > 0) {
        stop("Assays not found: ", paste(missing_assays, collapse = ", "))
    }

    if (!baseline %in% assays) {
        stop("Baseline assay must be included in 'assays': ", baseline)
    }

    ## validate diagnostics
    valid_diagnostics <- c("variance", "pca", "correlation")
    invalid_diag <- setdiff(diagnostics, valid_diagnostics)

    if (length(invalid_diag) > 0) {
        stop(
            "Invalid diagnostics: ",
            paste(invalid_diag, collapse = ", "),
            ". Valid options are: ",
            paste(valid_diagnostics, collapse = ", ")
        )
    }

    results <- list()

    ## -----------------------------
    ## Variance redistribution
    ## -----------------------------
    if ("variance" %in% diagnostics) {
        var_df <- varianceDelta(
            object,
            assays = assays,
            method = varianceMethod,
            baseline = baseline
        )

        .validatePlotDf(
            var_df,
            requiredCols = c("assay", "term", "delta"),
            fnName = "evaluateCorrections:varianceDelta",
            allowNA = FALSE
        )

        results$variance_delta <- var_df
    }

    ## -----------------------------
    ## PCA comparison
    ## -----------------------------
    if ("pca" %in% diagnostics) {
        pca_res <- comparePCA(
            object,
            assays = assays
        )

        ## optional: validate structure if standardised later
        results$pca_comparison <- pca_res
    }

    ## -----------------------------
    ## Correlation preservation
    ## -----------------------------
    if ("correlation" %in% diagnostics) {
        ## ensure at least one comparison exists
        if (length(setdiff(assays, baseline)) == 0) {
            stop(
                "At least one non-baseline assay is required for correlation comparison"
            )
        }

        ## compareCorrelations() requires the baseline to be among 'assays'
        cor_df <- compareCorrelations(
            object,
            baseline = baseline,
            assays = assays
        )

        ## enforce strict contract
        .validatePlotDf(
            cor_df,
            requiredCols = c("assay", "correlation_change"),
            fnName = "evaluateCorrections:correlation_change",
            allowNA = FALSE
        )

        results$correlation_change <- cor_df
    }


    ## -----------------------------
    ## Attach metadata (optional but powerful)
    ## -----------------------------
    attr(results, "baseline") <- baseline
    attr(results, "assays") <- assays
    attr(results, "diagnostics") <- diagnostics

    results
}

#'
#' Compute variance redistribution relative to baseline
#'
#' Calculates the change in variance explained by each model term
#' for one or more assays relative to a baseline assay.
#'
#' @param object BatchVariaData object
#' @param assays Character vector of assay names (default: all assays)
#' @param method Variance estimation method
#' @param formula Model formula to compare under. Required when the chosen
#'   method has results for more than one formula, since a delta between
#'   decompositions of different models is not interpretable.
#' @param baseline Baseline assay name. When \code{NULL} (the default) it
#'   is inferred from the correction ledger, falling back to \code{"raw"}.
#'
#' @return data.frame with columns: assay, method, term, delta
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")
#' varianceDelta(bv)
#'
#' @export

varianceDelta <- function(
    object,
    assays = NULL,
    method = "anova",
    formula = NULL,
    baseline = NULL
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    ## Deltas are only meaningful between results from the same model, so
    ## resolve one formula for every assay rather than letting each fall
    ## back to its own most-recent entry.
    formulaKey <- .resolveFormulaKey(
        S4Vectors::metadata(object)$variance_history,
        method,
        formula
    )

    available_assays <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available_assays
    }

    missing_assays <- setdiff(assays, available_assays)
    if (length(missing_assays) > 0) {
        stop("Assays not found: ", paste(missing_assays, collapse = ", "))
    }

    if (!is.null(baseline) && !baseline %in% assays) {
        stop("Baseline assay not found: ", baseline)
    }

    baseline <- .resolveBaseline(object, assays, baseline)

    ## Unlike varianceTable(), this function reports nothing but deltas, so
    ## a missing baseline is fatal rather than merely limiting.
    if (is.null(baseline)) {
        stop(
            "No baseline assay could be determined. Name one with ",
            "baseline =; the assays available are ",
            paste(assays, collapse = ", ")
        )
    }

    assays_no_base <- setdiff(assays, baseline)

    if (length(assays_no_base) == 0) {
        stop("At least one non-baseline assay is required")
    }

    ## ensure method is vector
    methods <- method

    ## -----------------------------
    ## Compute per-method deltas
    ## -----------------------------
    out_all <- lapply(methods, function(m) {
        ## baseline
        base_df <- .getVarianceResult(
            object,
            assayName = baseline,
            method = m,
            formulaKey = formulaKey
        )

        if (!all(c("term", "variance_fraction") %in% colnames(base_df))) {
            stop("Invalid variance result for baseline assay (method = ", m, ")")
        }

        base_df <- base_df[, c("term", "variance_fraction")]
        colnames(base_df)[2] <- "baseline"

        ## per assay
        out_list <- lapply(assays_no_base, function(a) {
            df <- .getVarianceResult(
                object,
                assayName = a,
                method = m,
                formulaKey = formulaKey
            )

            if (!all(c("term", "variance_fraction") %in% colnames(df))) {
                stop("Invalid variance result for assay: ", a, " (method = ", m, ")")
            }

            df <- df[, c("term", "variance_fraction")]

            merged <- merge(base_df, df, by = "term", all = FALSE)

            data.frame(
                assay = a,
                method = m,
                term = merged$term,
                delta = merged$variance_fraction - merged$baseline,
                stringsAsFactors = FALSE
            )
        })

        do.call(rbind, out_list)
    })

    df_out <- do.call(rbind, out_all)

    ## -----------------------------
    ## Final validation
    ## -----------------------------
    .validatePlotDf(
        df_out,
        requiredCols = c("assay", "method", "term", "delta"),
        fnName = "varianceDelta"
    )

    df_out
}

#'
#' Compare PCA variance structure across assays
#'
#' Computes the proportion of variance explained by principal components
#' for one or more assays.
#'
#' @param object BatchVariaData object
#' @param assays Character vector of assays (default: all assays)
#' @param nPCs Number of principal components to return (default: 3)
#'
#' @return data.frame with columns: assay, component, variance
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' comparePCA(bv)
#'
#' @export

comparePCA <- function(
    object,
    assays = NULL,
    nPCs = 3
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    available_assays <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available_assays
    }

    missing_assays <- setdiff(assays, available_assays)
    if (length(missing_assays) > 0) {
        stop("Assays not found: ", paste(missing_assays, collapse = ", "))
    }

    if (nPCs < 1) {
        stop("nPCs must be >= 1")
    }

    ## -----------------------------
    ## Compute PCA per assay
    ## -----------------------------
    res_list <- lapply(assays, function(a) {
        mat <- SummarizedExperiment::assay(object, a)
        mat <- as.matrix(mat)

        p <- stats::prcomp(t(mat), scale. = TRUE)

        imp <- summary(p)$importance

        ## safe PC selection
        max_pcs <- min(nPCs, ncol(imp))
        vars <- imp[2, seq_len(max_pcs)]

        data.frame(
            assay = a,
            component = paste0("PC", seq_len(max_pcs)),
            variance = as.numeric(vars),
            stringsAsFactors = FALSE
        )
    })

    df_out <- do.call(rbind, res_list)

    ## -----------------------------
    ## Final validation
    ## -----------------------------
    .validatePlotDf(
        df_out,
        requiredCols = c("assay", "component", "variance"),
        fnName = "comparePCA"
    )

    df_out
}

#'
#' Compare sample correlation structure across assays
#'
#' Computes the mean absolute difference in sample–sample correlation
#' matrices between a baseline assay and one or more comparison assays.
#'
#' @param object BatchVariaData object
#' @param baseline Baseline assay name. When \code{NULL} (the default) it
#'   is inferred from the correction ledger, falling back to \code{"raw"}.
#' @param assays Character vector of assays to compare (default: all assays)
#'
#' @return data.frame with columns: assay, correlation_change
#'
#' @details
#' Samples are the columns of each assay, so the correlation matrices are
#' computed over columns and are \code{ncol(assay)} by \code{ncol(assay)}.
#' The reported change is averaged over the unique sample pairs (the strict
#' upper triangle): the diagonal is identical in every assay by construction,
#' and each off-diagonal pair appears twice in a symmetric matrix.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' compareCorrelations(bv)
#'
#' @export

compareCorrelations <- function(
    object,
    baseline = "raw",
    assays = NULL
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    available_assays <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available_assays
    }

    missing_assays <- setdiff(assays, available_assays)
    if (length(missing_assays) > 0) {
        stop("Assays not found: ", paste(missing_assays, collapse = ", "))
    }

    if (!baseline %in% assays) {
        stop("Baseline assay must be included in 'assays': ", baseline)
    }

    assays_no_base <- setdiff(assays, baseline)

    if (length(assays_no_base) == 0) {
        stop("At least one non-baseline assay is required")
    }

    ## -----------------------------
    ## Baseline correlation
    ## -----------------------------
    base_mat <- SummarizedExperiment::assay(object, baseline)
    base_mat <- as.matrix(base_mat)

    ## Samples are the columns of the assay, so cor() on the matrix itself
    ## gives the sample-by-sample correlation structure.
    base_cor <- stats::cor(base_mat, use = "pairwise.complete.obs")

    ## Unique sample pairs only: the diagonal is identical across assays and
    ## each off-diagonal pair is duplicated in a symmetric matrix.
    sample_pairs <- upper.tri(base_cor)

    if (!any(sample_pairs)) {
        stop("At least two samples are required to compare correlation structure")
    }

    ## -----------------------------
    ## Compare assays
    ## -----------------------------
    res_list <- lapply(assays_no_base, function(a) {
        mat <- SummarizedExperiment::assay(object, a)
        mat <- as.matrix(mat)

        ## dimension check
        if (ncol(mat) != ncol(base_mat)) {
            stop("Sample count mismatch between baseline and assay: ", a)
        }

        cor_mat <- stats::cor(mat, use = "pairwise.complete.obs")

        diff <- mean(
            abs(base_cor[sample_pairs] - cor_mat[sample_pairs]),
            na.rm = TRUE
        )

        data.frame(
            assay = a,
            correlation_change = diff,
            stringsAsFactors = FALSE
        )
    })

    df_out <- do.call(rbind, res_list)

    ## -----------------------------
    ## Final validation
    ## -----------------------------
    .validatePlotDf(
        df_out,
        requiredCols = c("assay", "correlation_change"),
        fnName = "compareCorrelations"
    )

    df_out
}
