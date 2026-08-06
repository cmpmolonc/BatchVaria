#' Record variance decomposition results
#'
#' Stores the output of a variance decomposition method in the internal
#' variance ledger of a \code{BatchVariaData} object. Each entry captures
#' the assay, model formula, method used, and resulting variance summary.
#'
#' @param object A \code{BatchVariaData} object
#' @param assayName Character. Name of the assay associated with the result
#' @param formula Model formula used for variance decomposition
#' @param result Data frame containing variance decomposition results.
#'   Must include at least \code{source}, \code{term}, and
#'   \code{variance_fraction} columns. Results built by
#'   \code{\link{newVarianceSummary}}, and so every result from a
#'   registered engine, also carry \code{n_features}; a summary assembled
#'   by hand need not, and that column is then \code{NA} for the entry in
#'   \code{\link{varianceResults}} output. The \code{NA} records absent
#'   metadata rather than a failure, and nothing downstream reads it.
#' @param method Character. Name of the variance decomposition method
#'
#' @return Updated \code{BatchVariaData} object with appended variance record
#'
#' @details
#' Results are appended to \code{metadata(object)$variance_history} as a list of
#' entries. Each entry contains:
#' \itemize{
#'   \item \code{assay} – assay name
#'   \item \code{formula} – model formula
#'   \item \code{method} – variance decomposition method
#'   \item \code{result} – variance summary data frame
#' }
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' res <- data.frame(
#'     source = "anova",
#'     term = c("model", "residual"),
#'     variance_fraction = c(0.3, 0.7)
#' )
#' bv <- recordVariance(
#'     bv,
#'     assayName = "raw",
#'     formula = ~batch,
#'     result = res,
#'     method = "anova"
#' )
#' length(varianceHistory(bv))
#'
#' @export
recordVariance <- function(
    object,
    assayName,
    formula,
    result,
    method = "variancePartition"
) {
    stopifnot(is(object, "BatchVariaData"))

    ## recordVariance() is exported, so it is a way into the ledger that
    ## does not pass through profileVariance(). Enforce the result contract
    ## here rather than only in the caller, so that nothing reaches the
    ## ledger which the summary and plotting layers cannot read.
    .validateVarianceSummary(as.data.frame(result))

    vh <- metadata(object)$variance_history

    ## Keyed on the canonical formula rather than the formula object, so
    ## that ~ batch + group and ~ group + batch are recognised as the same
    ## model. Consumers resolve on (assay, method, formula) too, so
    ## distinct formulas coexist without colliding.
    formulaKey <- .formulaKey(formula)

    if (length(.matchLedger(
        vh,
        assayName = assayName,
        method = method,
        formulaKey = formulaKey
    )) > 0) {
        warning(
            "Variance result already recorded for assay '", assayName,
            "', method '", method, "', formula ", formulaKey,
            "; the newer result supersedes it"
        )
    }

    if (is.null(vh)) {
        vh <- list()
    }

    vh[[length(vh) + 1]] <- list(
        assay = assayName,
        formula = formula,
        formula_key = formulaKey,
        method = method,
        result = result,
        timestamp = Sys.time()
    )

    metadata(object)$variance_history <- vh

    object
}

#' Create variance summary tables
#'
#' Generates a variance percentage table and a delta-versus-baseline table
#' from the variance ledger, for a single method and model formula.
#'
#' @importFrom rlang .data
#' @param object BatchVariaData object
#' @param assays Character vector of assay names to include
#' @param method Variance method to summarise. Required when the ledger
#'   holds results for more than one method.
#' @param formula Model formula to summarise. Required when the chosen
#'   method has results for more than one formula.
#' @param baseline Assay that deltas are computed against. When
#'   \code{NULL} (the default) it is inferred from the correction ledger
#'   -- the assay the corrected assays were derived from -- falling back to
#'   \code{"raw"} if present. If no baseline can be determined the percent
#'   table is still returned and \code{delta} is \code{NULL}.
#' @param formatDelta Logical; if TRUE, format delta with signs
#' @param digits Number of digits for rounding
#' @return list with \code{percent}, \code{delta} and the \code{baseline}
#'   the delta was computed against. \code{delta} is \code{NULL} when no
#'   baseline could be determined.
#'
#' @details
#' A variance table describes one decomposition. Results from different
#' engines are not commensurable -- \code{pca} reports variance along
#' principal axes while \code{anova} and \code{variancePartition} report
#' variance attributable to model terms -- so stacking them in a single
#' table produces columns that sum to well over 100\%. Likewise, the same
#' assay profiled under two formulas yields two different decompositions.
#'
#' The \code{anova} engine contributes a \code{shared} row holding
#' variance that an unbalanced design cannot attribute to any single term.
#' It is zero under an orthogonal design and positive under confounding.
#' It is negative where terms explain more jointly than separately (a
#' suppression effect), which shows as a negative percentage; the column
#' still sums to 100 because the other terms are correspondingly larger.
#'
#' \code{method} and \code{formula} therefore select exactly one
#' decomposition. Either may be omitted when the ledger is unambiguous;
#' when it is not, the error lists what is available rather than choosing
#' silently. Where an assay has been profiled repeatedly under the same
#' key, the most recent result is used.
#'
#' @seealso \code{\link{assayVariance}} for the total variance that these
#'   percentages are expressed relative to
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' varianceTable(bv)
#'
#' @export
varianceTable <- function(
    object,
    assays = NULL,
    method = NULL,
    formula = NULL,
    baseline = NULL,
    formatDelta = FALSE,
    digits = 2
) {
    stopifnot(is(object, "BatchVariaData"))

    vh <- metadata(object)$variance_history

    if (length(vh) == 0) {
        stop("No variance history found. Run profileVariance() first.")
    }

    ## Resolve to exactly one decomposition before reading any results.
    method <- .resolveMethod(vh, method)
    formulaKey <- .resolveFormulaKey(vh, method, formula)

    profiled <- unique(vapply(
        .matchLedger(vh, method = method, formulaKey = formulaKey),
        function(x) x$assay,
        character(1)
    ))

    if (is.null(assays)) {
        assays <- profiled
    } else {
        missing_assays <- setdiff(assays, profiled)
        if (length(missing_assays) > 0) {
            stop(
                "No variance results recorded for assay(s) ",
                paste(missing_assays, collapse = ", "),
                " with method '", method, "' and formula ", formulaKey
            )
        }
    }

    ## An explicitly named baseline that is absent is a caller error. A
    ## baseline that simply cannot be inferred is not: the percent table
    ## does not depend on one, so report it and omit the delta.
    if (!is.null(baseline) && !baseline %in% assays) {
        stop(
            "Baseline assay '", baseline,
            "' is not among the assays being summarised: ",
            paste(assays, collapse = ", ")
        )
    }

    baseline <- .resolveBaseline(object, assays, baseline)

    ## One result per assay: .getVarianceResult() resolves repeat
    ## profiling of the same key by recency, so each (assay, term) pair
    ## appears exactly once and the pivot below cannot produce list-columns.
    df_long <- do.call(rbind, lapply(assays, function(a) {
        res <- as.data.frame(
            .getVarianceResult(object, a, method, formulaKey)
        )
        data.frame(
            term = res$term,
            assay = a,
            percent = res$variance_fraction * 100,
            stringsAsFactors = FALSE
        )
    }))

    duplicated_keys <- duplicated(df_long[, c("term", "assay")])
    if (any(duplicated_keys)) {
        stop(
            "Variance engine '", method, "' returned repeated terms for a ",
            "single assay: ",
            paste(unique(df_long$term[duplicated_keys]), collapse = ", ")
        )
    }

    # Percent table
    percent_table <- df_long |>
        tidyr::pivot_wider(
            names_from = "assay",
            values_from = "percent"
        )

    percent_table <- percent_table |>
        dplyr::select("term", dplyr::all_of(assays))

    colnames(percent_table)[1] <- "component"

    # Delta table
    if (is.null(baseline)) {
        warning(
            "No baseline assay could be determined, so no delta table was ",
            "produced. Name one with baseline =; the assays available are ",
            paste(assays, collapse = ", "),
            call. = FALSE
        )
        delta_table <- NULL
    } else {
        baseline_vals <- percent_table[[baseline]]
        delta_table <- percent_table

        for (assayName in setdiff(assays, baseline)) {
            delta_table[[assayName]] <-
                percent_table[[assayName]] - baseline_vals
        }

        delta_table <- delta_table |>
            dplyr::select(-dplyr::all_of(baseline))

        # Optional formatting
        if (formatDelta) {
            format_fun <- function(x) {
                if (is.na(x)) {
                    return(NA_character_)
                }
                if (round(x, digits) == 0) {
                    return("0")
                }
                if (x > 0) {
                    return(paste0("+", round(x, digits)))
                }
                ## remaining case: x < 0
                paste0("-", round(abs(x), digits))
            }
            ## column 1 is 'component', so format the assay columns only
            for (i in seq_len(ncol(delta_table))[-1L]) {
                delta_table[[i]] <- vapply(
                    delta_table[[i]], format_fun, character(1)
                )
            }
        } else {
            delta_table <- delta_table |>
                dplyr::mutate(dplyr::across(-"component", \(x) round(x, digits)))
        }
    }

    percent_table <- percent_table |>
        dplyr::mutate(dplyr::across(-"component", \(x) round(x, digits)))

    list(
        percent = percent_table,
        delta = delta_table,
        baseline = baseline
    )
}

#'
#' Summarise variance change vs raw for a model term
#'
#' @param object BatchVariaData object
#' @param term Model term to evaluate (e.g. "batch")
#' @param assays Assays to include
#' @param method Variance method to summarise. Required when the ledger
#'   holds results for more than one method.
#' @param formula Model formula to summarise. Required when the chosen
#'   method has results for more than one formula.
#' @param baseline Assay that changes are computed against. Inferred from
#'   the correction ledger when \code{NULL}.
#'
#' @return data.frame of variance change for a single model term versus
#'   the baseline assay
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' varianceChange(bv, term = "residual")
#'
#' @export
varianceChange <- function(
    object,
    term,
    assays = NULL,
    method = NULL,
    formula = NULL,
    baseline = NULL
) {
    vt <- varianceTable(
        object,
        assays = assays,
        method = method,
        formula = formula,
        baseline = baseline
    )
    delta <- vt$delta

    ## Unlike the percent table, a change is meaningless without something
    ## to change relative to.
    if (is.null(delta)) {
        stop(
            "No baseline assay could be determined, so variance change ",
            "cannot be computed. Name one with baseline ="
        )
    }

    if (!term %in% delta$component) {
        stop("Term not found in variance table: ", term)
    }

    row <- delta[delta$component == term, ]

    assays <- colnames(delta)[-1]

    data.frame(
        assay = assays,
        delta_variance = as.numeric(row[1, assays])
    )
}
