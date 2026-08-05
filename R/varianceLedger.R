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
#'   \code{variance_fraction} columns
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
    method = "variancePartition") {
    vh <- metadata(object)$variance_history

    if (any(vapply(vh, function(x) {
        identical(x$assay, assayName) &&
            identical(x$formula, formula) &&
            identical(x$method, method)
    }, logical(1)))) {
        warning("Variance result already recorded for this assay/formula/method")
    }

    if (is.null(vh)) {
        vh <- list()
    }

    vh[[length(vh) + 1]] <- list(
        assay = assayName,
        formula = formula,
        method = method,
        result = result,
        timestamp = Sys.time()
    )

    metadata(object)$variance_history <- vh

    object
}

#' Create variance summary tables
#'
#' Generates variance percentage table and delta vs raw table
#' from the variance ledger.
#'
#' @importFrom rlang .data
#' @param object BatchVariaData object
#' @param assays Character vector of assay names to include
#' @param formatDelta Logical; if TRUE, format delta with arrows
#' @param digits Number of digits for rounding
#' @return list with percent table and delta table
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
    formatDelta = FALSE,
    digits = 2
) {
    stopifnot(is(object, "BatchVariaData"))

    vh <- metadata(object)$variance_history

    if (length(vh) == 0) {
        stop("No variance history found. Run profileVariance() first.")
    }

    # Build long table
    df_list <- lapply(vh, function(entry) {
        df <- as.data.frame(entry$result)
        df$assay <- entry$assay
        df
    })

    df_long <- do.call(rbind, df_list)

    if (!is.null(assays)) {
        df_long <- df_long[df_long$assay %in% assays, ]
    } else {
        assays <- unique(df_long$assay)
    }

    # Percent table
    percent_table <- df_long |>
        dplyr::select("term", "assay", "variance_fraction") |>
        dplyr::mutate(percent = .data$variance_fraction * 100) |>
        dplyr::select(-"variance_fraction") |>
        tidyr::pivot_wider(
            names_from = "assay",
            values_from = "percent"
        )

    percent_table <- percent_table |>
        dplyr::select("term", dplyr::all_of(assays))

    colnames(percent_table)[1] <- "component"

    # Delta table
    if (!"raw" %in% assays) {
        stop("Delta calculation requires 'raw' assay")
    }

    raw_vals <- percent_table$raw
    delta_table <- percent_table

    for (assayName in assays) {
        if (assayName != "raw") {
            delta_table[[assayName]] <- percent_table[[assayName]] - raw_vals
        }
    }

    delta_table <- delta_table |>
        dplyr::select(-"raw")

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
            delta_table[[i]] <- vapply(delta_table[[i]], format_fun, character(1))
        }
    } else {
        delta_table <- delta_table |>
            dplyr::mutate(dplyr::across(-"component", \(x) round(x, digits)))
    }

    percent_table <- percent_table |>
        dplyr::mutate(dplyr::across(-"component", \(x) round(x, digits)))

    list(
        percent = percent_table,
        delta   = delta_table
    )
}

#'
#' Summarise variance change vs raw for a model term
#'
#' @param object BatchVariaData object
#' @param term Model term to evaluate (e.g. "batch")
#' @param assays Assays to include
#'
#' @return data.frame of variance change for a single model term vs raw
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' varianceChange(bv, term = "residual")
#'
#' @export
varianceChange <- function(object, term, assays = NULL) {
    vt <- varianceTable(object, assays = assays)
    delta <- vt$delta

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
