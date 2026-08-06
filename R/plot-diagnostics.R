#'
#' Plot PCA of expression data
#'
#' Visualise global sample structure using Principal Component Analysis (PCA).
#' This function extracts an expression assay from a \code{BatchVariaData}
#' object, performs PCA on the sample-wise expression matrix, and plots the
#' first two principal components. Samples can optionally be coloured by a
#' variable in \code{colData}.
#'
#' @param object A \code{BatchVariaData} object.
#' @param assays assays to use for PCA (default \code{"raw"}).
#' @param colourBy Optional column name from \code{colData(object)} used to
#' colour samples.
#' @param reference Assay whose PCA basis every panel is projected onto.
#'   Inferred from the correction lineage when \code{NULL}.
#' @param ... Further arguments, required by the generic; unused.
#'
#' @return A list of  \code{ggplot2} objects (one per assay).
#'
#' @details
#' PCA is computed on the **transposed expression matrix** so that samples
#' represent observations and genes represent variables. The first two
#' principal components are plotted.
#'
#' All panels share one basis, fitted on \code{reference} and used to
#' project every assay. Fitting a PCA separately per panel would make
#' "PC1" a different direction in each, so panels placed side by side
#' could not be compared even though their axes carry the same labels.
#' See \code{\link{comparePCA}} for the shared-basis convention used
#' throughout the package.
#'
#' If \code{colourBy} is specified, the corresponding variable from
#' \code{colData(object)} is used to colour points in the plot.
#'
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#'
#' plotPCA(
#'     bv,
#'     assays = "raw",
#'     colourBy = "batch"
#' )
#'
#' @importFrom BiocGenerics plotPCA
#' @export
setMethod("plotPCA", "BatchVariaData", function(
    object,
    assays = "raw",
    colourBy = "batch",
    reference = NULL,
    ...
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    missingAssays <- setdiff(assays, SummarizedExperiment::assayNames(object))
    if (length(missingAssays) > 0) {
        stop("Assays not found: ", paste(missingAssays, collapse = ", "))
    }

    if (!colourBy %in% colnames(SummarizedExperiment::colData(object))) {
        stop("colourBy not found in colData: ", colourBy)
    }

    if (ncol(object) < 2) {
        stop("At least two samples are required for PCA")
    }

    ## Panels sit side by side and invite comparison, so they share one
    ## basis: a per-panel fit would label different directions "PC1".
    reference <- .resolveBaseline(object, assays, reference)
    if (is.null(reference)) {
        reference <- assays[1]
    }
    basis <- .referenceBasis(object, reference)

    plots <- lapply(assays, function(a) {
        scores <- .projectOntoBasis(object, a, basis)

        df <- data.frame(
            PC1 = scores[, 1],
            PC2 = scores[, 2],
            colour = SummarizedExperiment::colData(object)[[colourBy]]
        )

        ggplot2::ggplot(df, ggplot2::aes(
            x = .data$PC1, y = .data$PC2,
            colour = .data$colour
        )) +
            ggplot2::geom_point(size = 2) +
            ggplot2::theme_minimal() +
            ggplot2::ggtitle(paste("PCA:", a))
    })

    plots
})

#'
#' Plot variance composition across model terms
#'
#' Visualise the proportion of expression variance attributed to each model
#' term for a given assay. Variance fractions are obtained from the
#' \code{variance_history} ledger recorded by \code{profileVariance()}.
#'
#' The result is displayed as a stacked bar chart showing the relative
#' contribution of batch effects, biological covariates, and residual
#' variance.
#'
#' @param object A \code{BatchVariaData} object.
#' @param assays Character vector of assays to plot. Defaults to every assay
#'   profiled with \code{method}.
#' @param method Variance engine used to compute the results (default
#' \code{"anova"}).
#' @param formula Model formula to plot. Required when the chosen method
#'   has results for more than one formula, since a single chart can only
#'   show one decomposition per assay.
#'
#' @return A \code{ggplot2} stacked bar chart representing the variance
#' fractions attributed to each term in the variance model.
#'
#' @details
#' This function retrieves variance attribution results from
#' \code{metadata(object)$variance_history}. One bar is drawn per assay and
#' the segments within it are the model terms, so each bar sums to 100 per
#' cent. No baseline is required: this is an absolute view of where variance
#' sits, whereas \code{\link{plotVarianceDelta}} shows the signed change in
#' those fractions relative to a baseline assay.
#'
#' Where an assay has been profiled more than once with the same method, the
#' most recent record is used.
#'
#' Stacked bars assume non-negative segments. The \code{anova} engine's
#' \code{shared} term is negative under suppression, where model terms
#' explain more jointly than separately; the bar then extends below zero
#' and the positive segments overshoot 100 per cent. This is warned about,
#' and the numbers should be read from \code{\link{varianceTable}}.
#'
#' @seealso \code{\link{plotVarianceDelta}} for the change-vs-baseline view.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' plotVarianceComposition(bv)
#'
#' @export

plotVarianceComposition <- function(
    object,
    assays = NULL,
    method = "anova",
    formula = NULL
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    vh <- metadata(object)$variance_history

    if (length(vh) == 0) {
        stop("No variance history found. Run profileVariance() first.")
    }

    ## A composition plot shows one decomposition per assay, so resolve a
    ## single formula rather than letting each assay fall back to its own
    ## most-recent entry.
    formulaKey <- .resolveFormulaKey(vh, method, formula)

    profiled <- unique(vapply(
        .matchLedger(vh, method = method, formulaKey = formulaKey),
        function(x) x$assay,
        character(1)
    ))

    if (length(profiled) == 0) {
        stop("No variance results recorded for method '", method, "'")
    }

    if (is.null(assays)) {
        assays <- profiled
    } else {
        missingAssays <- setdiff(assays, profiled)
        if (length(missingAssays) > 0) {
            stop(
                "No variance results recorded for assay(s) ",
                paste(missingAssays, collapse = ", "),
                " with method '", method, "' and formula ", formulaKey
            )
        }
    }

    ## -----------------------------
    ## Absolute composition per assay
    ## -----------------------------
    parts <- lapply(assays, function(a) {
        res <- .getVarianceResult(
            object,
            assayName = a,
            method = method,
            formulaKey = formulaKey
        )

        data.frame(
            assay = a,
            term = res$term,
            percent = res$variance_fraction * 100,
            stringsAsFactors = FALSE
        )
    })

    df <- do.call(rbind, parts)

    .validatePlotDf(
        df,
        requiredCols = c("assay", "term", "percent"),
        fnName = "plotVarianceComposition"
    )

    ## A stacked bar assumes non-negative segments. The 'shared' term can
    ## be negative under suppression, and geom_col() then draws below the
    ## axis while the positive segments overshoot 100 per cent, so the bar
    ## no longer reads as a composition even though it still sums to one.
    if (any(df$percent < 0)) {
        negatives <- unique(df$term[df$percent < 0])
        warning(
            "Negative variance percentages for term(s) ",
            paste(negatives, collapse = ", "),
            ". The stacked bars extend below zero and above 100 per cent; ",
            "read the values from varianceTable() rather than the chart",
            call. = FALSE
        )
    }

    ## keep the requested assay order stable across corrections
    df$assay <- factor(df$assay, levels = assays)

    ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = .data$assay,
            y = .data$percent,
            fill = .data$term
        )
    ) +
        ggplot2::geom_col(position = "stack") +
        ggplot2::theme_minimal() +
        ggplot2::ylab("Variance explained (%)") +
        ggplot2::xlab("Assay") +
        ggplot2::labs(fill = "Model term") +
        ggplot2::ggtitle(paste0("Variance composition (", method, ")"))
}


#'
#' Plot correlation preservation after correction
#'
#' Visualises the mean absolute change in sample–sample correlation
#' structure across assays following batch correction.
#'
#' @param evalResult Output from evaluateCorrections()
#'
#' @return ggplot object
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")
#' res <- evaluateCorrections(bv, baseline = "raw")
#' plotCorrelationChange(res)
#'
#' @export

plotCorrelationChange <- function(evalResult) {
    ## -----------------------------
    ## Defensive checks
    ## -----------------------------
    if (is.null(evalResult$correlation_change)) {
        stop(
            "plotCorrelationChange(): 'correlation_change' not found",
            "Did you run evaluateCorrections()?"
        )
    }

    df <- evalResult$correlation_change

    ## -----------------------------
    ## Validate plotting contract
    ## -----------------------------
    .validatePlotDf(
        df,
        requiredCols = c("assay", "correlation_change"),
        fnName = "plotCorrelationChange"
    )

    ## enforce factor ordering (stable plotting)
    df$assay <- factor(df$assay, levels = unique(df$assay))

    ## -----------------------------
    ## Plot
    ## -----------------------------
    ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = .data$assay,
            y = .data$correlation_change
        )
    ) +
        ggplot2::geom_col() +
        ggplot2::theme_minimal() +
        ggplot2::ylab("Mean absolute correlation change") +
        ggplot2::xlab("Correction assay") +
        ggplot2::ggtitle("Sample relationship preservation")
}

#'
#' Plot variance change relative to baseline
#'
#' @param evalResult output from evaluateCorrections()
#'
#' @return ggplot object
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")
#' evalResult <- evaluateCorrections(bv)
#' plotVarianceDelta(evalResult)
#' @export

plotVarianceDelta <- function(evalResult) {
    df <- evalResult$variance_delta

    .validatePlotDf(
        df,
        requiredCols = c("assay", "method", "term", "delta"),
        fnName = "plotVarianceDelta"
    )

    ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = .data$term,
            y = .data$delta,
            fill = .data$assay
        )
    ) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::facet_wrap(~method) +
        ggplot2::theme_minimal() +
        ggplot2::ylab("Variance change vs baseline") +
        ggplot2::xlab("Model term") +
        ggplot2::ggtitle("Variance redistribution after correction")
}

#'
#' Plot sample distance heatmap
#'
#' Computes pairwise sample distances and visualises them as a heatmap.
#'
#' @importFrom tidyselect all_of
#' @param object BatchVariaData object
#' @param assayName Character string specifying assay name
#'   (default \code{"raw"})
#'
#' @return A \code{ggplot} object
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' plotSampleDistance(bv, assayName = "raw")
#' @export
plotSampleDistance <- function(object, assayName = "raw") {
    ## Validate object
    stopifnot(is(object, "BatchVariaData"))

    if (!assayName %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay not found: ", assayName)
    }

    ## Extract assay
    mat <- SummarizedExperiment::assay(object, assayName)
    mat <- as.matrix(mat)

    if (!is.matrix(mat)) {
        stop("Assay must be coercible to a matrix")
    }

    ## Ensure sample names exist
    if (is.null(colnames(mat))) {
        colnames(mat) <- paste0("Sample", seq_len(ncol(mat)))
    }

    ## Compute distance matrix
    dm <- as.matrix(stats::dist(t(mat)))

    if (!is.matrix(dm)) {
        stop("Distance computation failed: expected matrix output")
    }

    ## Convert to tidy format (safe)
    df <- as.data.frame(dm)
    df$sample1 <- rownames(df)

    df_long <- tidyr::pivot_longer(
        df,
        cols = -all_of("sample1"),
        names_to = "sample2",
        values_to = "distance"
    )

    ## Remove self-distances (diagonal)
    df_long <- df_long[df_long$sample1 != df_long$sample2, ]

    ## Enforce consistent ordering (important for comparisons)
    sample_levels <- colnames(mat)
    df_long$sample1 <- factor(df_long$sample1, levels = sample_levels)
    df_long$sample2 <- factor(df_long$sample2, levels = sample_levels)

    ## Plot
    ggplot2::ggplot(
        df_long,
        ggplot2::aes(
            x = .data$sample1,
            y = .data$sample2,
            fill = .data$distance
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(name = "Distance") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            title = paste("Sample distance:", assayName),
            x = NULL,
            y = NULL
        )
}

#'
#' Radar plot of variance redistribution
#'
#' Visualises change in variance contribution (relative to baseline)
#' across model terms using a polar (radar) representation.
#'
#' @param object BatchVariaData object
#' @param assays Character vector of assays (default: all)
#' @param method Character vector of variance methods
#' @param baseline Baseline assay
#'
#' @return ggplot object
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, methods = "anova")
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' bv <- profileVariance(bv, ~batch, assayName = "raw_combat", methods = "anova")
#' plotVarianceRadar(bv)
#'
#' @export
plotVarianceRadar <- function(
    object,
    assays = NULL,
    method = "anova",
    baseline = "raw"
) {
    ## -----------------------------
    ## Compute variance deltas
    ## -----------------------------
    vd <- varianceDelta(
        object,
        assays = assays,
        method = method,
        baseline = baseline
    )

    ## -----------------------------
    ## Validate structure
    ## -----------------------------
    .validatePlotDf(
        vd,
        requiredCols = c("assay", "method", "term", "delta"),
        fnName = "plotVarianceRadar"
    )

    ## -----------------------------
    ## Ensure meaningful comparison
    ## -----------------------------
    if (length(unique(vd$assay)) < 1) {
        stop("No comparison assays available for radar plot")
    }

    ## -----------------------------
    ## Plot
    ## -----------------------------
    ggplot2::ggplot(
        vd,
        ggplot2::aes(
            x = .data$term,
            y = .data$delta,
            group = .data$assay,
            colour = .data$assay
        )
    ) +
        ggplot2::geom_polygon(fill = NA) +
        ggplot2::geom_point() +
        ggplot2::coord_polar() +
        ggplot2::facet_wrap(~method) +
        ggplot2::theme_minimal() +
        ggplot2::ylab("Variance change vs baseline") +
        ggplot2::xlab("Model term") +
        ggplot2::ggtitle("Variance redistribution across corrections")
}

#'
#' Plot batch mixing entropy
#'
#' Computes local neighbourhood entropy based on k-nearest neighbours
#' in PCA space to assess batch mixing.
#'
#' @param object BatchVariaData object
#' @param assayName Assay name
#' @param batchVar Column in colData representing batch
#' @param k Number of neighbours
#'
#' @return ggplot object
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' plotBatchEntropy(bv, assayName = "raw", batchVar = "batch")
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' plotBatchEntropy(bv, assayName = "raw_combat", batchVar = "batch")
#'
#' @export

plotBatchEntropy <- function(
    object,
    assayName = "raw",
    batchVar = "batch",
    k = 10
) {
    ## -----------------------------
    ## Defensive checks
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    if (!assayName %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay not found: ", assayName)
    }

    sampleData <- SummarizedExperiment::colData(object)

    if (!batchVar %in% colnames(sampleData)) {
        stop("batchVar not found in colData: ", batchVar)
    }

    batch <- as.vector(sampleData[[batchVar]])

    if (any(is.na(batch))) {
        stop("batchVar contains NA values")
    }

    ## -----------------------------
    ## Data extraction
    ## -----------------------------
    mat <- SummarizedExperiment::assay(object, assayName)
    mat <- as.matrix(mat)

    nSamples <- ncol(mat)

    if (k >= nSamples) {
        stop("k must be less than the number of samples (", nSamples, ")")
    }

    ## -----------------------------
    ## PCA (robust PC selection)
    ## -----------------------------
    pca <- stats::prcomp(t(mat), scale. = TRUE)

    nPCs <- min(10, ncol(pca$x))
    pcs <- pca$x[, seq_len(nPCs), drop = FALSE]

    ## -----------------------------
    ## kNN computation
    ## -----------------------------
    nn <- FNN::get.knn(pcs, k = k)$nn.index

    ## -----------------------------
    ## Entropy calculation
    ## -----------------------------
    entropy <- apply(nn, 1, function(idx) {
        b <- batch[idx]

        p <- table(b) / length(b)

        ## numerical safety (avoid log(0))
        p <- p[p > 0]

        -sum(p * log(p))
    })

    ## -----------------------------
    ## Plot
    ## -----------------------------
    df <- data.frame(
        entropy = entropy,
        batch = batch
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$entropy)) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(paste("Batch mixing entropy:", assayName)) +
        ggplot2::xlab("Neighbourhood entropy")
}

#'
#' Plot PCA correction trajectories
#'
#' @param object BatchVariaData object
#' @param assayBefore baseline assay
#' @param assayAfter corrected assay
#' @param colourBy variable for colouring
#' @param reference Assay whose PCA basis both assays are projected onto.
#'   Defaults to \code{assayBefore}.
#'
#' @return ggplot object
#'
#' @details
#' Both assays are projected onto a basis fitted once on \code{reference},
#' so the arrows measure movement relative to fixed directions and the
#' 'before' positions are identical in every trajectory drawn from the same
#' object. Fitting a PCA jointly over the pair being plotted would instead
#' make those positions depend on which correction they were compared with,
#' so two trajectories from one object could not be overlaid.
#'
#' @seealso \code{\link{comparePCA}} for the shared-basis convention and
#'   its trade-off.
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- profileVariance(bv, ~batch, assayName = "raw")
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' bv <- profileVariance(bv, ~batch, assayName = "raw_combat")
#' plotPCATrajectory(bv, assayBefore = "raw", assayAfter = "raw_combat")
#' @export

plotPCATrajectory <- function(
    object,
    assayBefore = "raw",
    assayAfter = "raw_combat",
    colourBy = "batch",
    reference = NULL
) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (!is(object, "BatchVariaData")) {
        stop("object must be a BatchVariaData instance")
    }

    missingAssays <- setdiff(
        c(assayBefore, assayAfter),
        SummarizedExperiment::assayNames(object)
    )
    if (length(missingAssays) > 0) {
        stop("Assays not found: ", paste(missingAssays, collapse = ", "))
    }

    if (!colourBy %in% colnames(SummarizedExperiment::colData(object))) {
        stop("colourBy not found in colData: ", colourBy)
    }

    if (ncol(object) < 2) {
        stop("At least two samples are required for PCA")
    }

    ## Both assays are projected onto a basis fitted once on the reference,
    ## so the 'before' positions are the same in every trajectory drawn from
    ## this object. Fitting jointly over the chosen pair instead would make
    ## those positions depend on which correction they were compared with,
    ## and two trajectories from one object could not be overlaid.
    if (is.null(reference)) {
        reference <- assayBefore
    }

    if (!reference %in% SummarizedExperiment::assayNames(object)) {
        stop("Reference assay not found: ", reference)
    }

    basis <- .referenceBasis(object, reference)

    before_coords <- .projectOntoBasis(object, assayBefore, basis)[
        , seq_len(2), drop = FALSE
    ]
    after_coords <- .projectOntoBasis(object, assayAfter, basis)[
        , seq_len(2), drop = FALSE
    ]

    df <- data.frame(
        x1 = before_coords[, 1],
        y1 = before_coords[, 2],
        x2 = after_coords[, 1],
        y2 = after_coords[, 2],
        colour = SummarizedExperiment::colData(object)[[colourBy]]
    )

    ggplot2::ggplot(df) +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = .data$x1, y = .data$y1,
                xend = .data$x2, yend = .data$y2,
                colour = .data$colour
            ),
            arrow = grid::arrow(length = grid::unit(0.15, "cm"))
        ) +
        ggplot2::geom_point(
            ggplot2::aes(x = .data$x1, y = .data$y1),
            colour = "grey40"
        ) +
        ggplot2::geom_point(
            ggplot2::aes(x = .data$x2, y = .data$y2),
            colour = "black"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(
            paste(
                "PCA correction trajectory:",
                assayBefore, "->", assayAfter
            )
        )
}
