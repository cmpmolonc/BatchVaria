## A single, shared basis for every cross-assay PCA view.
##
## Fitting a PCA independently per assay makes "PC1" a different direction
## in every column: the label is stable, the axis it names is not. Fitting
## jointly over a chosen pair of assays is better, but the basis then
## depends on which pair was chosen, so the same samples land in different
## positions in different plots and two views of one object cannot be
## overlaid.
##
## Projecting every assay onto a basis fitted once on a designated
## reference makes the axes mean the same thing everywhere, and makes a
## figure reproducible from the object alone rather than from the object
## plus the set of assays that happened to be passed in. The cost is that
## the reference is privileged, which is explicit rather than incidental.
.referenceBasis <- function(object, reference) {
    mat <- as.matrix(SummarizedExperiment::assay(object, reference))

    featureVar <- .rowVars(mat)
    flat <- !is.finite(featureVar) | featureVar <= 0

    if (sum(!flat) < 2L) {
        stop(
            "Assay '", reference, "' has fewer than two features with ",
            "non-zero variance; cannot fit a PCA basis"
        )
    }

    ## prcomp(scale. = TRUE) cannot rescale a constant feature. Dropping
    ## such features would remove them from the coordinate system entirely,
    ## so variance a correction later introduces there would be invisible
    ## rather than counted as lying off the basis - and a correction
    ## putting signal into a previously flat feature is exactly the
    ## pathology worth catching. Instead they are kept on their original
    ## scale: they carry no variance in the reference, so no component
    ## loads on them, but they remain part of the space being measured.
    scaleVec <- sqrt(featureVar)
    scaleVec[flat] <- 1

    centre <- rowMeans(mat)

    scaled <- scale(t(mat), center = centre, scale = scaleVec)
    fit <- stats::prcomp(scaled, center = FALSE, scale. = FALSE)

    list(
        reference = reference,
        features = rownames(mat),
        flat = rownames(mat)[flat],
        center = centre,
        scale = scaleVec,
        rotation = fit$rotation
    )
}

## Scores of one assay in the reference basis, using the reference's own
## centring and scaling so that the coordinates are directly comparable.
.projectOntoBasis <- function(object, assayName, basis) {
    mat <- as.matrix(SummarizedExperiment::assay(object, assayName))

    missingFeatures <- setdiff(basis$features, rownames(mat))
    if (length(missingFeatures) > 0) {
        stop(
            "Assay '", assayName, "' is missing ", length(missingFeatures),
            " features present in the reference assay '", basis$reference, "'"
        )
    }

    scaled <- scale(
        t(mat[basis$features, , drop = FALSE]),
        center = basis$center,
        scale = basis$scale
    )

    scaled %*% basis$rotation
}

## Variance carried along each reference axis, as a fraction of the
## assay's own total variance in the reference's scaled space.
##
## For the reference assay these fractions sum to one. For any other assay
## they sum to less than one whenever correction has moved variance into
## directions the reference basis does not span - the shortfall is
## informative, and is why these are not a partition.
.basisVarianceFractions <- function(object, assayName, basis) {
    scores <- .projectOntoBasis(object, assayName, basis)

    mat <- as.matrix(SummarizedExperiment::assay(object, assayName))
    scaled <- scale(
        t(mat[basis$features, , drop = FALSE]),
        center = basis$center,
        scale = basis$scale
    )

    total <- sum(.rowVars(t(scaled)))

    if (!is.finite(total) || total <= 0) {
        stop("Assay '", assayName, "' has no variance to project")
    }

    .rowVars(t(scores)) / total
}

#' Sample coordinates on a reference assay's principal axes
#'
#' Projects each assay onto a basis fitted once on a reference assay and
#' returns the sample coordinates.
#'
#' @param object A \code{SummarizedExperiment} meeting BatchVaria's
#'   requirements. See \link{BatchVaria-requirements}.
#' @param assays Character vector of assays to project (default: all).
#' @param nPCs Number of components to return (default: 2).
#' @param reference Assay supplying the basis. Inferred from the correction
#'   lineage when \code{NULL}, falling back to the first assay.
#'
#' @return A data.frame with one row per assay/sample combination and
#'   columns \code{assay}, \code{sample}, \code{reference} and one column
#'   per requested component, named \code{PC1}, \code{PC2} and so on.
#'
#' @details
#' Every assay is projected onto one basis, using the reference's own
#' centring and scaling, so coordinates are directly comparable across
#' assays: a sample that has not moved has the same position in every
#' assay. Fitting a PCA separately per assay would make \code{PC1} a
#' different direction in each, so the same sample would land in different
#' positions depending only on which assay it was drawn from. Fitting over
#' a chosen pair of assays would make the position depend on what the assay
#' was compared against.
#'
#' Because the axes are the reference's, they are principal components of
#' the reference and not of the assay being projected. Coordinates for a
#' corrected assay describe where its samples sit in the original
#' structure, which is the comparison of interest; they are not that
#' assay's own principal components and the components are not ordered by
#' its variance. \code{\link{comparePCA}} reports how much of each assay's
#' variance each axis carries, and \code{\link{basisRetention}} sums that
#' over the whole basis.
#'
#' This is the function \code{\link{plotBasisProjection}} draws, exported
#' so the coordinates are available without going through a plot.
#'
#' @seealso \code{\link{basisRetention}} for the variance summary,
#'   \code{\link{comparePCA}} for the per-component breakdown,
#'   \code{\link{plotBasisProjection}} to plot the coordinates.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' head(basisProjection(bv, assays = c("raw", "raw_combat")))
#'
#' @export
basisProjection <- function(object, assays = NULL, nPCs = 2,
                            reference = NULL) {
    .check_se(object)

    available <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available
    }

    missingAssays <- setdiff(assays, available)
    if (length(missingAssays) > 0) {
        stop("Assays not found: ", paste(missingAssays, collapse = ", "))
    }

    if (!is.numeric(nPCs) || length(nPCs) != 1L || is.na(nPCs) || nPCs < 1) {
        stop("'nPCs' must be a single number >= 1")
    }

    if (ncol(object) < 2) {
        stop("At least two samples are required for a projection")
    }

    reference <- .resolveBaseline(object, assays, reference)
    if (is.null(reference)) {
        reference <- assays[1]
    }

    if (!reference %in% available) {
        stop("Reference assay not found: ", reference)
    }

    basis <- .referenceBasis(object, reference)

    sampleNames <- colnames(object)
    if (is.null(sampleNames)) {
        sampleNames <- as.character(seq_len(ncol(object)))
    }

    rows <- lapply(assays, function(a) {
        scores <- .projectOntoBasis(object, a, basis)

        keep <- min(nPCs, ncol(scores))
        coords <- as.data.frame(scores[, seq_len(keep), drop = FALSE])
        colnames(coords) <- paste0("PC", seq_len(keep))

        out <- data.frame(
            assay = a,
            sample = sampleNames,
            reference = reference,
            stringsAsFactors = FALSE
        )

        cbind(out, coords)
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#' Variance retained on a reference assay's principal axes
#'
#' Reports how much of each assay's variance lies along the principal
#' directions of a reference assay, and how much has moved off them.
#'
#' @param object A \code{SummarizedExperiment} meeting BatchVaria's
#'   requirements. See \link{BatchVaria-requirements}.
#' @param assays Character vector of assays to measure (default: all).
#' @param reference Assay supplying the basis. Inferred from the correction
#'   lineage when \code{NULL}, falling back to the first assay.
#'
#' @return A data.frame with one row per assay and columns \code{assay},
#'   \code{reference}, \code{retention}, \code{off_basis} and
#'   \code{n_features}.
#'
#' @details
#' \code{retention} is the fraction of an assay's total variance that the
#' reference's principal axes account for; \code{off_basis} is the
#' remainder. For the reference itself retention is one by construction.
#'
#' The interpretation is a single distinction that no independently fitted
#' PCA can draw: a correction that redistributes variance *within* the
#' structure already present keeps retention near one, while a correction
#' that introduces new structure drives it down. On the example data ComBat
#' retains about 0.99 while an assay distorted with added noise retains
#' about 0.41 - both look like changed variance in a per-assay
#' decomposition, but only the second has moved variance into directions
#' the original data did not occupy.
#'
#' Features with no variance in the reference are kept in the basis on
#' their original scale rather than dropped, so variance a correction
#' introduces in a previously flat feature is counted as off-basis instead
#' of disappearing from the measurement.
#'
#' @seealso \code{\link{comparePCA}} for the per-component breakdown.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' bv <- runCorrection(bv, method = "combat", batch = "batch")
#' basisRetention(bv, assays = c("raw", "raw_combat", "raw_noise"))
#'
#' @export
basisRetention <- function(object, assays = NULL, reference = NULL) {
    .check_se(object)

    available <- SummarizedExperiment::assayNames(object)

    if (is.null(assays)) {
        assays <- available
    }

    missingAssays <- setdiff(assays, available)
    if (length(missingAssays) > 0) {
        stop("Assays not found: ", paste(missingAssays, collapse = ", "))
    }

    reference <- .resolveBaseline(object, assays, reference)
    if (is.null(reference)) {
        reference <- assays[1]
    }

    if (!reference %in% available) {
        stop("Reference assay not found: ", reference)
    }

    basis <- .referenceBasis(object, reference)

    tol <- sqrt(.Machine$double.eps)

    rows <- lapply(assays, function(a) {
        retention <- sum(.basisVarianceFractions(object, a, basis))

        ## an assay fully inside the basis returns 1 give or take rounding;
        ## report that as an exact 1 rather than as a negative off_basis
        if (abs(retention - 1) < tol) {
            retention <- 1
        }

        data.frame(
            assay = a,
            reference = reference,
            retention = retention,
            off_basis = 1 - retention,
            n_features = length(basis$features),
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}
