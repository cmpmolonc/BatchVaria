#' BatchVaria: variance-aware evaluation of batch correction
#'
#' Batch correction is usually judged by whether a batch effect got
#' smaller. BatchVaria judges it by where the variance went. It
#' decomposes variance before and after a correction, records each
#' decomposition alongside the total variance it is relative to and the
#' share the design cannot attribute to any single term, and reports what
#' the correction did to the structure of the data as well as to the
#' batch term.
#'
#' @section Every function takes a ' \code{SummarizedExperiment}: expression matrices as assays, sample
#' annotation in \code{colData}, and the two analysis ledgers in
#' \code{metadata()}. The ledgers are optional and are written by
#' \code{\link{runCorrection}} and \code{\link{profileVariance}}, so an
#' ordinary \code{SummarizedExperiment} is a valid starting point. What
#' the package requires of an object, and checks on entry to every public
#' function, is set out in \link{BatchVaria-requirements}.
#'
#' @section Entry points:
#' \describe{
#'   \item{Correcting}{\code{\link{runCorrection}} applies a registered
#'     method and records what it did.
#'     \code{\link{availableCorrectionMethods}} lists them;
#'     \code{\link{registerCorrectionMethod}} adds one.}
#'   \item{Profiling}{\code{\link{profileVariance}} decomposes variance
#'     with one or more engines.
#'     \code{\link{availableVarianceMethods}} lists them;
#'     \code{\link{registerVarianceEngine}} adds one.}
#'   \item{Reading}{\code{\link{varianceTable}},
#'     \code{\link{varianceChange}} and \code{\link{varianceDelta}}
#'     summarise one decomposition. \code{\link{assayVariance}} supplies
#'     the total variance a fraction is relative to, without which
#'     fractions cannot be compared across assays.}
#'   \item{Structure}{\code{\link{basisProjection}} and
#'     \code{\link{basisRetention}} project every assay onto a single
#'     basis fitted on a reference assay, so \code{PC1} means the same
#'     direction in each.}
#'   \item{Evaluating}{\code{\link{evaluateCorrections}} runs the
#'     diagnostics in one call. \code{\link{provenance}} reports what has
#'     been done to an object.}
#'   \item{Plotting}{\code{\link{plotVarianceComposition}},
#'     \code{\link{plotVarianceDelta}}, \code{\link{plotVarianceRadar}},
#'     \code{\link{plotBasisProjection}},
#'     \code{\link{plotPCATrajectory}},
#'     \code{\link{plotSampleDistance}},
#'     \code{\link{plotBatchEntropy}} and
#'     \code{\link{plotCorrelationChange}}.}
#' }
#'
#' @section Why fractions need a denominator:
#' Variance fractions are compositional: they sum to one within an assay.
#' Removing variance associated with one term therefore inflates the
#' fraction attributed to every other, so a rising fraction cannot on its
#' own be read as preserved signal. Multiplying a fraction by its assay's
#' total variance gives an absolute quantity that is immune to this and
#' is the appropriate basis for comparing assays.
#'
#' @section Citation:
#' Moir N, Sherwood K, Simpson TI (2026). BatchVaria: a variance-aware
#' framework for evaluating batch correction in high-dimensional omics
#' data. \emph{bioRxiv}.
#' \doi{10.64898/2026.05.07.721996}
#'
#' @seealso \link{BatchVaria-requirements} for what an object must carry,
#'   \code{\link{exampleBatchVaria}} for a dataset to try it on, and
#'   \code{vignette("BatchVaria-quickstart", package = "BatchVaria")}.
#'
#' @keywords internal
"_PACKAGE"
