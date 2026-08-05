#'
#' BatchVariaData container
#'
#' An extension of SummarizedExperiment for batch-effect exploration,
#' correction, and variance decomposition analysis.
#'
#' Expression matrices (raw and corrected) are stored as assays.
#' Analysis provenance, including correction parameters and variance
#' decomposition history, is stored in metadata().
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @return A \code{BatchVariaData} object, an S4 class extending
#'   \code{SummarizedExperiment}.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100)
#' class(bv)
#' SummarizedExperiment::assayNames(bv)
#'
#' @export
setClass("BatchVariaData", contains = "SummarizedExperiment")
