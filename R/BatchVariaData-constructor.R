#' Construct a BatchVariaData object
#'
#' Creates a BatchVariaData container built on a SummarizedExperiment.
#'
#' @param se A SummarizedExperiment containing expression data.
#'
#' @return A BatchVariaData object
#'
#' @seealso \code{\link{BatchVariaData-class}}
#'
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'     assays = list(raw = matrix(rnorm(200),
#'         nrow = 20,
#'         dimnames = list(
#'             paste0("g", seq_len(20)),
#'             paste0("s", seq_len(10))
#'         )
#'     )),
#'     colData = DataFrame(
#'         batch = factor(rep(c("A", "B"), 5)),
#'         row.names = paste0("s", seq_len(10))
#'     )
#' )
#' bv <- BatchVariaData(se)
#' bv
#' @export
BatchVariaData <- function(se = SummarizedExperiment::SummarizedExperiment()) {
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object")
    }

    out <- new("BatchVariaData") # empty, valid subclass
    for (nm in methods::slotNames("SummarizedExperiment")) {
        methods::slot(out, nm) <- methods::slot(se, nm)
    } # copy assays/colData/rowData/metadata/NAMES

    validObject(out)
    out
}
