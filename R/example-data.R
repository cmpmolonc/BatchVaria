#' Example BatchVaria dataset with batch and biological effects
#'
#' Generates a synthetic SummarizedExperiment and converts it to
#' a BatchVariaData object. Includes several pre-generated corrected
#' assays for vignette demonstrations.
#'
#' @examples
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 100, nSamples = 8)
#' bv
#' SummarizedExperiment::assayNames(bv)
#' @param nGenes Number of genes. Must be at least 10, so that each effect
#'   block (group, batch, both) contains at least one gene.
#' @param nSamples Number of samples. Must be a multiple of 4 and at least
#'   4, giving a balanced 2 x 2 design of two batches by two groups.
#'
#' @return BatchVariaData object
#'
#' @details
#' The expression matrix and the pre-generated \code{raw_noise} assay are
#' drawn with \code{\link[stats]{rnorm}}, so repeated calls return different
#' data. Seeding is deliberately left to the caller rather than done inside
#' the function: call \code{set.seed()} beforehand whenever you need
#' reproducible output.
#'
#' @seealso \code{\link{set.seed}}
#'
#' @export
exampleBatchVaria <- function(
    nGenes = 800,
    nSamples = 20
) {
    # -----------------------------
    # Validate the design
    # -----------------------------
    if (!is.numeric(nSamples) || length(nSamples) != 1 || is.na(nSamples)) {
        stop("'nSamples' must be a single number")
    }

    if (nSamples < 4 || nSamples %% 4 != 0) {
        nearest <- unique(c(
            max(4, 4 * floor(nSamples / 4)),
            max(4, 4 * ceiling(nSamples / 4))
        ))
        stop(
            "'nSamples' must be a multiple of 4 and at least 4: the example data ",
            "is a balanced 2 x 2 design (2 batches x 2 groups). Nearest valid ",
            "value(s) to ", nSamples, ": ", paste(nearest, collapse = " or ")
        )
    }

    if (!is.numeric(nGenes) || length(nGenes) != 1 || is.na(nGenes)) {
        stop("'nGenes' must be a single number")
    }

    if (nGenes < 10) {
        stop(
            "'nGenes' must be at least 10 so that each effect block ",
            "(group, batch, both) contains at least one gene"
        )
    }

    # -----------------------------
    # Sample metadata
    # -----------------------------
    batch <- rep(c("B1", "B2"), each = nSamples / 2)
    group <- rep(rep(c("Ctrl", "Trt"), each = nSamples / 4), 2)

    coldata <- S4Vectors::DataFrame(
        sample = paste0("S", seq_len(nSamples)),
        batch = batch,
        group = group
    )
    rownames(coldata) <- coldata$sample

    # -----------------------------
    # Base expression matrix
    # -----------------------------
    expr <- matrix(
        stats::rnorm(nGenes * nSamples, mean = 0, sd = 1),
        nrow = nGenes
    )

    # Gene subsets: first 20% respond to group, next 20% to batch,
    # next 10% to both. Integer counts keep indexing exact.
    n_group <- round(0.20 * nGenes)
    n_batch <- round(0.20 * nGenes)
    n_both <- round(0.10 * nGenes)

    g_group <- seq_len(n_group) # 1 .. n_group
    g_batch <- n_group + seq_len(n_batch) # next block
    g_both <- n_group + n_batch + seq_len(n_both) # next block

    # Biological effect
    expr[g_group, group == "Trt"] <-
        expr[g_group, group == "Trt"] + 2

    # Batch effect
    expr[g_batch, batch == "B2"] <-
        expr[g_batch, batch == "B2"] + 2

    # Both effects
    expr[g_both, group == "Trt"] <-
        expr[g_both, group == "Trt"] + 1.5

    expr[g_both, batch == "B2"] <-
        expr[g_both, batch == "B2"] + 1.5

    rownames(expr) <- paste0("Gene", seq_len(nGenes))
    colnames(expr) <- coldata$sample

    # -----------------------------
    # Create SummarizedExperiment
    # -----------------------------
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = expr),
        colData = coldata
    )

    # -----------------------------
    # Convert to BatchVariaData
    # -----------------------------
    bv <- BatchVariaData(se)

    # -----------------------------
    # Add pre-generated corrections
    # -----------------------------
    bv <- .addExampleCorrections(bv)

    bv
}

.addExampleCorrections <- function(bv) {
    raw <- SummarizedExperiment::assay(bv, "raw")

    # Over-correction
    raw_center <- t(scale(t(raw), center = TRUE, scale = FALSE))

    # Under-correction
    raw_scale <- raw * 0.9

    # Bad correction (noise)
    raw_noise <- raw + matrix(
        stats::rnorm(length(raw), sd = 1.5),
        nrow = nrow(raw)
    )

    SummarizedExperiment::assay(bv, "raw_center") <- raw_center
    SummarizedExperiment::assay(bv, "raw_scale") <- raw_scale
    SummarizedExperiment::assay(bv, "raw_noise") <- raw_noise

    bv
}
