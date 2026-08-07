## Registry of variance profiling engines.
##
## Held in an environment rather than a switch() so that the set of engines
## is data rather than code: another package, or a script, can add one
## without editing BatchVaria. The built-ins are registered in .onLoad().
.varianceEngines <- new.env(parent = emptyenv())

## Arguments profileVariance() supplies to every engine.
.varianceEngineArgs <- c("assayMatrix", "formula", "sampleData")

#' Register a variance profiling engine
#'
#' Adds a variance decomposition engine to the set that
#' \code{\link{profileVariance}} can dispatch to, so that the engines are
#' extensible rather than fixed.
#'
#' @param name Character. Name callers will pass as \code{methods}.
#' @param engine A function implementing the engine contract below.
#' @param overwrite Logical. Replace an engine of the same name.
#'
#' @return The registered name, invisibly.
#'
#' @details
#' An engine is a function with the signature
#'
#' \preformatted{function(assayMatrix, formula, sampleData, ...)}
#'
#' where \code{assayMatrix} is features by samples with zero-variance
#' features already excluded, \code{formula} is the model, and
#' \code{sampleData} is the \code{colData} of the object. The \code{...} is
#' required: \code{profileVariance()} forwards engine-specific arguments
#' through it, and an engine must tolerate arguments meant for another.
#'
#' The return value must be a data.frame with one row per model term and
#' the columns
#' \itemize{
#'   \item \code{source} -- the engine's own name
#'   \item \code{term} -- the model term the variance is attributed to
#'   \item \code{variance_fraction} -- that term's share of total variance
#'   \item \code{n_features} -- features the decomposition was computed over
#' }
#'
#' Extra columns are allowed and are carried through to
#' \code{\link{varianceResults}}, where they are \code{NA} for engines that
#' do not supply them; the built-in \code{anova} engine uses this to record
#' its sums-of-squares type and weighting. \code{variance_fraction} must be
#' non-negative except for the reserved term \code{"shared"}, which holds
#' variance an unbalanced design cannot attribute to any single term and
#' may be negative under suppression.
#'
#' Two term names are reserved, because they mean the same thing whatever
#' the engine: \code{"residual"} for variance no model term accounts for,
#' and \code{"shared"} as above. An engine whose underlying implementation
#' spells these differently must normalise before returning -- the
#' built-in \code{variancePartition} engine renames \code{"Residuals"} on
#' the way out. Harmonising the columns but not the vocabulary would
#' produce results that look comparable and cannot be joined. Every other
#' term an engine may name as its model does.
#'
#' Use \code{\link{newVarianceSummary}} to build a conforming result rather
#' than assembling the data.frame by hand.
#'
#' @seealso \code{\link{availableVarianceMethods}},
#'   \code{\link{unregisterVarianceEngine}},
#'   \code{\link{newVarianceSummary}}
#'
#' @examples
#' ## an engine that attributes everything to the first model term
#' flatEngine <- function(assayMatrix, formula, sampleData, ...) {
#'     newVarianceSummary(
#'         source = "flat",
#'         term = c("everything", "residual"),
#'         varianceFraction = c(0.4, 0.6),
#'         nFeatures = nrow(assayMatrix)
#'     )
#' }
#'
#' registerVarianceEngine("flat", flatEngine)
#' availableVarianceMethods()
#'
#' set.seed(1)
#' bv <- exampleBatchVaria(nGenes = 50)
#' bv <- profileVariance(bv, ~batch, assays = "raw", methods = "flat")
#' varianceResults(bv, method = "flat")
#'
#' unregisterVarianceEngine("flat")
#'
#' @export
registerVarianceEngine <- function(name, engine, overwrite = FALSE) {
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !nzchar(name)) {
        stop("'name' must be a single non-empty string")
    }

    if (!is.function(engine)) {
        stop("'engine' must be a function")
    }

    formalNames <- names(formals(engine))

    if (!"..." %in% formalNames) {
        stop(
            "A variance engine must accept '...': profileVariance() ",
            "forwards engine-specific arguments to every engine, so each ",
            "must tolerate arguments intended for another"
        )
    }

    missingArgs <- setdiff(.varianceEngineArgs, formalNames)
    if (length(missingArgs) > 0) {
        stop(
            "A variance engine must accept ",
            paste(.varianceEngineArgs, collapse = ", "),
            "; '", name, "' is missing ",
            paste(missingArgs, collapse = ", ")
        )
    }

    if (!overwrite && exists(name, envir = .varianceEngines, inherits = FALSE)) {
        stop(
            "A variance engine named '", name, "' is already registered. ",
            "Pass overwrite = TRUE to replace it"
        )
    }

    assign(name, engine, envir = .varianceEngines)

    invisible(name)
}

#' Remove a registered variance profiling engine
#'
#' @param name Character. Name of the engine to remove.
#'
#' @return The removed name, invisibly.
#'
#' @seealso \code{\link{registerVarianceEngine}}
#'
#' @examples
#' dummy <- function(assayMatrix, formula, sampleData, ...) {
#'     newVarianceSummary(
#'         source = "dummy", term = "all",
#'         varianceFraction = 1, nFeatures = nrow(assayMatrix)
#'     )
#' }
#' registerVarianceEngine("dummy", dummy)
#' unregisterVarianceEngine("dummy")
#' availableVarianceMethods()
#'
#' @export
unregisterVarianceEngine <- function(name) {
    if (!exists(name, envir = .varianceEngines, inherits = FALSE)) {
        stop(
            "No variance engine named '", name, "' is registered. ",
            "Registered engines: ",
            paste(availableVarianceMethods(), collapse = ", ")
        )
    }

    rm(list = name, envir = .varianceEngines)

    invisible(name)
}

.getVarianceEngine <- function(method) {
    if (!exists(method, envir = .varianceEngines, inherits = FALSE)) {
        stop(
            "Unknown variance method: ", method,
            ". Available methods are: ",
            paste(availableVarianceMethods(), collapse = ", ")
        )
    }

    get(method, envir = .varianceEngines, inherits = FALSE)
}
