## Progress output from registered engines.
##
## The split this enforces: 'verbose' governs what the pipeline says about
## its progress, not what it says about problems it found. Everything
## BatchVaria derives itself - no-op detection, negative 'shared',
## fingerprint staleness, dropped features, engine failure - is raised
## with warning() and is untouched here.
##
## Applied at the dispatch site rather than inside the built-in engines, so
## engines added through registerCorrectionMethod() and
## registerVarianceEngine() inherit the policy without implementing it.
.withEngineOutput <- function(verbose, expr) {
    if (isTRUE(verbose)) {
        return(expr)
    }

    suppressMessages(expr)
}
