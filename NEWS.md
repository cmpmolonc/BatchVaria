# BatchVaria 0.99.0

Initial Bioconductor submission.

## Container and provenance

* `BatchVariaData` extends `SummarizedExperiment`, holding the raw and
  corrected matrices as assays and every analysis decision in `metadata()`.

* `runCorrection()` applies a batch correction and appends an entry to the
  correction ledger recording the method, the input and output assays, the
  batch variable, any variables preserved during correction, and the total
  variance before and after. `preserve` names variables whose variation is
  to be kept; note that `limma::removeBatchEffect()` has an argument spelled
  `covariates` that means the opposite.

* A correction that returns an assay indistinguishable from its input warns
  and is recorded in the ledger as `no_op`. A method can succeed and change
  nothing -- `limma::removeBatchEffect()` does exactly this when batch is
  aliased with the preserved design, and reports it through a message rather
  than a trappable condition. Undetected, the result flows downstream as a
  correction and every comparison against the input assay becomes a
  comparison of an assay with itself. The flag is stored on the object, so
  the condition outlives the console.

* `summary()` shows the variables preserved by each correction and the
  `no_op` flag, so two runs of one method on one assay are distinguishable
  in the printed history. Downstream functions read that lineage to work out which assay a
  comparison should be made against, so no assay has to be named `raw`.

* `profileVariance()` records every decomposition in a variance ledger keyed
  on assay, method, formula and term. Formula keys are canonical, so
  `~ batch + group` and `~ group + batch` are recognised as one model while
  genuinely different models coexist without colliding.

* `varianceHistory()` and `varianceResults()` expose the ledger, the latter
  in tidy form with the canonical formula and each engine's own diagnostic
  columns.

## Variance profiling

* Two engines share a common schema: `anova`, which fits a linear model per
  feature, and `variancePartition`, which fits mixed models. Both attribute
  variance to model terms, so their results are directly comparable. On a
  balanced design they agree to three decimal places.

* The `anova` engine uses Type II sums of squares and emits one row per model
  term plus `residual` and `shared`. `shared` carries variance an unbalanced
  design cannot attribute to any single term: zero when the design is
  orthogonal, growing with confounding, and absorbing everything two terms
  explain when they are completely aliased. Feature-averaged and pooled
  weightings are both available, and the choice is recorded with the result.

* Features with zero or non-finite variance are excluded once per assay and
  the count reported. These are routine in unfiltered count data and would
  otherwise cause engine-specific failures.

* Each assay and method is attempted independently. One engine failing warns
  and is skipped rather than discarding the engines that succeeded; an error
  is raised only if every attempt fails.

* Engines build their own design matrices, so a formula written for one is
  never forced through another. `variancePartition` accepts random-effects
  notation such as `~ (1 | batch)`; `anova` models fixed effects only and
  declines such formulas by name.

## Reading variance correctly

* `assayVariance()` reports the total variance of each assay. Variance
  fractions are compositional, so removing variance associated with one term
  inflates the fraction attributed to every other. Multiplying a fraction by
  its assay's total variance gives an absolute quantity that is comparable
  across assays and immune to that effect.

* `varianceTable()`, `varianceChange()` and `varianceDelta()` summarise
  exactly one decomposition. When the ledger holds more than one method or
  formula they require `method` or `formula` to be named, listing what is
  available, rather than silently mixing results whose terms are not
  commensurable. The baseline is inferred from the correction lineage, and
  when none can be determined the composition table is still returned.

## Embeddings

* `comparePCA()`, `plotPCA()` and `plotPCATrajectory()` project every assay
  onto a single basis fitted on a reference assay. A per-assay fit would make
  `PC1` a different direction in each, and a fit over a chosen pair would move
  the same samples depending on what they were compared against.

* `basisRetention()` reports how much of an assay's variance lies on the
  reference's principal axes. A correction that redistributes variance within
  the existing structure retains close to all of it; one that introduces new
  structure does not.

* PCA is not among the variance engines. It is unsupervised, so its output is
  not attributable to a covariate, and its terms are latent axes rather than
  model terms.

## Extensibility

* `registerVarianceEngine()` and `unregisterVarianceEngine()` add and remove
  variance engines, so the set is extensible rather than fixed. The built-in
  engines register through the same entry point.

* `registerCorrectionMethod()` and `unregisterCorrectionMethod()` do the same
  for batch correction, behind a contract of the same shape. A method takes a
  matrix, the name of a single batch column, the sample data and an optional
  set of variables to preserve, and returns a matrix whose dimnames are
  identical to the input's. Sample order is part of the contract because
  every downstream comparison is keyed on it. Two methods ship: `combat`
  and `limma`.

* Correction methods build their own design matrices from `preserve`, so a
  design assembled for one is never forced through another.

* `newVarianceSummary()` builds and validates a conforming result. The engine
  contract is documented on `registerVarianceEngine()`.

## Evaluation and diagnostics

* `evaluateCorrections()` runs variance, PCA and correlation diagnostics in
  one call. `compareCorrelations()` measures the mean absolute change in
  sample-sample correlation structure.

* Plotting: `plotVarianceComposition()`, `plotVarianceDelta()`,
  `plotVarianceRadar()`, `plotBatchEntropy()`, `plotSampleDistance()`,
  `plotCorrelationChange()`, `plotPCA()` and `plotPCATrajectory()`.

## Example data

* `exampleBatchVaria()` generates a two-batch, two-group design. The
  `confounding` argument slides it from orthogonal to fully aliased, so the
  behaviour of `shared`, and the divergence between the two engines under
  confounding, can be demonstrated directly.
