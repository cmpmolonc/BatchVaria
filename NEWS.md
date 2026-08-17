# BatchVaria 0.99.0

Initial Bioconductor submission.

## Container and provenance

* BatchVaria works on a plain `SummarizedExperiment`: raw and corrected
  matrices as assays, every analysis decision in `metadata()`. There is no
  container class. Earlier drafts, including the preprint, describe a
  `BatchVariaData` class extending `SummarizedExperiment`. It has been
  removed, with its constructor and its `show` and `summary` methods.
  `exampleBatchVaria()` returns a `SummarizedExperiment`, and
  `provenance()` replaces `summary()`.

* What the original class's validity method checked is now checked at the point of
  use, by an internal `.check_se()` that every public function calls:
  ledger entries well formed, and the assays and `colData` columns they
  name present. This is strictly stronger than the validity method was.
  S4 validity runs at construction and on an explicit `validObject()` and
  at no other point - `[`, `assay<-`, `colData<-` and `metadata<-` do not
  invoke it - so it could only ever see the object as handed to the
  constructor, never as it is when a function is asked to use it. Both
  ledgers remain optional, so an ordinary `SummarizedExperiment` that
  nothing has been run on is a valid input. The requirements are
  documented at `?"BatchVaria-requirements"`.

* Referential integrity and currency are different questions, and the
  second is answered by a fingerprint recorded with each ledger entry and
  checked when the entry is read. Subsetting preserves every assay name
  and `colData` column, so an entry valid before `bv[, 1:6]` is still
  valid after it - and it previously left `varianceTable()` returning
  percentages that described samples the object no longer held.
  `varianceTable()` and
  `varianceResults()` now warn, naming the assays and how they changed.
  They still return the result: a stale decomposition is evidence of what
  was computed, it has only stopped describing the current object. The
  fingerprint covers content as well as dimensions, so an assay replaced
  in place by a matrix of the same shape is detected too.

* `runCorrection()` applies a batch correction and appends an entry to the
  correction ledger recording the method, the input and output assays, the
  batch variable, any variables preserved during correction, and the total
  variance before and after. `preserve` names variables whose variation is
  to be kept; note that `limma::removeBatchEffect()` has an argument named 
  `covariates` that means the opposite.

* A correction that returns an assay indistinguishable from its input warns
  and is recorded in the ledger as `no_op`. A method can succeed and change
  nothing - `limma::removeBatchEffect()` does exactly this when batch is
  aliased with the preserved design, and reports it through a message rather
  than a trappable condition. Undetected, the result flows downstream as a
  correction and every comparison against the input assay becomes a
  comparison of an assay with itself. The flag is stored on the object, so
  the condition outlives the console.

* `provenance()` reports what has been done to an object: its assays, the
  corrections that produced each one, and the variance decompositions
  recorded against them. It shows the variables preserved by each
  correction and the `no_op` flag, so two runs of one method on one assay
  are distinguishable in the printed history. Downstream functions read
  that lineage to work out which assay a comparison should be made
  against, so no assay has to be named `raw`.

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

* `profileVariance()` runs `anova` alone unless told otherwise. Profiling
  with every registered engine would leave the ledger holding several
  decompositions after one call, so the shortest path through the package
  ended at `varianceTable()` asking which was meant. It would also make
  `registerVarianceEngine()` change the behaviour of calls that never named
  the new engine. Comparing engines is a deliberate act: name them.

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

* `basisProjection()`, `comparePCA()`, `plotBasisProjection()` and
  `plotPCATrajectory()` project every assay onto a single basis fitted on a
  reference assay. A per-assay fit would make `PC1` a different direction in
  each, and a fit over a chosen pair would move the same samples depending on
  what they were compared against.

* `basisProjection()` returns the sample coordinates as a data.frame.
  Previously the projection was reachable only through a plot, so the one
  thing in this layer with no equivalent elsewhere - a shared reference basis
  rather than a per-assay fit - could not be obtained without unpacking a
  ggplot. `basisRetention()` published a summary computed from coordinates
  nobody could get at. `plotBasisProjection()` now draws what
  `basisProjection()` returns rather than recomputing it.

* `plotPCA()` is renamed `plotBasisProjection()` and is a plain function
  rather than a method on `BiocGenerics::plotPCA()`. The name described the
  technique and collided with other packages' `plotPCA()`; what distinguishes
  this one is the basis, and after projection the axes are the reference's
  principal components rather than those of the assay being drawn.

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

* The engine contract reserves two term names, `residual` and `shared`,
  because they mean the same thing whatever the engine. An engine whose
  implementation spells them differently normalises before returning: the
  built-in `variancePartition` engine renames `Residuals` on the way out.
  Harmonising the columns but not the vocabulary would produce results that
  look comparable and cannot be joined. Every other term an engine names as
  its model does.

* `newVarianceSummary()` builds and validates a conforming result. The engine
  contract is documented on `registerVarianceEngine()`.

* `runCorrection()` and `profileVariance()` take `verbose`, default `FALSE`,
  which governs progress output from the engine they call and nothing else.
  Everything the package derives itself - no-op detection, negative `shared`,
  fingerprint staleness, dropped features, engine failure - is a warning and
  is raised whatever `verbose` is set to. The separation is structural rather
  than a convention to maintain: progress is suppressed by silencing message
  conditions, and no diagnostic in the package is a message.

  Suppression is applied where the engine is invoked, so methods and engines
  added through the registries inherit it without implementing it. Messages
  only, never `capture.output()`: `sva::ComBat()` reports progress through
  message conditions but reports one finding about the data - features with
  uniform expression within a batch - on stdout, and capturing stdout would
  delete it. That notice is passed through as sva's own rather than re-emitted
  as a BatchVaria warning, since re-broadcasting a dependency's finding means
  owning its wording across versions.

## Evaluation and diagnostics

* `evaluateCorrections()` runs variance, PCA and correlation diagnostics in
  one call. `compareCorrelations()` measures the mean absolute change in
  sample-sample correlation structure.

* Plotting: `plotVarianceComposition()`, `plotVarianceDelta()`,
  `plotVarianceRadar()`, `plotBatchEntropy()`, `plotSampleDistance()`,
  `plotCorrelationChange()`, `plotBasisProjection()` and
  `plotPCATrajectory()`.

## Example data

* `exampleBatchVaria()` generates a two-batch, two-group design and returns a
  `SummarizedExperiment`. The `confounding` argument slides it from orthogonal
  to fully aliased, so the behaviour of `shared`, and the divergence between
  the two engines under confounding, can be demonstrated directly.

## Documentation

* `?BatchVaria` gives an overview: what the package is for, the entry points
  grouped by what they answer, why a variance fraction needs its denominator,
  and the citation. It links `?"BatchVaria-requirements"`, which is where a
  reader arriving at "bring your own SummarizedExperiment" finds out what
  that object has to carry.

* The vignette gains a diagnostics section covering the plotting layer -
  composition, radar, basis projection, PCA trajectory, sample distance,
  correlation change and batch entropy - each with what it answers and how to
  read it. It sits after the numerical sections rather than interleaved,
  because every claim the vignette makes is a number, and the plots draw those
  numbers rather than standing in for them.

* The vignette also shows `provenance()`, and states plainly that the object
  is an ordinary `SummarizedExperiment` throughout.
