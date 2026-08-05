CHANGES IN VERSION 0.99.0
-------------------------

NEW FEATURES

    o Initial Bioconductor submission.

    o BatchVariaData container extending SummarizedExperiment, recording
      correction and variance-decomposition provenance in metadata().

    o run_correction() applies batch correction and appends each step to a
      correction ledger.

    o profile_variance() harmonises the PCA, ANOVA and variancePartition
      engines into a common tidy variance schema.

    o varianceTable(), variance_delta() and evaluate_corrections() compare
      variance composition, PCA structure and correlation preservation
      before and after correction.

    o Diagnostic plots for variance composition, variance change, batch
      entropy, sample distance and PCA trajectory.
