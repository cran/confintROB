# confintROB 1.2

* The parametric bootstrap (`boot.type = "parametric"`) for `lmer` and `rlmer`
  fits now generates both the random effects and the errors anew in every
  bootstrap sample, as described in Mason, Cantoni & Ghisletta (2021). In
  earlier versions the random effects were held fixed at their predicted
  values and only the errors were redrawn, which gave intervals that were too
  narrow, most notably for the fixed effects and the random effects' standard
  deviations. This affects `method = "boot"` and `method = "BCa"`. The default
  `boot.type = "wild"`, Wald intervals and `varComprob` fits are not affected.

* The wild bootstrap (`boot.type = "wild"`) and the jackknife for
  `method = "BCa"` of `lmer` and `rlmer` fits now use independent clusters
  when the model has several grouping factors. Earlier versions used the
  grouping factor with the most levels, which for nested designs such as
  `(1 | batch/cask)` is the innermost one (`cask:batch`): the wild bootstrap
  gave correlated casks independent weights and understated the variability
  between batches, and `method = "BCa"` failed and returned percentile
  intervals with a warning. The clusters are now the coarsest grouping factor
  that contains all the others (`batch`), the same rule as
  `simulate(type = "wild")` in robustlmm. Models with crossed grouping factors
  have no independent clusters and now stop with an error unless `clusterID`
  is given. `clusterID`, the name of a grouping factor or a vector with one
  value per observation, is no longer ignored for `lmer` and `rlmer` fits with
  `boot.type = "wild"` or `method = "BCa"`. Models with a single grouping
  factor give the same results as before (apart from the BCa jackknife fix
  below).

* The jackknife that estimates the acceleration of `method = "BCa"` for `lmer`
  and `rlmer` fits now leaves out one whole cluster at a time. Earlier
  versions split the observations by the cluster sizes instead of by cluster,
  so each jackknife sample left out a single observation. This changes all
  BCa intervals for `lmer` and `rlmer` fits; `varComprob` fits are not
  affected.
