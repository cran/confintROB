# confintROB

Contains a set of R functions to compute percentile or BCa confidence intervals based on bootstrap schemes, for linear mixed models estimated with both classical and robust methods. 

The scripts allow for balanced and unbalanced data, using [`rlmer`](https://cran.r-project.org/package=robustlmm) ([Koller, 2016](https://doi.org/10.18637/jss.v075.i06)) for robust estimation or [`lmer`](https://cran.r-project.org/package=lme4) ([Bates et al., 2015](https://doi.org/10.18637/jss.v067.i01)) for classical estimation. 

Two bootstrap schemes are implemented: the wild bootstrap and the parametric bootstrap. 

Both are adapted for object of class `lmerMod` (with [`lmer`](https://cran.r-project.org/package=lme4)), `rlmerMod` (with [`rlmer`](https://cran.r-project.org/package=robustlmm)) or `varComprob` ([`varComprob`](https://cran.r-project.org/package=robustvarComp)) ([Agostinelli & Yohai, 2016](https://doi.org/10.1080/01621459.2015.1115358)). 

A comparison of performance with balanced datasets is available in [Mason, Cantoni and Ghisletta (2024)](https://psycnet.apa.org/record/2024-54080-001).  

### Contents

The main function is `confintROB`; it is implemented as an extension of the popular `confint` (`confint.merMod`) function of the [lme4 package](https://cran.r-project.org/package=lme4).



