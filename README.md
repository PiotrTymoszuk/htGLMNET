# htGLMNET
High Throughput Light Weight Regularized Regression Modeling for Molecular Data

## Summary
The package enables for fitting and basic evaluation of regularized regression models via the interface of the _glmnet_ package [^1]. 
The pre-processing, modeling, evalaution, and prediction tools were fine-tuned specifically for use with genomic, transcriptomic, proteomic and pharmaco-molecular data.

## Instalation

The package requires one non-CRAN package (_microViz_, https://github.com/PiotrTymoszuk/microViz), please make sure to install it first:

```
devtools::install_github('PiotrTymoszuk/microViz') ## dependency microViz

devtools::install_github('PiotrTymoszuk/htGLMNET')

```
## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/biggrExtra/blob/main/LICENSE).


## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Acknowledgements

Many thanks to the authors of the packages [_glmnet_](https://glmnet.stanford.edu/articles/glmnet.html), [_sva_](https://bioconductor.org/packages/sva), [_tidyverse_](https://www.tidyverse.org/), [_rlang_](https://rlang.r-lib.org/), [_furrr_](https://furrr.futureverse.org/), , [ggrepel](https://ggrepel.slowkow.com/), and [_caret_](https://topepo.github.io/caret/). 

Special thanks to Danielle Maese, Robert Gruener and Stephanie Huang, the authors of the _oncoPredict_ and _pRRophetic_ for the idea of using regularized linear regression for prediction of anti-cancer drug response.

## References

[^1]: Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. J Stat Softw (2010) 33:1–22. doi:10.18637/jss.v033.i01


