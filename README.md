
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lodi

## Overview

`lodi` is a package that implements censored likelihood multiple
imputation (CLMI) for single pollutant models with exposure biomarkers
below their respective detection limits. Additionally, implementations
for standard methods such as single imputation with a constant and
complete-case analysis are provided, although we do not recommend these
methods for datasets with a relatively high percent below their
respective detection limits (say \>25%).

You can learn more about how to use CLMI by working through the
example provided in `browseVignettes("lodi")`.

## Installation
`lodi` requires `rlang >= 0.3.0` to be installed, so you may want to update `rlang` before installing. 
```r
install.packages("lodi")
```
## Development version

To get a bug fix, or use a feature from the development version, you can
install lodi from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("umich-cphds/lodi", build_opts = c())
```

## Getting help

If you encounter a clear bug, please file a minimal reproducible example
on [github](https://github.com/umich-cphds/lodi/issues). For questions,
please email Jonathan Boss at <bossjona@umich.edu>.

## References

Boss J, Mukherjee B, Ferguson KK, et al. Estimating outcome-exposure
associations when exposure biomarker detection limits vary across
batches. *Epidemiology*. 2019. Epub ahead of print.
[10.1097/EDE.0000000000001052](https://doi.org/10.1097/EDE.0000000000001052)
