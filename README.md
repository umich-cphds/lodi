
<!-- README.md is generated from README.Rmd. Please edit that file -->

`lodi` is a R package that implements censored likelihood multiple
imputation (CLMI) for single pollutant models with exposure biomarkers
below their respective detection limits. `lodi` also contains implementations
for standard methods such as single imputation with a constant and
complete-case analysis, although those methods are primarily designed for comparison with `clmi`.
## Installation
`lodi` requires `rlang >= 0.3.0` to be installed, so you may want to install /update `rlang` before installing `lodi`.

The package can be installed from CRAN 
```r
install.packages("lodi")
```
Or from Github
``` r
# install.packages("devtools")
devtools::install_github("umich-cphds/lodi", build_opts = c())
```
The github version may contain bugfixes not yet present on CRAN, so if you are experiencing issues, you may want to try the Github version of the package.

## Example
Once `lodi` is installed, you can load up R and type 

```
vignete("lodi")
```
to learn how to use the method.

## Bugs
If you encounter a bug, please open an issue in the 'Issues' tab on Github. 
## Contact
For questions ot feedback, please email Jonathan Boss at <bossjona@umich.edu> or Alexander Rix <alexrix@umich.edu>.
## References
Boss J, Mukherjee B, Ferguson KK, et al. Estimating outcome-exposure
associations when exposure biomarker detection limits vary across batches.
*Epidemiology*. 2019;30(5):746-755.
[10.1097/EDE.0000000000001052](https://doi.org/10.1097/EDE.0000000000001052)
