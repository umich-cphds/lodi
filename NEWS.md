# lodi v1.0.0

## Major changes
* lodi now requires R >= 3.5.0 due to object serialization changes
 
# lodi v0.9.2

## Minor changes
* Added badges to README.Rmd / README.md to show the CRAN version of lodi, Github version of lodi, and a Travis-CI badge showing whether or not lodi tests are passing. If the Github version and CRAN versions differ, the badges appear as different colors.
* Added some unit tests to be run when updating the package / before CRAN submission.
* Added check to ensure exposure data that is below the `lod` is coded as `NA`. Users might input an already log transformed exposure in clmi, while not transforming the `lod` variable. If the transformation is logarithmic, this will catch it
* Changed the optimization algorithm uses in the MLE procedure from `L-BFGS-B` to `BFGS` with reparameterized variance (ie  `exp(var)`) to ensure non-negativity instead of box constraints. Using `L-BFGS-B` could cause non-finite values of the objective function due to large gradients when evaluating the initial value of the function.

## Bug fixes
* Fixed a subtle bug when calling lodi inside of a function when using a custom exposure function. lodi was internally using the parent frame to enclose the transformation, which is the lodi package namespace, instead of the caller's frame.

# lodi v0.9.1

## Minor changes
* Added additional checks to `clmi` to protect against improper variable ordering.
* Added `verbose` option to `clmi` to print out parsing information.
* Updated `clmi` documentation to provide additional exposition on proper use of the package.

## Bug fixes
* Corrected calculation of censored likelihood maximum likelihood estimates, and named them.
* `clmi` now returns imputed data.frames that have the same row ordering as the original.

# lodi v0.9.0
Initial CRAN release.
