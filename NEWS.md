# lodi v0.9.2

## Minor changes
* Added check to see there is exposure data not code as `NA` that is below the given `lod` variable in `clmi`. This is good check in and of itself, but in particular it is supposed to help when users input an already log transformed exposure in clmi, while inputing an untransformed `lod` variable. This isn't perfect, but it hopefuly should help.
* Changed the optimization algorithm uses in the MLE procedure from `L-BFGS-B` to `BFGS`. Using the boxed constrained version of `BFGS` could cause non-finite values of the objective function due to large gradients when evaluating the initial value of the function. These large gradients could be due to the slack variables being activate initially, but it is unclear. At any rate, `clmi` now uses `BFGS` with a reparameterized variance (`:= exp(var)`) to ensure non-negativity, with the appropriate transformations to the fisher information matrix to ensure it gives the same answer as before.

# lodi v0.9.1

## Minor changes
* Added additional checks to `clmi` to protect against improper variable ordering.
* Added `verbose` option to `clmi` to print out parsing information.
* Updated `clmi` documentation to provide additional exposition on proper use of the package.

## Bug fixes
* Corrected calculation of censored likelihood maximum likelihood estimates, and named them.
* `clmi` now returns imputed data.frames that have the same row ordering as the original.

# lodi v0.9.0
Initial release.
